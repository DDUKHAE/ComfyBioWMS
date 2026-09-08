"""STAR (Spliced Transcripts Alignment to a Reference) node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binaries: STAR
Galaxy wrapper: galaxyproject/tools-iuc tools/rgrnastar/rg_rnaStar.xml
"""

import shlex
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from .common import discover_samples
except Exception:
    try:
        from nodes.class_2.common import discover_samples
    except Exception:
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, Path(rev_input).expanduser().resolve() if rev_input.strip() else None)]


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and str(custom_dir).strip():
        out = Path(custom_dir).expanduser().resolve()
    else:
        try:
            folder_paths = __import__("folder_paths")
            base = Path(folder_paths.get_output_directory())
        except Exception:
            base = Path.cwd() / "ComfyUI" / "output"
        out = base / node_name
    return out


def _dir(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"{label} is not a directory: {path}")
    return path


def _run(argv: list[str], cwd: Path) -> None:
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"STAR exited with code {e.returncode}: {shlex.join(argv)}") from e


class STARGenomeGenerate:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("star_index_dir",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "genome_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "gtf_file": ("STRING", {"default": ""}),
                "sjdb_overhang": ("INT", {"default": 100, "min": 1, "max": 1000}),
                "genome_sa_index_nbases": ("INT", {"default": 14, "min": 1, "max": 18}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        genome_fasta: str,
        index_dir: str = "",
        gtf_file: str = "",
        sjdb_overhang: int = 100,
        genome_sa_index_nbases: int = 14,
        threads: int = 4,
        extra_command: str = "",
    ):
        fa_path = _file(genome_fasta, "Genome FASTA")
        gtf_path = _file(gtf_file, "GTF Annotation") if gtf_file.strip() else None

        executable = shutil.which("STAR")
        if not executable:
            raise RuntimeError("STAR executable not found on PATH; install bioconda package star")

        out = _output_dir("STARGenomeGenerate", index_dir)
        out.mkdir(parents=True, exist_ok=True)

        argv = [
            executable,
            "--runMode", "genomeGenerate",
            "--genomeDir", str(out),
            "--genomeFastaFiles", str(fa_path),
            "--runThreadN", str(threads),
            "--genomeSAindexNbases", str(genome_sa_index_nbases),
        ]
        if gtf_path:
            argv.extend(["--sjdbGTFfile", str(gtf_path), "--sjdbOverhang", str(sjdb_overhang)])

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        sa_file = out / "SA"
        if not sa_file.is_file() or sa_file.stat().st_size == 0:
            raise RuntimeError(f"STAR failed to create valid index in {out}")

        return (str(out),)


class STARAlignReads:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = (
        "aligned_bam",
        "transcriptome_bam",
        "gene_counts_tab",
        "splice_junctions",
        "log_final",
        "unmapped_reads_fwd",
    )

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "star_index_dir": ("STRING", {"default": ""}),
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "quant_mode": (
                    ["Both (GeneCounts + TranscriptomeSAM)", "GeneCounts", "TranscriptomeSAM", "None"],
                    {"default": "Both (GeneCounts + TranscriptomeSAM)"},
                ),
                "twopass_mode": (["Basic", "None"], {"default": "Basic"}),
                "out_unmapped_fastx": ("BOOLEAN", {"default": True}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        star_index_dir: str,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        quant_mode: str = "Both (GeneCounts + TranscriptomeSAM)",
        twopass_mode: str = "Basic",
        out_unmapped_fastx: bool = True,
        extra_command: str = "",
    ):
        idx_dir = _dir(star_index_dir, "STAR Index Directory")
        samples = discover_samples(reads_fwd, reads_rev)

        executable = shutil.which("STAR")
        if not executable:
            raise RuntimeError("STAR executable not found on PATH; install bioconda package star")

        base_out = _output_dir("STARAlignReads", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(reads_fwd).strip()).is_file()
        aligned_bams, tx_bams, gene_counts, sjs, logs, unmapped = [], [], [], [], [], []

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)
            prefix = f"{sample_out}/"

            read_files = [str(fwd_path)]
            if rev_path:
                read_files.append(str(rev_path))

            argv = [
                executable,
                "--runMode", "alignReads",
                "--genomeDir", str(idx_dir),
                "--readFilesIn", *read_files,
                "--runThreadN", str(threads),
                "--outFileNamePrefix", prefix,
                "--outSAMtype", "BAM", "SortedByCoordinate",
            ]

            if str(fwd_path).endswith(".gz"):
                argv.extend(["--readFilesCommand", "zcat"])
            elif str(fwd_path).endswith(".bz2"):
                argv.extend(["--readFilesCommand", "bzcat"])

            quant_map = {
                "Both (GeneCounts + TranscriptomeSAM)": ["TranscriptomeSAM", "GeneCounts"],
                "GeneCounts": ["GeneCounts"],
                "TranscriptomeSAM": ["TranscriptomeSAM"],
            }
            if quant_mode in quant_map:
                argv.extend(["--quantMode", *quant_map[quant_mode]])

            if twopass_mode == "Basic":
                argv.extend(["--twopassMode", "Basic"])

            if out_unmapped_fastx:
                argv.extend(["--outSAMunmapped", "Within", "--outReadsUnmapped", "Fastx"])

            if extra_command.strip():
                argv.extend(shlex.split(extra_command))
            _run(argv, sample_out)

            aligned_bam = sample_out / "Aligned.sortedByCoord.out.bam"
            if not aligned_bam.is_file() or aligned_bam.stat().st_size == 0:
                raise RuntimeError(f"STAR failed to generate sorted BAM for {sample_id} in {sample_out}")

            t_bam = sample_out / "Aligned.toTranscriptome.out.bam"
            g_cnt = sample_out / "ReadsPerGene.out.tab"
            sj = sample_out / "SJ.out.tab"
            lg = sample_out / "Log.final.out"
            unm = sample_out / "Unmapped.out.mate1"

            aligned_bams.append(str(aligned_bam))
            tx_bams.append(str(t_bam) if t_bam.is_file() else "")
            gene_counts.append(str(g_cnt) if g_cnt.is_file() else "")
            sjs.append(str(sj) if sj.is_file() else "")
            logs.append(str(lg) if lg.is_file() else "")
            unmapped.append(str(unm) if unm.is_file() else "")

        if is_single:
            return (aligned_bams[0], tx_bams[0], gene_counts[0], sjs[0], logs[0], unmapped[0])
        return (
            ",".join(aligned_bams),
            ",".join(tx_bams),
            ",".join(gene_counts),
            ",".join(sjs),
            ",".join(logs),
            ",".join(unmapped),
        )


NODE_CLASS_MAPPINGS = {
    "STARGenomeGenerate": STARGenomeGenerate,
    "STARAlignReads": STARAlignReads,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "STARGenomeGenerate": "STAR: Generate Genome Index",
    "STARAlignReads": "STAR: Spliced Read Alignment",
}


# Backward compatibility aliases
STARGenomeGenerateNode = STARGenomeGenerate
STARAlignReadsNode = STARAlignReads

__all__ = [
    "STARGenomeGenerate",
    "STARGenomeGenerateNode",
    "STARAlignReads",
    "STARAlignReadsNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
