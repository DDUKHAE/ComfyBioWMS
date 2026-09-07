"""Salmon (Accurate, fast, and bias-aware transcript expression quantification) node.

Python packages: none
External binaries: salmon
Galaxy wrapper: galaxyproject/tools-iuc tools/salmon/salmon.xml
"""

import shlex
import shutil
import subprocess
import sys
from pathlib import Path


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
        raise RuntimeError(f"Salmon exited with code {e.returncode}: {shlex.join(argv)}") from e


def _map_strandedness_to_salmon(strandedness: str, paired: bool) -> str:
    s = strandedness.strip().lower()
    if s in ("a", "auto", "automatic"):
        return "A"
    if paired:
        if s in ("reverse", "stranded_reverse", "rf", "isr"):
            return "ISR"
        if s in ("forward", "stranded_forward", "fr", "isf"):
            return "ISF"
        return "IU"
    else:
        if s in ("reverse", "stranded_reverse", "sr"):
            return "SR"
        if s in ("forward", "stranded_forward", "sf"):
            return "SF"
        return "U"


class SalmonIndex:
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("salmon_index_dir",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "transcripts_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "decoys_file": ("STRING", {"default": ""}),
                "kmer_len": ("INT", {"default": 31, "min": 11, "max": 31}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        transcripts_fasta: str,
        index_dir: str = "",
        decoys_file: str = "",
        kmer_len: int = 31,
        threads: int = 4,
        extra_command: str = "",
    ):
        tx_path = _file(transcripts_fasta, "Transcripts FASTA")
        decoy_path = _file(decoys_file, "Decoy File") if decoys_file.strip() else None

        executable = shutil.which("salmon")
        if not executable:
            raise RuntimeError("Salmon executable not found on PATH; install bioconda package salmon")

        out = _output_dir("SalmonIndex", index_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[Salmon] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        argv = [
            executable,
            "index",
            "-t",
            str(tx_path),
            "-i",
            str(out),
            "-k",
            str(kmer_len),
            "-p",
            str(threads),
        ]
        if decoy_path:
            argv += ["-d", str(decoy_path)]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        info_file = out / "info.json"
        if not info_file.is_file() or info_file.stat().st_size == 0:
            raise RuntimeError(f"Salmon failed to generate valid index in {out}")

        return (str(out),)


class SalmonQuantReads:
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("quant_sf", "salmon_out_dir")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "salmon_index_dir": ("STRING", {"default": ""}),
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "strandedness": ("STRING", {"default": "auto"}),
                "validate_mappings": ("BOOLEAN", {"default": True}),
                "gc_bias": ("BOOLEAN", {"default": True}),
                "seq_bias": ("BOOLEAN", {"default": True}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        salmon_index_dir: str,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        strandedness: str = "auto",
        validate_mappings: bool = True,
        gc_bias: bool = True,
        seq_bias: bool = True,
        threads: int = 4,
        extra_command: str = "",
    ):
        idx_dir = _dir(salmon_index_dir, "Salmon Index Directory")
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None

        executable = shutil.which("salmon")
        if not executable:
            raise RuntimeError("Salmon executable not found on PATH; install bioconda package salmon")

        out = _output_dir("SalmonQuantReads", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[Salmon] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        lib_type = _map_strandedness_to_salmon(strandedness, paired=bool(rev_path))

        argv = [
            executable,
            "quant",
            "-i",
            str(idx_dir),
            "-l",
            lib_type,
            "-o",
            str(out),
            "-p",
            str(threads),
        ]
        if rev_path:
            argv += ["-1", str(fwd_path), "-2", str(rev_path)]
        else:
            argv += ["-r", str(fwd_path)]

        if validate_mappings:
            argv.append("--validateMappings")
        if gc_bias:
            argv.append("--gcBias")
        if seq_bias:
            argv.append("--seqBias")

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        quant_file = out / "quant.sf"
        if not quant_file.is_file() or quant_file.stat().st_size == 0:
            raise RuntimeError(f"Salmon failed to produce quant.sf in {out}")

        return (str(quant_file), str(out))


class SalmonQuantAlignment:
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("quant_sf", "salmon_out_dir")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "transcripts_fasta": ("STRING", {"default": ""}),
                "transcriptome_bam": ("STRING", {"default": ""}),
            },
            "optional": {
                "strandedness": ("STRING", {"default": "auto"}),
                "sample_eq": ("BOOLEAN", {"default": False}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        transcripts_fasta: str,
        transcriptome_bam: str,
        output_dir: str = "",
        strandedness: str = "auto",
        sample_eq: bool = False,
        threads: int = 4,
        extra_command: str = "",
    ):
        tx_path = _file(transcripts_fasta, "Transcripts FASTA")
        bam_path = _file(transcriptome_bam, "Transcriptome BAM")

        executable = shutil.which("salmon")
        if not executable:
            raise RuntimeError("Salmon executable not found on PATH; install bioconda package salmon")

        out = _output_dir("SalmonQuantAlignment", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[Salmon] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        lib_type = _map_strandedness_to_salmon(strandedness, paired=True)

        argv = [
            executable,
            "quant",
            "-t",
            str(tx_path),
            "-l",
            lib_type,
            "-a",
            str(bam_path),
            "-o",
            str(out),
            "-p",
            str(threads),
        ]
        if sample_eq:
            argv.append("--sampleEq")

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        quant_file = out / "quant.sf"
        if not quant_file.is_file() or quant_file.stat().st_size == 0:
            raise RuntimeError(f"Salmon failed to produce quant.sf from BAM in {out}")

        return (str(quant_file), str(out))


NODE_CLASS_MAPPINGS = {
    "SalmonIndex": SalmonIndex,
    "SalmonQuantReads": SalmonQuantReads,
    "SalmonQuantAlignment": SalmonQuantAlignment,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "SalmonIndex": "Salmon: Transcriptome Index",
    "SalmonQuantReads": "Salmon: FastQ Quasi-mapping Quant",
    "SalmonQuantAlignment": "Salmon: Transcriptome BAM Quant",
}


# Backward compatibility aliases
SalmonIndexNode = SalmonIndex
SalmonQuantReadsNode = SalmonQuantReads
SalmonQuantAlignmentNode = SalmonQuantAlignment

__all__ = [
    "SalmonIndex",
    "SalmonIndexNode",
    "SalmonQuantReads",
    "SalmonQuantReadsNode",
    "SalmonQuantAlignment",
    "SalmonQuantAlignmentNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
