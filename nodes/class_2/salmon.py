"""Salmon (Accurate, fast, and bias-aware transcript expression quantification) node.

Supports single FASTQ/BAM inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binaries: salmon
Galaxy wrapper: galaxyproject/tools-iuc tools/salmon/salmon.xml
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
        import re
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, Path(rev_input).expanduser().resolve() if rev_input.strip() else None)]


try:
    from bioflow.runtime.command_runner import BioCommandRunner, resolve_tool_environment, validate_extra_command_conflicts
    from bioflow.runtime.artifacts import compute_input_fingerprint, get_run_output_dir
except Exception:
    try:
        import sys
        sys.path.append(str(Path(__file__).resolve().parents[2] / "engine" / "src"))
        from bioflow.runtime.command_runner import BioCommandRunner, resolve_tool_environment, validate_extra_command_conflicts
        from bioflow.runtime.artifacts import compute_input_fingerprint, get_run_output_dir
    except Exception:
        BioCommandRunner = None
        resolve_tool_environment = None
        compute_input_fingerprint = None


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


def _run(argv: list[str], cwd: Path, inputs=None, outputs=None, extra_command="", managed_flags=None) -> None:
    if BioCommandRunner is not None:
        BioCommandRunner.run(
            argv=argv,
            cwd=cwd,
            node_type="Salmon",
            inputs=inputs,
            outputs=outputs,
            extra_command=extra_command,
            managed_flags=managed_flags,
        )
    else:
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
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("salmon_index_dir",)

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        if compute_input_fingerprint is not None:
            return compute_input_fingerprint(**kwargs)
        return float("NaN")

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

        executable = None
        if resolve_tool_environment is not None:
            try:
                _, executable = resolve_tool_environment("salmon")
            except Exception:
                executable = shutil.which("salmon")
        else:
            executable = shutil.which("salmon")

        if not executable:
            raise RuntimeError("Salmon executable not found in any conda environment or PATH; install bioconda package salmon")

        out = _output_dir("SalmonIndex", index_dir)
        out.mkdir(parents=True, exist_ok=True)

        managed_flags = {"-t", "-i", "-k", "-p", "-d"}

        argv = [
            executable,
            "index",
            "-t", str(tx_path),
            "-i", str(out),
            "-k", str(kmer_len),
            "-p", str(threads),
        ]
        if decoy_path:
            argv += ["-d", str(decoy_path)]

        inputs = [tx_path]
        if decoy_path:
            inputs.append(decoy_path)

        try:
            _run(
                argv=argv,
                cwd=out,
                inputs=inputs,
                outputs=[out / "info.json"],
                extra_command=extra_command,
                managed_flags=managed_flags,
            )
        except TypeError:
            if extra_command.strip():
                argv.extend(shlex.split(extra_command))
            _run(argv, out)

        info_file = out / "info.json"
        if not info_file.is_file() or info_file.stat().st_size == 0:
            raise RuntimeError(f"Salmon failed to generate valid index in {out}")

        return (str(out),)


class SalmonQuantReads:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("quant_sf", "salmon_out_dir")

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        if compute_input_fingerprint is not None:
            return compute_input_fingerprint(**kwargs)
        return float("NaN")

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
        samples = discover_samples(reads_fwd, reads_rev)

        executable = None
        if resolve_tool_environment is not None:
            try:
                _, executable = resolve_tool_environment("salmon")
            except Exception:
                executable = shutil.which("salmon")
        else:
            executable = shutil.which("salmon")

        if not executable:
            raise RuntimeError("Salmon executable not found in any conda environment or PATH; install bioconda package salmon")

        base_out = _output_dir("SalmonQuantReads", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(reads_fwd.strip()).is_file()
        quant_files = []

        managed_flags = {
            "-i", "-l", "-o", "-p", "-1", "-2", "-r",
            "--validateMappings", "--gcBias", "--seqBias"
        }

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            lib_type = _map_strandedness_to_salmon(strandedness, paired=bool(rev_path))
            argv = [
                executable,
                "quant",
                "-i", str(idx_dir),
                "-l", lib_type,
                "-o", str(sample_out),
                "-p", str(threads),
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

            inputs = [fwd_path, idx_dir]
            if rev_path:
                inputs.append(rev_path)
            quant_file = sample_out / "quant.sf"
            outputs = [quant_file]

            try:
                _run(
                    argv=argv,
                    cwd=sample_out,
                    inputs=inputs,
                    outputs=outputs,
                    extra_command=extra_command,
                    managed_flags=managed_flags,
                )
            except TypeError:
                if extra_command.strip():
                    argv.extend(shlex.split(extra_command))
                _run(argv, sample_out)

            if not quant_file.is_file() or quant_file.stat().st_size == 0:
                raise RuntimeError(f"Salmon failed to produce quant.sf for {sample_id} in {sample_out}")

            quant_files.append(str(quant_file))

        if is_single:
            return quant_files[0], str(base_out)
        return ",".join(quant_files), str(base_out)


class SalmonQuantAlignment:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("quant_sf", "salmon_out_dir")

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        if compute_input_fingerprint is not None:
            return compute_input_fingerprint(**kwargs)
        return float("NaN")

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
        bams = []
        raw_bam = str(transcriptome_bam).strip()
        if "," in raw_bam:
            bams = [_file(b.strip(), "Transcriptome BAM") for b in raw_bam.split(",") if b.strip()]
        else:
            p = Path(raw_bam).expanduser().resolve()
            if p.is_file():
                bams = [p]
            elif p.is_dir():
                bams = sorted(list(p.glob("*.bam")) + list(p.glob("*/*.bam")))
                if not bams:
                    raise FileNotFoundError(f"No BAM files found in directory: {p}")
            else:
                raise FileNotFoundError(f"Transcriptome BAM is not a file or directory: {p}")

        executable = None
        if resolve_tool_environment is not None:
            try:
                _, executable = resolve_tool_environment("salmon")
            except Exception:
                executable = shutil.which("salmon")
        else:
            executable = shutil.which("salmon")

        if not executable:
            raise RuntimeError("Salmon executable not found in any conda environment or PATH; install bioconda package salmon")

        base_out = _output_dir("SalmonQuantAlignment", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(bams) == 1 and Path(raw_bam).is_file()
        quant_files = []
        used_sample_ids: set[str] = set()

        managed_flags = {"-t", "-l", "-a", "-o", "-p", "--sampleEq"}

        for idx, b in enumerate(bams):
            stem = b.stem
            for sfx in (".toTranscriptome.out", ".transcriptome", ".sorted", ".aligned", ".transcriptome.out"):
                if stem.endswith(sfx):
                    stem = stem[:-len(sfx)]
                    break
            sample_id = stem
            if sample_id in ("Aligned", "unaligned", "accepted_hits", "mapped") or not sample_id:
                parent_name = b.parent.name
                if parent_name and parent_name not in (".", "", "output", "STARAlignReads", "STAR", "SalmonQuantAlignment"):
                    sample_id = parent_name

            if sample_id in used_sample_ids:
                sample_id = f"{sample_id}_{idx + 1}"
            used_sample_ids.add(sample_id)

            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            lib_type = _map_strandedness_to_salmon(strandedness, paired=True)
            argv = [
                executable,
                "quant",
                "-t", str(tx_path),
                "-l", lib_type,
                "-a", str(b),
                "-o", str(sample_out),
                "-p", str(threads),
            ]
            if sample_eq:
                argv.append("--sampleEq")

            quant_file = sample_out / "quant.sf"
            inputs = [tx_path, b]
            outputs = [quant_file]

            try:
                _run(
                    argv=argv,
                    cwd=sample_out,
                    inputs=inputs,
                    outputs=outputs,
                    extra_command=extra_command,
                    managed_flags=managed_flags,
                )
            except TypeError:
                if extra_command.strip():
                    argv.extend(shlex.split(extra_command))
                _run(argv, sample_out)

            if not quant_file.is_file() or quant_file.stat().st_size == 0:
                raise RuntimeError(f"Salmon failed to produce quant.sf from BAM {b.name} in {sample_out}")

            quant_files.append(str(quant_file))

        if is_single:
            return quant_files[0], str(base_out)
        return ",".join(quant_files), str(base_out)

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
