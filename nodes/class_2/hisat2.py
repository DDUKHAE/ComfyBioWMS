"""HISAT2 (Hierarchical Indexing for Spliced Alignment of Transcripts 2) node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binaries: hisat2, hisat2-build
Galaxy wrapper: galaxyproject/tools-iuc tools/hisat2/hisat2.xml
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


def _run(argv: list[str], cwd: Path) -> None:
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"HISAT2 exited with code {e.returncode}: {shlex.join(argv)}") from e


class HISAT2Build:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("hisat2_index_prefix",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reference_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reference_fasta: str,
        index_dir: str = "",
        threads: int = 4,
        extra_command: str = "",
    ):
        ref_path = _file(reference_fasta, "Reference FASTA")
        executable = shutil.which("hisat2-build")
        if not executable:
            raise RuntimeError("hisat2-build executable not found on PATH; install bioconda package hisat2")

        out = _output_dir("HISAT2Build", index_dir)
        out.mkdir(parents=True, exist_ok=True)

        prefix = f"{out}/genome"
        argv = [
            executable,
            "-p", str(threads),
            str(ref_path),
            prefix,
        ]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        idx_files = list(out.glob("genome.*.ht2"))
        if not idx_files:
            raise RuntimeError(f"hisat2-build failed to produce index files in {out}")

        return (prefix,)


class HISAT2Align:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("aligned_sam", "summary_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "hisat2_index_prefix": ("STRING", {"default": ""}),
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "strandedness": ("STRING", {"default": "auto"}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        hisat2_index_prefix: str,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        strandedness: str = "auto",
        threads: int = 4,
        extra_command: str = "",
    ):
        samples = discover_samples(reads_fwd, reads_rev)

        executable = shutil.which("hisat2")
        if not executable:
            raise RuntimeError("hisat2 executable not found on PATH; install bioconda package hisat2")

        base_out = _output_dir("HISAT2Align", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(reads_fwd).strip()).is_file()
        sams, summaries = [], []

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            out_sam = sample_out / f"{sample_id}_hisat2.sam"
            summary_txt = sample_out / f"{sample_id}_hisat2_summary.txt"

            argv = [
                executable,
                "-x", hisat2_index_prefix,
                "-p", str(threads),
                "-S", str(out_sam),
                "--summary-file", str(summary_txt),
            ]
            if rev_path:
                argv += ["-1", str(fwd_path), "-2", str(rev_path)]
                s = strandedness.strip().lower()
                if "rev" in s or "rf" in s or "isr" in s:
                    argv += ["--rna-strandness", "RF"]
                elif "fwd" in s or "fr" in s or "isf" in s:
                    argv += ["--rna-strandness", "FR"]
            else:
                argv += ["-U", str(fwd_path)]
                s = strandedness.strip().lower()
                if "rev" in s or "r" in s:
                    argv += ["--rna-strandness", "R"]
                elif "fwd" in s or "f" in s:
                    argv += ["--rna-strandness", "F"]

            if extra_command.strip():
                argv.extend(shlex.split(extra_command))
            _run(argv, sample_out)

            if not out_sam.is_file() or out_sam.stat().st_size == 0:
                raise RuntimeError(f"hisat2 failed to produce SAM file for {sample_id} in {sample_out}")

            sams.append(str(out_sam))
            summaries.append(str(summary_txt) if summary_txt.is_file() else "")

        if is_single:
            return sams[0], summaries[0]
        return ",".join(sams), ",".join(summaries)


NODE_CLASS_MAPPINGS = {
    "HISAT2Build": HISAT2Build,
    "HISAT2Align": HISAT2Align,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "HISAT2Build": "HISAT2: Build Reference Index",
    "HISAT2Align": "HISAT2: Spliced Read Alignment",
}


# Backward compatibility aliases
HISAT2BuildNode = HISAT2Build
HISAT2AlignNode = HISAT2Align

__all__ = [
    "HISAT2Build",
    "HISAT2BuildNode",
    "HISAT2Align",
    "HISAT2AlignNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
