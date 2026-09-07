"""BBTools (BBSplit: Split sequencing reads across multiple references to remove contamination) node.

Python packages: none
External binaries: bbsplit.sh (from BBMap / BBTools)
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



def _run(argv: list[str], cwd: Path) -> None:
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"BBSplit exited with code {e.returncode}: {shlex.join(argv)}") from e


class BBSplit:
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("clean_reads_fwd", "clean_reads_rev", "refstats_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
                "contaminant_reference_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "max_memory_gb": ("INT", {"default": 8, "min": 1, "max": 512}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        contaminant_reference_fasta: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        max_memory_gb: int = 8,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None
        ref_path = _file(contaminant_reference_fasta, "Contaminant Reference FASTA")

        executable = shutil.which("bbsplit.sh")
        if not executable:
            raise RuntimeError(
                "bbsplit.sh executable not found on PATH; install bioconda package bbmap"
            )

        out = _output_dir("BBSplit", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[BBSplit] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        clean_fwd = out / f"clean_{fwd_path.name}"
        clean_rev = out / f"clean_{rev_path.name}" if rev_path else None
        stats_file = out / "refstats.txt"

        argv = [
            executable,
            f"-Xmx{max_memory_gb}g",
            f"threads={threads}",
            f"ref={ref_path}",
            f"refstats={stats_file}",
        ]
        if rev_path:
            argv += [
                f"in1={fwd_path}",
                f"in2={rev_path}",
                f"out_u1={clean_fwd}",
                f"out_u2={clean_rev}",
            ]
        else:
            argv += [
                f"in={fwd_path}",
                f"out_u={clean_fwd}",
            ]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not clean_fwd.is_file() or clean_fwd.stat().st_size == 0:
            raise RuntimeError(f"BBSplit failed to produce clean unmapped reads at {clean_fwd}")

        return (str(clean_fwd), str(clean_rev) if clean_rev and clean_rev.is_file() else "", str(stats_file))


NODE_CLASS_MAPPINGS = {"BBSplit": BBSplit}
NODE_DISPLAY_NAME_MAPPINGS = {"BBSplit": "BBTools: BBSplit Host/Contaminant Filter"}


# Backward compatibility aliases
BBSplitNode = BBSplit

__all__ = ["BBSplit",
    "BBSplitNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
