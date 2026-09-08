"""SortMeRNA ribosomal RNA filtering node.

Python packages: none
External binaries: sortmerna
Galaxy wrapper: galaxyproject/tools-iuc tools/sortmerna/sortmerna.xml
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
        raise RuntimeError(f"SortMeRNA exited with code {e.returncode}: {shlex.join(argv)}") from e


class SortMeRNA:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("non_rrna_reads_fwd", "non_rrna_reads_rev", "sortmerna_log")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
                "rrna_database_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        rrna_database_fasta: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None
        db_path = _file(rrna_database_fasta, "rRNA Database FASTA")

        executable = shutil.which("sortmerna")
        if not executable:
            raise RuntimeError(
                "sortmerna executable not found on PATH; install bioconda package sortmerna"
            )

        out = _output_dir("SortMeRNA", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        aligned_prefix = str(out / "rrna_aligned")
        other_prefix = str(out / "non_rrna")

        argv = [
            executable,
            "--ref",
            str(db_path),
            "--reads",
            str(fwd_path),
            "--workdir",
            str(out),
            "--aligned",
            aligned_prefix,
            "--other",
            other_prefix,
            "--fastx",
            "--threads",
            str(threads),
        ]
        if rev_path:
            argv += ["--reads", str(rev_path), "--paired_out"]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        # Output resolution
        clean_fwd = list(out.glob("non_rrna*.f*q*"))
        log_files = list(out.glob("*log*")) + list(out.glob("*.log"))
        log_file = str(log_files[0]) if log_files else ""

        if not clean_fwd:
            raise RuntimeError(f"SortMeRNA failed to produce non-rRNA output reads in {out}")

        fwd_file = str(clean_fwd[0])
        rev_file = str(clean_fwd[1]) if len(clean_fwd) > 1 else ""

        return (fwd_file, rev_file, log_file)


NODE_CLASS_MAPPINGS = {"SortMeRNA": SortMeRNA}
NODE_DISPLAY_NAME_MAPPINGS = {"SortMeRNA": "SortMeRNA: Ribosomal RNA Filter"}


# Backward compatibility aliases
SortMeRNANode = SortMeRNA

__all__ = ["SortMeRNA",
    "SortMeRNANode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
