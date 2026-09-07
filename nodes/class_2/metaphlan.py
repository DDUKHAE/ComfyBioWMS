"""MetaPhlAn taxonomic profiling node.

Python packages: none
External binaries: metaphlan
Galaxy wrapper: galaxyproject/tools-iuc tools/metaphlan/
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
        raise RuntimeError(f"MetaPhlAn exited with code {e.returncode}: {shlex.join(argv)}") from e


class MetaPhlAn:
    CATEGORY = "ComfyBIO/Microbiome"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("profiled_metagenome_tsv",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "input_type": (["fastq", "fasta", "bowtie2out", "sam"], {"default": "fastq"}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        input_type: str = "fastq",
        threads: int = 4,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads")
        rev_path = _file(reads_rev, "Reverse reads") if reads_rev.strip() else None

        executable = shutil.which("metaphlan")
        if not executable:
            raise RuntimeError("metaphlan executable not found on PATH; install bioconda package metaphlan")

        out = _output_dir("MetaPhlAn", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        out_tsv = out / f"{fwd_path.stem}_profile.tsv"

        reads_arg = f"{fwd_path},{rev_path}" if rev_path else str(fwd_path)

        argv = [
            executable,
            reads_arg,
            "--input_type",
            input_type,
            "--output_file",
            str(out_tsv),
            "--nproc",
            str(threads),
        ]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
            raise RuntimeError(f"MetaPhlAn failed to produce profile TSV at {out_tsv}")

        return (str(out_tsv),)


NODE_CLASS_MAPPINGS = {"MetaPhlAn": MetaPhlAn}
NODE_DISPLAY_NAME_MAPPINGS = {"MetaPhlAn": "MetaPhlAn: Taxonomic Profiler"}


# Backward compatibility aliases
MetaPhlAnNode = MetaPhlAn

__all__ = ["MetaPhlAn",
    "MetaPhlAnNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
