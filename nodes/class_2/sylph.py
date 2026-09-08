"""Sylph (Ultra-fast and memory-efficient metagenomic profiler) node.

Python packages: none
External binaries: sylph
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
        raise RuntimeError(f"Sylph exited with code {e.returncode}: {shlex.join(argv)}") from e


class SylphProfile:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Metagenomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("sylph_profile_tsv",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
                "sylph_db": ("STRING", {"default": ""}),
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
        sylph_db: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None
        db_path = _file(sylph_db, "Sylph Database (.syldb)")

        executable = shutil.which("sylph")
        if not executable:
            raise RuntimeError("sylph executable not found on PATH; install bioconda package sylph")

        out = _output_dir("SylphProfile", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        out_tsv = out / f"{fwd_path.stem}_sylph_profile.tsv"

        argv = [
            executable,
            "profile",
            str(db_path),
            str(fwd_path),
            "-o",
            str(out_tsv),
            "-t",
            str(threads),
        ]
        if rev_path:
            argv.append(str(rev_path))

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
            raise RuntimeError(f"Sylph failed to produce profile TSV at {out_tsv}")

        return (str(out_tsv),)


NODE_CLASS_MAPPINGS = {"SylphProfile": SylphProfile}
NODE_DISPLAY_NAME_MAPPINGS = {"SylphProfile": "Sylph: Fast Metagenomic Profiler"}


# Backward compatibility aliases
SylphProfileNode = SylphProfile

__all__ = ["SylphProfile",
    "SylphProfileNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
