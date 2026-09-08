"""Bakta rapid & standardized annotation of bacterial genomes.

Python packages: none
External binaries: bakta
Galaxy wrapper: galaxyproject/tools-iuc tools/bakta/
"""

import shlex
import shutil
import subprocess
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


class Bakta:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Microbiome"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("annotation_gff3", "proteins_faa")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "genome_fasta": ("STRING", {"default": ""}),
                "bakta_db": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, genome_fasta: str, bakta_db: str, output_dir: str = "", threads: int = 4, extra_command: str = ""):
        fasta_path = _file(genome_fasta, "Genome FASTA")
        db_path = _dir(bakta_db, "Bakta Database")

        executable = shutil.which("bakta")
        if not executable:
            raise RuntimeError("bakta executable not found on PATH; install bioconda package bakta")

        out = _output_dir("Bakta", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        prefix = f"{fasta_path.stem}_bakta"

        argv = [
            executable,
            "--db",
            str(db_path),
            "--output",
            str(out),
            "--prefix",
            prefix,
            "--threads",
            str(threads),
            "--force",
            str(fasta_path),
        ]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Bakta failed ({proc.returncode}): {proc.stderr}")

        gff = out / f"{prefix}.gff3"
        faa = out / f"{prefix}.faa"

        if not gff.is_file() or gff.stat().st_size == 0:
            raise RuntimeError(f"Bakta failed to produce annotation at {gff}")

        return (str(gff), str(faa) if faa.is_file() else "")


NODE_CLASS_MAPPINGS = {"Bakta": Bakta}
NODE_DISPLAY_NAME_MAPPINGS = {"Bakta": "Bakta: Bacterial Genome Annotation"}


# Backward compatibility aliases
BaktaNode = Bakta

__all__ = ["Bakta",
    "BaktaNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
