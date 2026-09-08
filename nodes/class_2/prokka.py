"""Prokka rapid prokaryotic genome annotation node.

Python packages: none
External binaries: prokka
Galaxy wrapper: galaxyproject/tools-iuc tools/prokka/
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
        raise RuntimeError(f"Prokka exited with code {e.returncode}: {shlex.join(argv)}") from e


class Prokka:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Microbiome"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("annotation_gff", "protein_faa", "transcript_ffn")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "genome_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "kingdom": (["Bacteria", "Archaea", "Mitochondria", "Viruses"], {"default": "Bacteria"}),
                "genus": ("STRING", {"default": ""}),
                "species": ("STRING", {"default": ""}),
                "cpus": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        genome_fasta: str,
        output_dir: str = "",
        kingdom: str = "Bacteria",
        genus: str = "",
        species: str = "",
        cpus: int = 4,
        extra_command: str = "",
    ):
        fasta_path = _file(genome_fasta, "Genome FASTA")

        executable = shutil.which("prokka")
        if not executable:
            raise RuntimeError("prokka executable not found on PATH; install bioconda package prokka")

        out = _output_dir("Prokka", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        prefix = f"{fasta_path.stem}_prokka"

        argv = [
            executable,
            "--outdir",
            str(out),
            "--prefix",
            prefix,
            "--kingdom",
            kingdom,
            "--cpus",
            str(cpus),
            "--force",
        ]
        if genus.strip():
            argv += ["--genus", genus.strip()]
        if species.strip():
            argv += ["--species", species.strip()]

        argv.append(str(fasta_path))
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        gff_out = out / f"{prefix}.gff"
        faa_out = out / f"{prefix}.faa"
        ffn_out = out / f"{prefix}.ffn"

        if not gff_out.is_file() or gff_out.stat().st_size == 0:
            raise RuntimeError(f"Prokka failed to generate GFF annotation at {gff_out}")

        return (str(gff_out), str(faa_out) if faa_out.is_file() else "", str(ffn_out) if ffn_out.is_file() else "")


NODE_CLASS_MAPPINGS = {"Prokka": Prokka}
NODE_DISPLAY_NAME_MAPPINGS = {"Prokka": "Prokka: Rapid Prokaryotic Genome Annotation"}


# Backward compatibility aliases
ProkkaNode = Prokka

__all__ = ["Prokka",
    "ProkkaNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
