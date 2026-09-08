"""Bracken (Bayesian Re-estimation of Abundance with Kraken) node.

Python packages: none
External binaries: bracken
Galaxy wrapper: galaxyproject/tools-iuc tools/bracken/bracken.xml
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
        raise RuntimeError(f"Bracken exited with code {e.returncode}: {shlex.join(argv)}") from e


class Bracken:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Metagenomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("bracken_abundance_tsv", "bracken_report_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "kraken_report": ("STRING", {"default": ""}),
                "bracken_db": ("STRING", {"default": ""}),
            },
            "optional": {
                "taxonomic_level": (["S (Species)", "G (Genus)", "F (Family)", "O (Order)", "C (Class)", "P (Phylum)", "D (Domain)"], {"default": "S (Species)"}),
                "read_length": ("INT", {"default": 100, "min": 25, "max": 1000}),
                "threshold": ("INT", {"default": 10, "min": 0, "max": 100000}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        kraken_report: str,
        bracken_db: str,
        output_dir: str = "",
        taxonomic_level: str = "S (Species)",
        read_length: int = 100,
        threshold: int = 10,
        extra_command: str = "",
    ):
        report_path = _file(kraken_report, "Kraken Report")
        db_path = _dir(bracken_db, "Bracken Database Directory")

        executable = shutil.which("bracken")
        if not executable:
            raise RuntimeError("bracken executable not found on PATH; install bioconda package bracken")

        out = _output_dir("Bracken", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        level = taxonomic_level[0]  # S, G, F, etc.
        out_abundance = out / f"{report_path.stem}_bracken_abund_{level}.tsv"
        out_report = out / f"{report_path.stem}_bracken_report_{level}.txt"

        argv = [
            executable,
            "-d",
            str(db_path),
            "-i",
            str(report_path),
            "-o",
            str(out_abundance),
            "-w",
            str(out_report),
            "-r",
            str(read_length),
            "-l",
            level,
            "-t",
            str(threshold),
        ]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not out_abundance.is_file() or out_abundance.stat().st_size == 0:
            raise RuntimeError(f"Bracken failed to produce abundance file at {out_abundance}")

        return (str(out_abundance), str(out_report) if out_report.is_file() else "")


NODE_CLASS_MAPPINGS = {"Bracken": Bracken}
NODE_DISPLAY_NAME_MAPPINGS = {"Bracken": "Bracken: Re-estimate Abundance"}


# Backward compatibility aliases
BrackenNode = Bracken

__all__ = ["Bracken",
    "BrackenNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
