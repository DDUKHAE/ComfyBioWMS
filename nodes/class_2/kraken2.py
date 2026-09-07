"""Kraken2 taxonomic classification node.

Python packages: none
External binaries: kraken2
Galaxy wrapper: galaxyproject/tools-iuc tools/kraken2/kraken2.xml
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
        raise RuntimeError(f"Kraken2 exited with code {e.returncode}: {shlex.join(argv)}") from e


class Kraken2Classify:
    CATEGORY = "ComfyBIO/Metagenomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("kraken_report_txt", "classified_output_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
                "kraken2_db": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "confidence": ("FLOAT", {"default": 0.0, "min": 0.0, "max": 1.0, "step": 0.05}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        kraken2_db: str,
        output_dir: str = "",
        reads_rev: str = "",
        confidence: float = 0.0,
        threads: int = 4,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None
        db_path = _dir(kraken2_db, "Kraken2 Database Directory")

        executable = shutil.which("kraken2")
        if not executable:
            raise RuntimeError("kraken2 executable not found on PATH; install bioconda package kraken2")

        out = _output_dir("Kraken2Classify", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[Kraken2] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        report_txt = out / f"{fwd_path.stem}_kraken2_report.txt"
        classified_txt = out / f"{fwd_path.stem}_kraken2_output.txt"

        argv = [
            executable,
            "--db",
            str(db_path),
            "--report",
            str(report_txt),
            "--output",
            str(classified_txt),
            "--threads",
            str(threads),
            "--confidence",
            str(confidence),
        ]
        if fwd_path.name.endswith(".gz"):
            argv.append("--gzip-compressed")
        elif fwd_path.name.endswith(".bz2"):
            argv.append("--bzip2-compressed")

        if rev_path:
            argv.append("--paired")
            argv += [str(fwd_path), str(rev_path)]
        else:
            argv.append(str(fwd_path))

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not report_txt.is_file() or report_txt.stat().st_size == 0:
            raise RuntimeError(f"Kraken2 failed to produce report file at {report_txt}")

        return (str(report_txt), str(classified_txt) if classified_txt.is_file() else "")


NODE_CLASS_MAPPINGS = {"Kraken2Classify": Kraken2Classify}
NODE_DISPLAY_NAME_MAPPINGS = {"Kraken2Classify": "Kraken2: Taxonomic Classification"}


# Backward compatibility aliases
Kraken2ClassifyNode = Kraken2Classify

__all__ = ["Kraken2Classify",
    "Kraken2ClassifyNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
