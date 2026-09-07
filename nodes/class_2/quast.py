"""QUAST genome assembly quality assessment node.

Python packages: none
External binary: quast.py==5.3.0 (Conda: bioconda::quast=5.3.0; Apt: quast)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/quast/quast.xml and tools/quast/macros.xml (QUAST 5.3.0+galaxy1)
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



def _run(argv, cwd):
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"quast.py exited with code {e.returncode}: {shlex.join(argv)}") from e


def _nonempty(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"quast.py did not create a nonempty {label}: {path}")
    return str(path)


def _validate_report_tsv(path: Path) -> str:
    _nonempty(path, "report.tsv")
    content = path.read_text(errors="replace")
    if not any("contig" in line.lower() or "length" in line.lower() for line in content.splitlines()):
        raise RuntimeError(f"report.tsv does not contain valid QUAST metrics: {path}")
    return str(path)


def _validate_report_html(path: Path) -> str:
    _nonempty(path, "report.html")
    content = path.read_text(errors="replace")
    if "<html" not in content.lower() and "<!doctype html" not in content.lower():
        raise RuntimeError(f"report.html is not a valid HTML document: {path}")
    return str(path)


class Quast:
    CATEGORY = "ComfyBIO/Genome Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("report_tsv", "report_html")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "assembly_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "reference": ("STRING", {"default": ""}),
                "features": ("STRING", {"default": ""}),
                "min_contig": ("INT", {"default": 500, "min": 0, "max": 1000000}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "large": ("BOOLEAN", {"default": False}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        assembly_fasta,
        output_dir: str = "",
        reference="",
        features="",
        min_contig=500,
        threads=4,
        large=False,
        extra_command="",
    ):
        assembly_path = _file(assembly_fasta, "Assembly FASTA")
        ref_path = _file(reference, "Reference") if reference else None
        feat_path = _file(features, "Features") if features else None

        executable = shutil.which("quast.py") or shutil.which("quast")
        if not executable:
            raise RuntimeError(
                "quast.py executable not found on PATH; install conda package quast=5.3.0"
            )

        out = _output_dir("Quast", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(
                f"[quast] ignored node-managed extra options: {' '.join(ignored)}",
                file=sys.stderr,
            )

        argv = [
            executable,
            str(assembly_path),
            "-o",
            str(out),
            "-m",
            str(min_contig),
            "-t",
            str(threads),
        ]
        if ref_path:
            argv += ["-r", str(ref_path)]
        if feat_path:
            argv += ["-g", str(feat_path)]
        if large:
            argv.append("--large")

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        report_tsv = out / "report.tsv"
        report_html = out / "report.html"

        return (
            _validate_report_tsv(report_tsv),
            _validate_report_html(report_html),
        )


NODE_CLASS_MAPPINGS = {"Quast": Quast}
NODE_DISPLAY_NAME_MAPPINGS = {"Quast": "QUAST: Evaluate Assembly"}


# Backward compatibility aliases
QuastNode = Quast

__all__ = ["Quast",
    "QuastNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
