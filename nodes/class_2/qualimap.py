"""Qualimap RNA-seq quality control node.

Python packages: none
External binaries: qualimap
Galaxy wrapper: galaxyproject/tools-iuc tools/qualimap/
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
        raise RuntimeError(f"Qualimap exited with code {e.returncode}: {shlex.join(argv)}") from e


class QualimapRNASeq:
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("report_html", "results_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
                "gtf_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "strandedness": ("STRING", {"default": "non-strand-specific"}),
                "java_mem_size": ("STRING", {"default": "4G"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        gtf_file: str,
        output_dir: str = "",
        strandedness: str = "non-strand-specific",
        java_mem_size: str = "4G",
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM")
        gtf_path = _file(gtf_file, "Annotation GTF")

        executable = shutil.which("qualimap")
        if not executable:
            raise RuntimeError(
                "qualimap executable not found on PATH; install bioconda package qualimap"
            )

        out = _output_dir("QualimapRNASeq", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[Qualimap] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        strand_protocol = "non-strand-specific"
        s = strandedness.strip().lower()
        if "rev" in s or "rf" in s or "isr" in s:
            strand_protocol = "strand-specific-reverse"
        elif "fwd" in s or "fr" in s or "isf" in s:
            strand_protocol = "strand-specific-forward"

        argv = [
            executable,
            "rnaseq",
            "-bam",
            str(bam_path),
            "-gtf",
            str(gtf_path),
            "-outdir",
            str(out),
            "-outfile",
            "qualimapReport.html",
            "-outformat",
            "html",
            "-p",
            strand_protocol,
            f"--java-mem-size={java_mem_size}",
        ]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        html_file = out / "qualimapReport.html"
        results_txt = out / "rnaseq_qc_results.txt"

        if not html_file.is_file() or html_file.stat().st_size == 0:
            raise RuntimeError(f"Qualimap failed to create HTML report at {html_file}")

        return (str(html_file), str(results_txt) if results_txt.is_file() else "")


NODE_CLASS_MAPPINGS = {"QualimapRNASeq": QualimapRNASeq}
NODE_DISPLAY_NAME_MAPPINGS = {"QualimapRNASeq": "Qualimap: RNA-seq Quality Control"}


# Backward compatibility aliases
QualimapRNASeqNode = QualimapRNASeq

__all__ = ["QualimapRNASeq",
    "QualimapRNASeqNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
