"""StringTie transcript assembly and quantification node.

Python packages: none
External binaries: stringtie
Galaxy wrapper: galaxyproject/tools-iuc tools/stringtie/stringtie.xml
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
        raise RuntimeError(f"StringTie exited with code {e.returncode}: {shlex.join(argv)}") from e


class StringTie:
    CATEGORY = "ComfyBIO/Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("assembled_gtf", "gene_abundances_tsv")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "guide_gtf": ("STRING", {"default": ""}),
                "strandedness": ("STRING", {"default": "auto"}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        output_dir: str = "",
        guide_gtf: str = "",
        strandedness: str = "auto",
        threads: int = 4,
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM")
        gtf_path = _file(guide_gtf, "Guide GTF") if guide_gtf.strip() else None

        executable = shutil.which("stringtie")
        if not executable:
            raise RuntimeError(
                "stringtie executable not found on PATH; install bioconda package stringtie"
            )

        out = _output_dir("StringTie", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[StringTie] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        out_gtf = out / f"{bam_path.stem}_transcripts.gtf"
        out_abund = out / f"{bam_path.stem}_gene_abundances.tsv"

        argv = [
            executable,
            str(bam_path),
            "-o",
            str(out_gtf),
            "-A",
            str(out_abund),
            "-p",
            str(threads),
        ]
        if gtf_path:
            argv += ["-G", str(gtf_path)]

        s = strandedness.strip().lower()
        if s in ("reverse", "rf", "isr"):
            argv.append("--rf")
        elif s in ("forward", "fr", "isf"):
            argv.append("--fr")

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not out_gtf.is_file() or out_gtf.stat().st_size == 0:
            raise RuntimeError(f"StringTie failed to create assembled GTF at {out_gtf}")

        return (str(out_gtf), str(out_abund) if out_abund.is_file() else "")


NODE_CLASS_MAPPINGS = {"StringTie": StringTie}
NODE_DISPLAY_NAME_MAPPINGS = {"StringTie": "StringTie: Transcript Assembly & Quant"}


# Backward compatibility aliases
StringTieNode = StringTie

__all__ = ["StringTie",
    "StringTieNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
