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
    OUTPUT_NODE = True
    OUPUT_NODE = True
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
        raw = str(bam_file).strip()
        p = Path(raw).expanduser().resolve()
        if p.is_dir():
            bams = sorted(p.rglob("*.bam"))
            if not bams:
                raise FileNotFoundError(f"No BAM files found in directory: {p}")
        elif "," in raw:
            bams = [Path(x.strip()).expanduser().resolve() for x in raw.split(",") if x.strip()]
        else:
            bams = [_file(raw, "Input BAM")]

        gtf_path = _file(guide_gtf, "Guide GTF") if guide_gtf.strip() else None

        executable = shutil.which("stringtie")
        if not executable:
            raise RuntimeError(
                "stringtie executable not found on PATH; install bioconda package stringtie"
            )

        base_out = _output_dir("StringTie", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(bams) == 1 and p.is_file()
        gtfs, abunds = [], []

        for bam_path in bams:
            sample_out = base_out if is_single else (base_out / bam_path.stem)
            sample_out.mkdir(parents=True, exist_ok=True)

            out_gtf = sample_out / f"{bam_path.stem}_transcripts.gtf"
            out_abund = sample_out / f"{bam_path.stem}_gene_abundances.tsv"

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
            _run(argv, sample_out)

            if not out_gtf.is_file() or out_gtf.stat().st_size == 0:
                raise RuntimeError(f"StringTie failed to create assembled GTF at {out_gtf}")

            gtfs.append(str(out_gtf))
            abunds.append(str(out_abund) if out_abund.is_file() else "")

        if is_single:
            return (gtfs[0], abunds[0])
        return (",".join(gtfs), ",".join(abunds))


NODE_CLASS_MAPPINGS = {"StringTie": StringTie}
NODE_DISPLAY_NAME_MAPPINGS = {"StringTie": "StringTie: Transcript Assembly & Quant"}


# Backward compatibility aliases
StringTieNode = StringTie

__all__ = ["StringTie",
    "StringTieNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
