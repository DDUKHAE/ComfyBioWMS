"""BWA-MEM2 indexing and alignment nodes.

Python packages: none
External binary: bwa-mem2==2.3 (Conda: bioconda::bwa-mem2=2.3;
Apt: bwa-mem2 where provided by the distribution)
Galaxy wrappers: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/bwa_mem2/bwa-mem2-idx.xml and tools/bwa_mem2/bwa-mem2.xml
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



def _executable() -> str:
    executable = shutil.which("bwa-mem2")
    if not executable:
        raise RuntimeError(
            "bwa-mem2 executable not found on PATH; install conda package bwa-mem2=2.3"
        )
    return executable


def _run_logs(argv, cwd):
    try:
        subprocess.run(argv, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bwa-mem2 exited {e.returncode}: {shlex.join(argv)}") from e


def _run_alignment(argv, output_sam, cwd):
    try:
        with output_sam.open("wb") as output:
            subprocess.run(argv, cwd=cwd, stdout=output, check=True)
    except subprocess.CalledProcessError as e:
        output_sam.unlink(missing_ok=True)
        raise RuntimeError(f"bwa-mem2 exited {e.returncode}: {shlex.join(argv)}") from e


class BwaMem2Index:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("indexed_reference",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reference_fasta": ("STRING", {"default": ""}),
            },
            "optional": {
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, reference_fasta, output_dir: str = "", extra_command=""):
        source = _file(reference_fasta, "Reference FASTA")
        executable = _executable()
        out = _output_dir("BwaMem2Index", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        reference = out / source.name
        if source.resolve() != reference.resolve():
            shutil.copy2(source, reference)
        argv = [executable, "index", *shlex.split(extra_command), str(reference)]
        _run_logs(argv, out)
        for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
            sidecar = Path(str(reference) + suffix)
            if not sidecar.is_file() or sidecar.stat().st_size == 0:
                raise RuntimeError(f"bwa-mem2 did not create index sidecar: {sidecar}")
        return (str(reference),)


class BwaMem2Align:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("sam_file",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "indexed_reference": ("STRING", {"default": ""}),
                "read1": ("STRING", {"default": ""}),
            },
            "optional": {
                "read2": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "preset": (
                    ["illumina", "pacbio", "ont2d", "intractg"],
                    {"default": "illumina"},
                ),
                "minimum_seed_length": (
                    "INT",
                    {"default": 19, "min": 1, "max": 1000000},
                ),
                "band_width": ("INT", {"default": 100, "min": 1, "max": 1000000}),
                "minimum_score": ("INT", {"default": 30, "min": 0, "max": 1000000}),
                "read_group": ("STRING", {"default": ""}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        indexed_reference,
        read1,
        output_sam: str = "",
        read2="",
        threads=1,
        preset="illumina",
        minimum_seed_length=19,
        band_width=100,
        minimum_score=30,
        read_group="",
        extra_command="",
    ):
        reference = _file(indexed_reference, "Indexed reference")
        read1_path = _file(read1, "Read 1")
        read2_path = _file(read2, "Read 2") if read2 else None
        executable = _executable()
        output = Path(output_sam).expanduser().resolve() if (output_sam and str(output_sam).strip()) else (_output_dir("BwaMem2Align") / f"{Path(read1).stem}.sam")
        output.parent.mkdir(parents=True, exist_ok=True)
        
        argv = [executable, "mem", "-t", str(threads)]
        if preset != "illumina":
            argv += ["-x", preset]
        argv += [
            "-k",
            str(minimum_seed_length),
            "-w",
            str(band_width),
            "-T",
            str(minimum_score),
        ]
        if read_group:
            argv += ["-R", read_group]
        argv += kept + [str(reference), str(read1_path)]
        if read2_path:
            argv.append(str(read2_path))
        _run_alignment(argv, output, output.parent)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"bwa-mem2 did not create a nonempty SAM file: {output}")
        with output.open(encoding="utf-8") as handle:
            if not handle.readline().startswith("@"):
                raise RuntimeError(f"bwa-mem2 output is missing a SAM header: {output}")
        return (str(output),)


NODE_CLASS_MAPPINGS = {
    "BwaMem2Index": BwaMem2Index,
    "BwaMem2Align": BwaMem2Align,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "BwaMem2Index": "BWA-MEM2: Index Reference",
    "BwaMem2Align": "BWA-MEM2: Align Reads",
}


# Backward compatibility aliases
BwaMem2IndexNode = BwaMem2Index
BwaMem2AlignNode = BwaMem2Align

__all__ = [
    "BwaMem2Index",
    "BwaMem2IndexNode",
    "BwaMem2Align",
    "BwaMem2AlignNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
