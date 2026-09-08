"""bwa-mem2 alignment node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binary: bwa-mem2==2.3 (Conda: bioconda::bwa-mem2=2.3; Apt: bwa-mem2)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/bwa_mem2/bwa_mem2.xml (bwa-mem2 2.2.1+galaxy1; executable dependency 2.3)
"""

import shlex
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from .common import discover_samples
except Exception:
    try:
        from nodes.class_2.common import discover_samples
    except Exception:
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, Path(rev_input).expanduser().resolve() if rev_input.strip() else None)]


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


def _executable():
    executable = shutil.which("bwa-mem2")
    if not executable:
        raise RuntimeError(
            "bwa-mem2 executable not found on PATH; install conda package bwa-mem2=2.3"
        )
    return executable


def _run_logs(argv, cwd):
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"bwa-mem2 exited with code {e.returncode}: {shlex.join(argv)}"
        ) from e


def _run_alignment(argv, output_file: Path, cwd: Path):
    try:
        with output_file.open("wb") as handle:
            subprocess.run(argv, cwd=str(cwd), stdout=handle, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"bwa-mem2 exited with code {e.returncode}: {shlex.join(argv)}"
        ) from e


class BwaMem2Index:
    OUTPUT_NODE = True
    OUPUT_NODE = True
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
    OUTPUT_NODE = True
    OUPUT_NODE = True
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
        samples = discover_samples(read1, read2)
        executable = _executable()

        base_out = _output_dir("BwaMem2Align")
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(read1).strip()).is_file()
        sam_files = []

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)
            output = sample_out / (f"{sample_id}.sam" if not is_single else (Path(output_sam).name if output_sam.strip() else f"{fwd_path.stem}.sam"))

            argv = [executable, "mem", "-t", str(threads)]
            if preset != "illumina":
                argv += ["-x", preset]
            argv += [
                "-k", str(minimum_seed_length),
                "-w", str(band_width),
                "-T", str(minimum_score),
            ]
            if read_group:
                argv += ["-R", read_group]
            if extra_command.strip():
                argv.extend(shlex.split(extra_command))

            argv += [str(reference), str(fwd_path)]
            if rev_path:
                argv.append(str(rev_path))

            _run_alignment(argv, output, sample_out)
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f"bwa-mem2 did not create a nonempty SAM file for {sample_id}: {output}")
            with output.open(encoding="utf-8") as handle:
                if not handle.readline().startswith("@"):
                    raise RuntimeError(f"bwa-mem2 output is missing a SAM header: {output}")

            sam_files.append(str(output))

        if is_single:
            return (sam_files[0],)
        return (",".join(sam_files),)


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
