"""deepTools computeMatrix and plotProfile node.

Python packages: deeptools
External binaries: computeMatrix, plotProfile
Galaxy wrapper: galaxyproject/tools-iuc tools/deeptools/
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



class DeeptoolsMatrix:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Epigenomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("matrix_gz", "profile_plot_png")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bigwig_file": ("STRING", {"default": ""}),
                "regions_bed": ("STRING", {"default": ""}),
            },
            "optional": {
                "reference_point": (["TSS", "TES", "center"], {"default": "TSS"}),
                "upstream": ("INT", {"default": 2000, "min": 0, "max": 50000}),
                "downstream": ("INT", {"default": 2000, "min": 0, "max": 50000}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bigwig_file: str,
        regions_bed: str,
        output_dir: str = "",
        reference_point: str = "TSS",
        upstream: int = 2000,
        downstream: int = 2000,
        threads: int = 4,
        extra_command: str = "",
    ):
        bw_path = _file(bigwig_file, "BigWig file")
        bed_path = _file(regions_bed, "Regions BED file")

        executable = shutil.which("computeMatrix")
        if not executable:
            raise RuntimeError("computeMatrix executable not found on PATH; install bioconda package deeptools")

        plot_exec = shutil.which("plotProfile")

        out = _output_dir("DeeptoolsMatrix", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        matrix_gz = out / f"{bw_path.stem}_matrix.gz"
        plot_png = out / f"{bw_path.stem}_profile.png"

        argv = [
            executable,
            "reference-point",
            "--referencePoint",
            reference_point,
            "-b",
            str(upstream),
            "-a",
            str(downstream),
            "-R",
            str(bed_path),
            "-S",
            str(bw_path),
            "-o",
            str(matrix_gz),
            "-p",
            str(threads),
        ]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"computeMatrix failed ({proc.returncode}): {proc.stderr}")

        if plot_exec and matrix_gz.is_file():
            subprocess.run([plot_exec, "-m", str(matrix_gz), "-out", str(plot_png)], capture_output=True)

        return (str(matrix_gz), str(plot_png) if plot_png.is_file() else "")


NODE_CLASS_MAPPINGS = {"DeeptoolsMatrix": DeeptoolsMatrix}
NODE_DISPLAY_NAME_MAPPINGS = {"DeeptoolsMatrix": "deepTools: Compute Matrix & Profile"}


# Backward compatibility aliases
DeeptoolsMatrixNode = DeeptoolsMatrix

__all__ = ["DeeptoolsMatrix",
    "DeeptoolsMatrixNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
