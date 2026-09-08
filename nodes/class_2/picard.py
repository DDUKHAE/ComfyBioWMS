"""Picard MarkDuplicates node for identifying PCR and optical duplicates.

Python packages: none
External binaries: picard (or java -jar picard.jar)
Galaxy wrapper: galaxyproject/tools-iuc tools/picard/picard_MarkDuplicates.xml
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
        raise RuntimeError(f"Picard MarkDuplicates exited with code {e.returncode}: {shlex.join(argv)}") from e


class PicardMarkDuplicates:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("marked_bam", "metrics_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_bam": ("STRING", {"default": ""}),
            },
            "optional": {
                "remove_duplicates": ("BOOLEAN", {"default": False}),
                "optical_duplicate_pixel_distance": ("INT", {"default": 100, "min": 0, "max": 2500}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        input_bam: str,
        output_dir: str = "",
        remove_duplicates: bool = False,
        optical_duplicate_pixel_distance: int = 100,
        extra_command: str = "",
    ):
        bam_path = _file(input_bam, "Input BAM")

        executable = shutil.which("picard")
        if not executable:
            raise RuntimeError("picard executable not found on PATH; install bioconda package picard")

        out = _output_dir("PicardMarkDuplicates", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        marked_bam = out / f"markdup_{bam_path.name}"
        metrics_file = out / f"{bam_path.stem}_markdup_metrics.txt"

        argv = [
            executable,
            "MarkDuplicates",
            "-I",
            str(bam_path),
            "-O",
            str(marked_bam),
            "-M",
            str(metrics_file),
            f"--REMOVE_DUPLICATES={str(remove_duplicates).lower()}",
            f"--OPTICAL_DUPLICATE_PIXEL_DISTANCE={optical_duplicate_pixel_distance}",
        ]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not marked_bam.is_file() or marked_bam.stat().st_size == 0:
            raise RuntimeError(f"Picard MarkDuplicates failed to create BAM at {marked_bam}")

        return (str(marked_bam), str(metrics_file) if metrics_file.is_file() else "")


NODE_CLASS_MAPPINGS = {"PicardMarkDuplicates": PicardMarkDuplicates}
NODE_DISPLAY_NAME_MAPPINGS = {"PicardMarkDuplicates": "Picard: Mark Duplicates"}


# Backward compatibility aliases
PicardMarkDuplicatesNode = PicardMarkDuplicates

__all__ = ["PicardMarkDuplicates",
    "PicardMarkDuplicatesNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
