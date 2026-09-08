"""Preseq library complexity prediction node.

Python packages: none
External binaries: preseq
Galaxy wrapper: galaxyproject/tools-iuc tools/preseq/
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



class Preseq:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("preseq_estimates_txt",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "paired_end": ("BOOLEAN", {"default": True}),
                "extrapolate_yield": ("INT", {"default": 100000000, "min": 1000000, "max": 1000000000}),
                "step_size": ("INT", {"default": 1000000, "min": 10000, "max": 50000000}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        output_dir: str = "",
        paired_end: bool = True,
        extrapolate_yield: int = 100000000,
        step_size: int = 1000000,
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM")

        executable = shutil.which("preseq")
        if not executable:
            raise RuntimeError("preseq executable not found on PATH; install bioconda package preseq")

        out = _output_dir("Preseq", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        out_txt = out / f"{bam_path.stem}_preseq_lc_extrap.txt"

        argv = [
            executable,
            "lc_extrap",
            "-bam",
            "-o",
            str(out_txt),
            "-e",
            str(extrapolate_yield),
            "-s",
            str(step_size),
        ]
        if paired_end:
            argv.append("-pe")
        argv.append(str(bam_path))

        proc = subprocess.run(argv + kept, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"preseq failed ({proc.returncode}): {proc.stderr}")

        if not out_txt.is_file() or out_txt.stat().st_size == 0:
            raise RuntimeError(f"preseq failed to generate output at {out_txt}")

        return (str(out_txt),)


NODE_CLASS_MAPPINGS = {"Preseq": Preseq}
NODE_DISPLAY_NAME_MAPPINGS = {"Preseq": "Preseq: Library Complexity Estimation"}


# Backward compatibility aliases
PreseqNode = Preseq

__all__ = ["Preseq",
    "PreseqNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
