"""Mosdepth fast BAM/CRAM depth and coverage calculation node.

Python packages: none
External binaries: mosdepth
Galaxy wrapper: galaxyproject/tools-iuc tools/mosdepth/
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
        raise RuntimeError(f"Mosdepth exited with code {e.returncode}: {shlex.join(argv)}") from e


class Mosdepth:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Genomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("summary_txt", "global_dist_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "window_size": ("INT", {"default": 500, "min": 0, "max": 100000}),
                "no_per_base": ("BOOLEAN", {"default": True}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 64}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        output_dir: str = "",
        window_size: int = 500,
        no_per_base: bool = True,
        threads: int = 4,
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM/CRAM")

        executable = shutil.which("mosdepth")
        if not executable:
            raise RuntimeError("mosdepth executable not found on PATH; install bioconda package mosdepth")

        out = _output_dir("Mosdepth", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        prefix = str(out / f"{bam_path.stem}_depth")

        argv = [executable, "-t", str(threads)]
        if window_size > 0:
            argv += ["-b", str(window_size)]
        if no_per_base:
            argv.append("-n")

        argv += [prefix, str(bam_path)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        summary_file = Path(f"{prefix}.mosdepth.summary.txt")
        dist_file = Path(f"{prefix}.mosdepth.global.dist.txt")

        if not summary_file.is_file() or summary_file.stat().st_size == 0:
            raise RuntimeError(f"Mosdepth failed to produce summary at {summary_file}")

        return (str(summary_file), str(dist_file) if dist_file.is_file() else "")


NODE_CLASS_MAPPINGS = {"Mosdepth": Mosdepth}
NODE_DISPLAY_NAME_MAPPINGS = {"Mosdepth": "Mosdepth: Fast BAM Coverage & Depth"}


# Backward compatibility aliases
MosdepthNode = Mosdepth

__all__ = ["Mosdepth",
    "MosdepthNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
