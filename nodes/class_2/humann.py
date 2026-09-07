"""HUMAnN metabolic pathway and functional profiling node.

Python packages: none
External binaries: humann
Galaxy wrapper: galaxyproject/tools-iuc tools/humann/
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



class HUMAnN:
    CATEGORY = "ComfyBIO/Microbiome"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("genefamilies_tsv", "pathabundance_tsv")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, reads_file: str, output_dir: str = "", threads: int = 4, extra_command: str = ""):
        reads_path = _file(reads_file, "Reads File")

        executable = shutil.which("humann")
        if not executable:
            raise RuntimeError("humann executable not found on PATH; install bioconda package humann")

        out = _output_dir("HUMAnN", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        argv = [
            executable,
            "--input",
            str(reads_path),
            "--output",
            str(out),
            "--threads",
            str(threads),
        ]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"HUMAnN failed ({proc.returncode}): {proc.stderr}")

        gf_candidates = list(out.glob("*_genefamilies.tsv"))
        pa_candidates = list(out.glob("*_pathabundance.tsv"))

        if not gf_candidates or gf_candidates[0].stat().st_size == 0:
            raise RuntimeError(f"HUMAnN failed to produce genefamilies in {out}")

        return (str(gf_candidates[0]), str(pa_candidates[0]) if pa_candidates else "")


NODE_CLASS_MAPPINGS = {"HUMAnN": HUMAnN}
NODE_DISPLAY_NAME_MAPPINGS = {"HUMAnN": "HUMAnN: Functional Metabolic Profiler"}


# Backward compatibility aliases
HUMAnNNode = HUMAnN

__all__ = ["HUMAnN",
    "HUMAnNNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
