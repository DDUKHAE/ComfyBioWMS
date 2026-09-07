"""Kallisto near-optimal RNA-seq quantification node.

Python packages: none
External binaries: kallisto
Galaxy wrapper: galaxyproject/tools-iuc tools/kallisto/
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



class KallistoQuant:
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("abundance_tsv", "abundance_h5")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "kallisto_index": ("STRING", {"default": ""}),
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "fragment_length": ("FLOAT", {"default": 200.0, "min": 1.0, "max": 10000.0}),
                "sd": ("FLOAT", {"default": 20.0, "min": 0.1, "max": 1000.0}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        kallisto_index: str,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        fragment_length: float = 200.0,
        sd: float = 20.0,
        extra_command: str = "",
    ):
        idx_path = _file(kallisto_index, "Kallisto Index")
        fwd_path = _file(reads_fwd, "Forward reads")
        rev_path = _file(reads_rev, "Reverse reads") if reads_rev.strip() else None

        executable = shutil.which("kallisto")
        if not executable:
            raise RuntimeError("kallisto executable not found on PATH; install bioconda package kallisto")

        out = _output_dir("KallistoQuant", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        argv = [
            executable,
            "quant",
            "-i",
            str(idx_path),
            "-o",
            str(out),
            "-t",
            str(threads),
        ]
        if rev_path:
            argv += [str(fwd_path), str(rev_path)]
        else:
            argv += ["--single", "-l", str(fragment_length), "-s", str(sd), str(fwd_path)]

        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"kallisto quant failed ({proc.returncode}): {proc.stderr}")

        tsv_out = out / "abundance.tsv"
        h5_out = out / "abundance.h5"

        if not tsv_out.is_file() or tsv_out.stat().st_size == 0:
            raise RuntimeError(f"kallisto failed to produce abundance.tsv in {out}")

        return (str(tsv_out), str(h5_out) if h5_out.is_file() else "")


NODE_CLASS_MAPPINGS = {"KallistoQuant": KallistoQuant}
NODE_DISPLAY_NAME_MAPPINGS = {"KallistoQuant": "Kallisto: Pseudoalignment & Quant"}


# Backward compatibility aliases
KallistoQuantNode = KallistoQuant

__all__ = ["KallistoQuant",
    "KallistoQuantNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
