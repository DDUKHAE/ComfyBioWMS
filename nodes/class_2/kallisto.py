"""Kallisto pseudoalignment and quantification node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binaries: kallisto
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


class KallistoQuant:
    OUTPUT_NODE = True
    OUPUT_NODE = True
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
        samples = discover_samples(reads_fwd, reads_rev)

        executable = shutil.which("kallisto")
        if not executable:
            raise RuntimeError("kallisto executable not found on PATH; install bioconda package kallisto")

        base_out = _output_dir("KallistoQuant", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(reads_fwd).strip()).is_file()
        tsv_list, h5_list = [], []

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            argv = [
                executable,
                "quant",
                "-i", str(idx_path),
                "-o", str(sample_out),
                "-t", str(threads),
            ]
            if rev_path:
                argv += [str(fwd_path), str(rev_path)]
            else:
                argv += ["--single", "-l", str(fragment_length), "-s", str(sd), str(fwd_path)]

            if extra_command.strip():
                argv += shlex.split(extra_command.strip())

            proc = subprocess.run(argv, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"kallisto quant failed for {sample_id} ({proc.returncode}): {proc.stderr}")

            tsv_out = sample_out / "abundance.tsv"
            h5_out = sample_out / "abundance.h5"

            if not tsv_out.is_file() or tsv_out.stat().st_size == 0:
                raise RuntimeError(f"kallisto failed to produce abundance.tsv in {sample_out}")

            tsv_list.append(str(tsv_out))
            h5_list.append(str(h5_out) if h5_out.is_file() else "")

        if is_single:
            return tsv_list[0], h5_list[0]
        return ",".join(tsv_list), ",".join(h5_list)


NODE_CLASS_MAPPINGS = {"KallistoQuant": KallistoQuant}
NODE_DISPLAY_NAME_MAPPINGS = {"KallistoQuant": "Kallisto: Pseudoalignment & Quant"}


# Backward compatibility aliases
KallistoQuantNode = KallistoQuant

__all__ = ["KallistoQuant",
    "KallistoQuantNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
