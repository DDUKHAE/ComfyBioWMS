"""Trimmomatic read trimming node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

Python packages: none
External binaries: trimmomatic
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


class Trimmomatic:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Read Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("trimmed_reads_fwd", "trimmed_reads_rev")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "leading": ("INT", {"default": 3, "min": 0, "max": 40}),
                "trailing": ("INT", {"default": 3, "min": 0, "max": 40}),
                "minlen": ("INT", {"default": 36, "min": 1, "max": 1000}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        threads: int = 4,
        leading: int = 3,
        trailing: int = 3,
        minlen: int = 36,
        extra_command: str = "",
    ):
        samples = discover_samples(reads_fwd, reads_rev)

        executable = shutil.which("trimmomatic")
        if not executable:
            raise RuntimeError("trimmomatic executable not found on PATH; install bioconda package trimmomatic")

        base_out = _output_dir("Trimmomatic", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(reads_fwd).strip()).is_file()
        fwd_list, rev_list = [], []

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            paired = bool(rev_path)
            mode = "PE" if paired else "SE"

            out_fwd = sample_out / f"{sample_id}_trimmed_R1.fq.gz"
            out_fwd_unpaired = sample_out / f"{sample_id}_unpaired_R1.fq.gz"
            out_rev = sample_out / f"{sample_id}_trimmed_R2.fq.gz" if paired else None
            out_rev_unpaired = sample_out / f"{sample_id}_unpaired_R2.fq.gz" if paired else None

            argv = [executable, mode, "-threads", str(threads)]
            if paired:
                argv += [str(fwd_path), str(rev_path), str(out_fwd), str(out_fwd_unpaired), str(out_rev), str(out_rev_unpaired)]
            else:
                argv += [str(fwd_path), str(out_fwd)]

            argv += [f"LEADING:{leading}", f"TRAILING:{trailing}", f"MINLEN:{minlen}"]
            if extra_command.strip():
                argv += shlex.split(extra_command.strip())

            proc = subprocess.run(argv, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"Trimmomatic failed for {sample_id} ({proc.returncode}): {proc.stderr}")

            if not out_fwd.is_file() or out_fwd.stat().st_size == 0:
                raise RuntimeError(f"Trimmomatic failed to produce output for {sample_id} at {out_fwd}")

            fwd_list.append(str(out_fwd))
            if out_rev and out_rev.is_file():
                rev_list.append(str(out_rev))

        if is_single:
            return str(fwd_list[0]), (str(rev_list[0]) if rev_list else "")
        return ",".join(fwd_list), ",".join(rev_list)


NODE_CLASS_MAPPINGS = {"Trimmomatic": Trimmomatic}
NODE_DISPLAY_NAME_MAPPINGS = {"Trimmomatic": "Trimmomatic: Read Trimming"}


# Backward compatibility aliases
TrimmomaticNode = Trimmomatic

__all__ = ["Trimmomatic",
    "TrimmomaticNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
