"""Trimmomatic flexible read trimming node.

Python packages: none
External binaries: trimmomatic
Galaxy wrapper: galaxyproject/tools-iuc tools/trimmomatic/
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



class Trimmomatic:
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("trimmed_fwd", "trimmed_rev")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "leading": ("INT", {"default": 3, "min": 0, "max": 50}),
                "trailing": ("INT", {"default": 3, "min": 0, "max": 50}),
                "minlen": ("INT", {"default": 36, "min": 1, "max": 500}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        leading: int = 3,
        trailing: int = 3,
        minlen: int = 36,
        threads: int = 4,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None

        executable = shutil.which("trimmomatic")
        if not executable:
            raise RuntimeError("trimmomatic executable not found on PATH; install bioconda package trimmomatic")

        out = _output_dir("Trimmomatic", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        paired = bool(rev_path)
        mode = "PE" if paired else "SE"

        out_fwd = out / f"{fwd_path.stem}_trimmed.fq.gz"
        out_fwd_unpaired = out / f"{fwd_path.stem}_unpaired.fq.gz"
        out_rev = out / f"{rev_path.stem}_trimmed.fq.gz" if paired else None
        out_rev_unpaired = out / f"{rev_path.stem}_unpaired.fq.gz" if paired else None

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
            raise RuntimeError(f"Trimmomatic failed ({proc.returncode}): {proc.stderr}")

        if not out_fwd.is_file() or out_fwd.stat().st_size == 0:
            raise RuntimeError(f"Trimmomatic failed to produce output at {out_fwd}")

        return (str(out_fwd), str(out_rev) if out_rev and out_rev.is_file() else "")


NODE_CLASS_MAPPINGS = {"Trimmomatic": Trimmomatic}
NODE_DISPLAY_NAME_MAPPINGS = {"Trimmomatic": "Trimmomatic: Read Trimming"}


# Backward compatibility aliases
TrimmomaticNode = Trimmomatic

__all__ = ["Trimmomatic",
    "TrimmomaticNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
