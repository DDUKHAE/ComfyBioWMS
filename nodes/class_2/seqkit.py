"""SeqKit FASTA/FASTQ sequence manipulation node.

Python packages: none
External binaries: seqkit
Galaxy wrapper: galaxyproject/tools-iuc tools/seqkit/
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



class SeqKitStats:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Genomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("stats_tsv",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "sequence_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "all_stats": ("BOOLEAN", {"default": True}),
                "threads": ("INT", {"default": 2, "min": 1, "max": 64}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        sequence_file: str,
        output_dir: str = "",
        all_stats: bool = True,
        threads: int = 2,
        extra_command: str = "",
    ):
        seq_path = _file(sequence_file, "Sequence File (FASTA/Q)")

        executable = shutil.which("seqkit")
        if not executable:
            raise RuntimeError("seqkit executable not found on PATH; install bioconda package seqkit")

        out = _output_dir("SeqKitStats", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        out_tsv = out / f"{seq_path.stem}_seqkit_stats.tsv"

        argv = [executable, "stats", "-T", "-j", str(threads)]
        if all_stats:
            argv.append("-a")
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())
        argv.append(str(seq_path))

        with open(out_tsv, "w", encoding="utf-8") as wf:
            proc = subprocess.run(argv, stdout=wf, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            raise RuntimeError(f"seqkit stats failed ({proc.returncode}): {proc.stderr}")

        if not out_tsv.is_file() or out_tsv.stat().st_size == 0:
            raise RuntimeError(f"SeqKit failed to produce stats file at {out_tsv}")

        return (str(out_tsv),)


NODE_CLASS_MAPPINGS = {"SeqKitStats": SeqKitStats}
NODE_DISPLAY_NAME_MAPPINGS = {"SeqKitStats": "SeqKit: Sequence File Statistics"}


# Backward compatibility aliases
SeqKitStatsNode = SeqKitStats

__all__ = ["SeqKitStats",
    "SeqKitStatsNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
