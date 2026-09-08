"""Hifiasm de novo haplotype-resolved assembler for PacBio HiFi reads.

Python packages: none
External binaries: hifiasm
Galaxy wrapper: galaxyproject/tools-iuc tools/hifiasm/
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
        raise RuntimeError(f"Hifiasm exited with code {e.returncode}: {shlex.join(argv)}") from e


class HifiasmAssemble:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("primary_contigs_gfa", "alternate_contigs_gfa")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "hifi_reads_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, hifi_reads_file: str, output_dir: str = "", threads: int = 4, extra_command: str = ""):
        reads_path = _file(hifi_reads_file, "HiFi reads file")

        executable = shutil.which("hifiasm")
        if not executable:
            raise RuntimeError("hifiasm executable not found on PATH; install bioconda package hifiasm")

        out = _output_dir("HifiasmAssemble", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        prefix = str(out / f"{reads_path.stem}_asm")

        argv = [executable, "-o", prefix, "-t", str(threads), str(reads_path)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        primary_gfa = Path(f"{prefix}.bp.p_ctg.gfa")
        alt_gfa = Path(f"{prefix}.bp.a_ctg.gfa")

        if not primary_gfa.is_file() or primary_gfa.stat().st_size == 0:
            # Fallback check for alternate file naming
            primary_candidates = list(out.glob("*p_ctg.gfa"))
            if not primary_candidates:
                raise RuntimeError(f"Hifiasm failed to produce primary contigs GFA in {out}")
            primary_gfa = primary_candidates[0]

        return (str(primary_gfa), str(alt_gfa) if alt_gfa.is_file() else "")


NODE_CLASS_MAPPINGS = {"HifiasmAssemble": HifiasmAssemble}
NODE_DISPLAY_NAME_MAPPINGS = {"HifiasmAssemble": "Hifiasm: HiFi Haplotype-resolved Assembler"}


# Backward compatibility aliases
HifiasmAssembleNode = HifiasmAssemble

__all__ = ["HifiasmAssemble",
    "HifiasmAssembleNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
