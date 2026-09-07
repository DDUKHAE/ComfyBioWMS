"""SPAdes genome assembly node.

Python packages: none
External binary: spades.py==4.3.0 (Conda: bioconda::spades=4.3.0; Apt: spades)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/spades/spades.xml and tools/spades/macros.xml (SPAdes 4.3.0+galaxy1)
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



def _run(argv, cwd):
    try:
        subprocess.run(argv, cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"spades.py exited with code {e.returncode}: {shlex.join(argv)}") from e


def _nonempty(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"spades.py did not create a nonempty {label}: {path}")
    return str(path)


def _validate_contigs(path: Path) -> str:
    _nonempty(path, "contigs.fasta")
    first_line = ""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if stripped:
                first_line = stripped
                break
    if not first_line.startswith(">"):
        raise RuntimeError(f"contigs.fasta is not a valid FASTA file: {path}")
    return str(path)


class Spades:
    CATEGORY = "ComfyBIO/Genome Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("contigs", "scaffolds", "spades_log")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "read1": ("STRING", {"default": ""}),
            },
            "optional": {
                "read2": ("STRING", {"default": ""}),
                "unpaired_reads": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 16, "min": 1, "max": 256}),
                "memory_gb": ("INT", {"default": 16, "min": 1, "max": 1024}),
                "careful": ("BOOLEAN", {"default": False}),
                "sc": ("BOOLEAN", {"default": False}),
                "meta": ("BOOLEAN", {"default": False}),
                "isolate": ("BOOLEAN", {"default": False}),
                "only_assembler": ("BOOLEAN", {"default": False}),
                "cov_cutoff": ("STRING", {"default": "off"}),
                "kmers": ("STRING", {"default": "auto"}),
                "phred_offset": ("STRING", {"default": "auto"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        read1,
        output_dir: str = "",
        read2="",
        unpaired_reads="",
        threads=16,
        memory_gb=16,
        careful=False,
        sc=False,
        meta=False,
        isolate=False,
        only_assembler=False,
        cov_cutoff="off",
        kmers="auto",
        phred_offset="auto",
        extra_command="",
    ):
        read1_path = _file(read1, "Read 1")
        read2_path = _file(read2, "Read 2") if read2 else None
        unpaired_path = _file(unpaired_reads, "Unpaired reads") if unpaired_reads else None

        executable = shutil.which("spades.py")
        if not executable:
            raise RuntimeError(
                "spades.py executable not found on PATH; install conda package spades=4.3.0"
            )

        out = _output_dir("Spades", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(
                f"[spades] ignored node-managed extra options: {' '.join(ignored)}",
                file=sys.stderr,
            )

        argv = [
            executable,
            "-o",
            str(out),
        ]
        if read2_path:
            argv += ["-1", str(read1_path), "-2", str(read2_path)]
        else:
            argv += ["-s", str(read1_path)]
        if unpaired_path:
            argv += ["-s", str(unpaired_path)]

        argv += ["-t", str(threads), "-m", str(memory_gb)]

        if careful:
            argv.append("--careful")
        if sc:
            argv.append("--sc")
        if meta:
            argv.append("--meta")
        if isolate:
            argv.append("--isolate")
        if only_assembler:
            argv.append("--only-assembler")
        if cov_cutoff and str(cov_cutoff).strip().lower() != "off":
            argv += ["--cov-cutoff", str(cov_cutoff).strip()]
        if kmers and str(kmers).strip().lower() != "auto":
            argv += ["-k", str(kmers).strip()]
        if phred_offset and str(phred_offset).strip().lower() != "auto":
            argv += ["--phred-offset", str(phred_offset).strip()]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        contigs_path = out / "contigs.fasta"
        scaffolds_path = out / "scaffolds.fasta"
        log_path = out / "spades.log"

        scaffolds = (
            str(scaffolds_path)
            if scaffolds_path.is_file() and scaffolds_path.stat().st_size > 0
            else ""
        )
        return (
            _validate_contigs(contigs_path),
            scaffolds,
            _nonempty(log_path, "spades.log"),
        )


NODE_CLASS_MAPPINGS = {"Spades": Spades}
NODE_DISPLAY_NAME_MAPPINGS = {"Spades": "SPAdes: Assemble Genome"}


# Backward compatibility aliases
SpadesNode = Spades

__all__ = ["Spades",
    "SpadesNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
