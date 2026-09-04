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
import threading
from pathlib import Path


_MANAGED_OPTIONS = {
    "-o": 1,
    "--output-dir": 1,
    "-1": 1,
    "--pe1-1": 1,
    "--pe-1": 1,
    "-2": 1,
    "--pe1-2": 1,
    "--pe-2": 1,
    "-s": 1,
    "--pe1-s": 1,
    "--pe-s": 1,
    "--12": 1,
    "--pe1-12": 1,
    "--pe-12": 1,
    "-t": 1,
    "--threads": 1,
    "-m": 1,
    "--memory": 1,
    "--careful": 0,
    "--sc": 0,
    "--meta": 0,
    "--isolate": 0,
    "--only-assembler": 0,
    "--cov-cutoff": 1,
    "-k": 1,
    "--kmers": 1,
    "--phred-offset": 1,
}


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _filter_extra(text, managed):
    tokens, kept, ignored = shlex.split(text), [], []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        match = next(
            (
                flag
                for flag in managed
                if token == flag
                or token.startswith(flag + "=")
                or (
                    len(flag) == 2
                    and managed[flag] == 1
                    and token.startswith(flag)
                    and token != flag
                )
            ),
            None,
        )
        if match:
            ignored.append(token)
            i += 1 + (managed[match] if token == match else 0)
        else:
            kept.append(token)
            i += 1
    return kept, ignored


def _run(argv, cwd):
    process = subprocess.Popen(
        argv,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
    )
    stdout, stderr = [], []

    def drain(pipe, target, collected):
        for line in iter(pipe.readline, ""):
            collected.append(line)
            print(line, end="", file=target, flush=True)

    threads = [
        threading.Thread(target=drain, args=(process.stdout, sys.stdout, stdout)),
        threading.Thread(target=drain, args=(process.stderr, sys.stderr, stderr)),
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    code = process.wait()
    if code:
        raise RuntimeError(
            f"spades.py exited {code}: {shlex.join(argv)}\n{''.join(stderr)}"
        )
    return "".join(stdout), "".join(stderr)


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


class SpadesNode:
    CATEGORY = "ComfyBIO/Genome Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("contigs", "scaffolds", "spades_log")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "read1": ("STRING", {"default": ""}),
                "output_dir": ("STRING", {"default": "spades_out"}),
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
        output_dir,
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

        out = Path(output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)

        kept, ignored = _filter_extra(extra_command, _MANAGED_OPTIONS)
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

        _run(argv + kept, out)

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


NODE_CLASS_MAPPINGS = {"SpadesNode": SpadesNode}
NODE_DISPLAY_NAME_MAPPINGS = {"SpadesNode": "SPAdes: Assemble Genome"}

__all__ = ["SpadesNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
