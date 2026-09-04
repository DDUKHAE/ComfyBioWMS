"""BWA-MEM2 indexing and alignment nodes.

Python packages: none
External binary: bwa-mem2==2.3 (Conda: bioconda::bwa-mem2=2.3;
Apt: bwa-mem2 where provided by the distribution)
Galaxy wrappers: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/bwa_mem2/bwa-mem2-idx.xml and tools/bwa_mem2/bwa-mem2.xml
"""

import shlex
import shutil
import subprocess
import sys
import threading
from pathlib import Path


_ALIGN_MANAGED_OPTIONS = {
    "-t": 1,
    "-x": 1,
    "-k": 1,
    "-w": 1,
    "-T": 1,
    "-R": 1,
}


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _executable() -> str:
    executable = shutil.which("bwa-mem2")
    if not executable:
        raise RuntimeError(
            "bwa-mem2 executable not found on PATH; install conda package bwa-mem2=2.3"
        )
    return executable


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


def _run_logs(argv, cwd):
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
        raise RuntimeError(f"bwa-mem2 exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")


def _run_alignment(argv, output_sam, cwd):
    stderr = []
    with output_sam.open("wb") as output:
        process = subprocess.Popen(
            argv,
            cwd=cwd,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        for line in iter(process.stderr.readline, ""):
            stderr.append(line)
            print(line, end="", file=sys.stderr, flush=True)
        code = process.wait()
    if code:
        output_sam.unlink(missing_ok=True)
        raise RuntimeError(f"bwa-mem2 exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")


class BwaMem2IndexNode:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("indexed_reference",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reference_fasta": ("STRING", {"default": ""}),
                "output_dir": ("STRING", {"default": "bwa_mem2_index"}),
            },
            "optional": {
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, reference_fasta, output_dir, extra_command=""):
        source = _file(reference_fasta, "Reference FASTA")
        executable = _executable()
        out = Path(output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)
        reference = out / source.name
        if source.resolve() != reference.resolve():
            shutil.copy2(source, reference)
        argv = [executable, "index", *shlex.split(extra_command), str(reference)]
        _run_logs(argv, out)
        for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"):
            sidecar = Path(str(reference) + suffix)
            if not sidecar.is_file() or sidecar.stat().st_size == 0:
                raise RuntimeError(f"bwa-mem2 did not create index sidecar: {sidecar}")
        return (str(reference),)


class BwaMem2AlignNode:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("sam_file",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "indexed_reference": ("STRING", {"default": ""}),
                "read1": ("STRING", {"default": ""}),
                "output_sam": ("STRING", {"default": "alignment/reads.sam"}),
            },
            "optional": {
                "read2": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "preset": (
                    ["illumina", "pacbio", "ont2d", "intractg"],
                    {"default": "illumina"},
                ),
                "minimum_seed_length": (
                    "INT",
                    {"default": 19, "min": 1, "max": 1000000},
                ),
                "band_width": ("INT", {"default": 100, "min": 1, "max": 1000000}),
                "minimum_score": ("INT", {"default": 30, "min": 0, "max": 1000000}),
                "read_group": ("STRING", {"default": ""}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        indexed_reference,
        read1,
        output_sam,
        read2="",
        threads=1,
        preset="illumina",
        minimum_seed_length=19,
        band_width=100,
        minimum_score=30,
        read_group="",
        extra_command="",
    ):
        reference = _file(indexed_reference, "Indexed reference")
        read1_path = _file(read1, "Read 1")
        read2_path = _file(read2, "Read 2") if read2 else None
        executable = _executable()
        output = Path(output_sam).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        kept, ignored = _filter_extra(extra_command, _ALIGN_MANAGED_OPTIONS)
        if ignored:
            print(
                f"[bwa-mem2] ignored node-managed extra options: {' '.join(ignored)}",
                file=sys.stderr,
            )

        argv = [executable, "mem", "-t", str(threads)]
        if preset != "illumina":
            argv += ["-x", preset]
        argv += [
            "-k",
            str(minimum_seed_length),
            "-w",
            str(band_width),
            "-T",
            str(minimum_score),
        ]
        if read_group:
            argv += ["-R", read_group]
        argv += kept + [str(reference), str(read1_path)]
        if read2_path:
            argv.append(str(read2_path))
        _run_alignment(argv, output, output.parent)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"bwa-mem2 did not create a nonempty SAM file: {output}")
        with output.open(encoding="utf-8") as handle:
            if not handle.readline().startswith("@"):
                raise RuntimeError(f"bwa-mem2 output is missing a SAM header: {output}")
        return (str(output),)


NODE_CLASS_MAPPINGS = {
    "BwaMem2IndexNode": BwaMem2IndexNode,
    "BwaMem2AlignNode": BwaMem2AlignNode,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "BwaMem2IndexNode": "BWA-MEM2: Index Reference",
    "BwaMem2AlignNode": "BWA-MEM2: Align Reads",
}

__all__ = [
    "BwaMem2IndexNode",
    "BwaMem2AlignNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
