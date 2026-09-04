"""FastQC quality-report node.

Python packages: none
External binaries: FastQC==0.12.1 and Java runtime
Conda: bioconda::fastqc=0.12.1, conda-forge::openjdk; Apt: fastqc, default-jre
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/fastqc/rgFastQC.xml (FastQC 0.74+galaxy1; executable dependency 0.12.1)
"""

import shlex
import shutil
import subprocess
import sys
import threading
from pathlib import Path


_MANAGED_OPTIONS = {
    "-o": 1,
    "--outdir": 1,
    "-t": 1,
    "--threads": 1,
    "-c": 1,
    "--contaminants": 1,
    "-a": 1,
    "--adapters": 1,
    "-l": 1,
    "--limits": 1,
    "--nogroup": 0,
    "--min_length": 1,
    "-k": 1,
    "--kmers": 1,
    "-f": 1,
    "--format": 1,
    "--extract": 0,
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
        raise RuntimeError(f"FastQC exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")


def _stem(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".bz2"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".fastq", ".fq", ".bam", ".sam"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


class FastQCNode:
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("html_report", "zip_report")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_file": ("STRING", {"default": ""}),
                "output_dir": ("STRING", {"default": "fastqc"}),
            },
            "optional": {
                "threads": ("INT", {"default": 2, "min": 1, "max": 256}),
                "contaminants": ("STRING", {"default": ""}),
                "adapters": ("STRING", {"default": ""}),
                "limits": ("STRING", {"default": ""}),
                "nogroup": ("BOOLEAN", {"default": False}),
                "min_length": ("INT", {"default": 0, "min": 0, "max": 1000000}),
                "kmers": ("INT", {"default": 7, "min": 2, "max": 10}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        input_file,
        output_dir,
        threads=2,
        contaminants="",
        adapters="",
        limits="",
        nogroup=False,
        min_length=0,
        kmers=7,
        extra_command="",
    ):
        input_path = _file(input_file, "Read input")
        extra_files = [
            ("--contaminants", contaminants),
            ("--adapters", adapters),
            ("--limits", limits),
        ]
        validated_extra_files = [
            (flag, _file(value, flag)) for flag, value in extra_files if value
        ]
        executable = shutil.which("fastqc")
        if not executable:
            raise RuntimeError(
                "FastQC executable not found on PATH; install conda package fastqc=0.12.1"
            )

        out = Path(output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)
        kept, ignored = _filter_extra(extra_command, _MANAGED_OPTIONS)
        if ignored:
            print(f"[FastQC] ignored node-managed extra options: {' '.join(ignored)}", file=sys.stderr)

        file_format = "bam" if input_path.suffix.lower() == ".bam" else "sam" if input_path.suffix.lower() == ".sam" else "fastq"
        argv = [
            executable,
            "--outdir",
            str(out),
            "--threads",
            str(threads),
            "--kmers",
            str(kmers),
            "--format",
            file_format,
            "--extract",
        ]
        for flag, path in validated_extra_files:
            argv += [flag, str(path)]
        if nogroup:
            argv.append("--nogroup")
        if min_length:
            argv += ["--min_length", str(min_length)]
        _run(argv + kept + [str(input_path)], out)

        stem = _stem(input_path)
        html = out / f"{stem}_fastqc.html"
        archive = out / f"{stem}_fastqc.zip"
        for path in (html, archive):
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"FastQC did not create a nonempty report: {path}")
        return str(html), str(archive)


NODE_CLASS_MAPPINGS = {"FastQCNode": FastQCNode}
NODE_DISPLAY_NAME_MAPPINGS = {"FastQCNode": "FastQC: Read Quality Report"}

__all__ = ["FastQCNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
