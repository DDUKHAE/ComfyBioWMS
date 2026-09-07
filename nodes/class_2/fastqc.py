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
        raise RuntimeError(f"FastQC exited with code {e.returncode}: {shlex.join(argv)}") from e


def _stem(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".bz2"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".fastq", ".fq", ".bam", ".sam"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


class FastQC:
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("html_report", "zip_report")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_file": ("STRING", {"default": ""}),
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
        output_dir: str = "",
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

        out = _output_dir("FastQC", output_dir)
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


NODE_CLASS_MAPPINGS = {"FastQC": FastQC}
NODE_DISPLAY_NAME_MAPPINGS = {"FastQC": "FastQC: Read Quality Report"}


# Backward compatibility aliases
FastQCNode = FastQC

__all__ = ["FastQC",
    "FastQCNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
