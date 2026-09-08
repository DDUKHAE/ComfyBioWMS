"""FastQC quality-report node.

Supports single FASTQ/BAM/SAM inputs as well as directory inputs (or comma-separated paths)
to automatically process multiple files in a single node.

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

try:
    from .common import discover_samples
except Exception:
    try:
        from nodes.class_2.common import discover_samples
    except Exception:
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, None)]


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
    OUTPUT_NODE = True
    OUPUT_NODE = True
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
        raw_val = str(input_file).strip()
        p = Path(raw_val).expanduser().resolve()

        if p.is_file():
            input_files = [p]
        elif p.is_dir():
            exts = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".sam")
            input_files = sorted([f for f in p.rglob("*") if f.is_file() and any(f.name.lower().endswith(e) for e in exts)])
            if not input_files:
                raise FileNotFoundError(f"No FASTQ/BAM/SAM files found in directory: {p}")
        elif "," in raw_val:
            input_files = [Path(x.strip()).expanduser().resolve() for x in raw_val.split(",") if x.strip()]
            for f in input_files:
                if not f.is_file():
                    raise FileNotFoundError(f"Input file not found: {f}")
        else:
            raise FileNotFoundError(f"Read input is neither a file nor a directory: {raw_val}")

        executable = shutil.which("fastqc")
        if not executable:
            raise RuntimeError(
                "FastQC executable not found on PATH; install conda package fastqc=0.12.1"
            )

        out = _output_dir("FastQC", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        extra_files = [
            ("--contaminants", contaminants),
            ("--adapters", adapters),
            ("--limits", limits),
        ]
        validated_extra_files = [
            (flag, _file(value, flag)) for flag, value in extra_files if value
        ]

        html_reports = []
        zip_reports = []

        for curr_file in input_files:
            file_format = "bam" if curr_file.suffix.lower() == ".bam" else "sam" if curr_file.suffix.lower() == ".sam" else "fastq"
            argv = [
                executable,
                "--outdir", str(out),
                "--threads", str(threads),
                "--kmers", str(kmers),
                "--format", file_format,
                "--extract",
            ]
            for flag, path_arg in validated_extra_files:
                argv += [flag, str(path_arg)]
            if nogroup:
                argv.append("--nogroup")
            if min_length:
                argv += ["--min_length", str(min_length)]
            if extra_command.strip():
                argv.extend(shlex.split(extra_command))

            argv.append(str(curr_file))
            _run(argv, out)

            stem = _stem(curr_file)
            html = out / f"{stem}_fastqc.html"
            archive = out / f"{stem}_fastqc.zip"
            if not html.is_file() or html.stat().st_size == 0:
                # Check for extracted dir or partial match
                possible_html = list(out.glob(f"*{stem}*fastqc.html"))
                if possible_html:
                    html = possible_html[0]
            if not archive.is_file() or archive.stat().st_size == 0:
                possible_zip = list(out.glob(f"*{stem}*fastqc.zip"))
                if possible_zip:
                    archive = possible_zip[0]

            html_reports.append(str(html))
            zip_reports.append(str(archive))

        if len(input_files) == 1 and p.is_file():
            return html_reports[0], zip_reports[0]
        return ",".join(html_reports), ",".join(zip_reports)


NODE_CLASS_MAPPINGS = {"FastQC": FastQC}
NODE_DISPLAY_NAME_MAPPINGS = {"FastQC": "FastQC: Read Quality Report"}


# Backward compatibility aliases
FastQCNode = FastQC

__all__ = ["FastQC",
    "FastQCNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
