"""fastp read preprocessing node.

Python packages: none
External binary: fastp==1.3.6 (Conda: bioconda::fastp=1.3.6; Apt: fastp)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/fastp/fastp.xml and tools/fastp/macros.xml (fastp 1.3.6+galaxy0)
"""

import json
import shlex
import shutil
import subprocess
import sys
import threading
from pathlib import Path


_MANAGED_OPTIONS = {
    "-i": 1,
    "--in1": 1,
    "-I": 1,
    "--in2": 1,
    "-o": 1,
    "--out1": 1,
    "-O": 1,
    "--out2": 1,
    "-w": 1,
    "--thread": 1,
    "-q": 1,
    "--qualified_quality_phred": 1,
    "-u": 1,
    "--unqualified_percent_limit": 1,
    "-n": 1,
    "--n_base_limit": 1,
    "-l": 1,
    "--length_required": 1,
    "--detect_adapter_for_pe": 0,
    "-c": 0,
    "--correction": 0,
    "-j": 1,
    "--json": 1,
    "-h": 1,
    "--html": 1,
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
        raise RuntimeError(f"fastp exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")
    return "".join(stdout), "".join(stderr)


def _nonempty(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"fastp did not create a nonempty {label}: {path}")
    return str(path)


class FastpNode:
    CATEGORY = "ComfyBIO/Read Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = ("read1", "read2", "json_report", "html_report")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "read1": ("STRING", {"default": ""}),
                "output_dir": ("STRING", {"default": "fastp"}),
            },
            "optional": {
                "read2": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "qualified_quality_phred": (
                    "INT",
                    {"default": 15, "min": 0, "max": 93},
                ),
                "unqualified_percent_limit": (
                    "INT",
                    {"default": 40, "min": 0, "max": 100},
                ),
                "n_base_limit": ("INT", {"default": 5, "min": 0, "max": 1000000}),
                "length_required": ("INT", {"default": 15, "min": 0, "max": 1000000}),
                "detect_adapter_for_pe": ("BOOLEAN", {"default": False}),
                "correction": ("BOOLEAN", {"default": False}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        read1,
        output_dir,
        read2="",
        threads=4,
        qualified_quality_phred=15,
        unqualified_percent_limit=40,
        n_base_limit=5,
        length_required=15,
        detect_adapter_for_pe=False,
        correction=False,
        extra_command="",
    ):
        read1_path = _file(read1, "Read 1")
        read2_path = _file(read2, "Read 2") if read2 else None
        executable = shutil.which("fastp")
        if not executable:
            raise RuntimeError(
                "fastp executable not found on PATH; install conda package fastp=1.3.6"
            )

        out = Path(output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)
        out1 = out / "R1.fastq.gz"
        out2 = out / "R2.fastq.gz"
        report_json = out / "fastp.json"
        report_html = out / "fastp.html"
        kept, ignored = _filter_extra(extra_command, _MANAGED_OPTIONS)
        if ignored:
            print(f"[fastp] ignored node-managed extra options: {' '.join(ignored)}", file=sys.stderr)

        argv = [
            executable,
            "-i",
            str(read1_path),
            "-o",
            str(out1),
            "-w",
            str(threads),
            "-q",
            str(qualified_quality_phred),
            "-u",
            str(unqualified_percent_limit),
            "-n",
            str(n_base_limit),
            "-l",
            str(length_required),
            "-j",
            str(report_json),
            "-h",
            str(report_html),
        ]
        if read2_path:
            argv += ["-I", str(read2_path), "-O", str(out2)]
        if detect_adapter_for_pe:
            argv.append("--detect_adapter_for_pe")
        if correction:
            argv.append("--correction")
        _run(argv + kept, out)

        json.loads(report_json.read_text())
        return (
            _nonempty(out1, "read 1 output"),
            _nonempty(out2, "read 2 output") if read2_path else "",
            _nonempty(report_json, "JSON report"),
            _nonempty(report_html, "HTML report"),
        )


NODE_CLASS_MAPPINGS = {"FastpNode": FastpNode}
NODE_DISPLAY_NAME_MAPPINGS = {"FastpNode": "fastp: Trim and QC Reads"}

__all__ = ["FastpNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
