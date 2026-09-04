"""QUAST genome assembly quality assessment node.

Python packages: none
External binary: quast.py==5.3.0 (Conda: bioconda::quast=5.3.0; Apt: quast)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/quast/quast.xml and tools/quast/macros.xml (QUAST 5.3.0+galaxy1)
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
    "-r": 1,
    "--reference": 1,
    "-g": 1,
    "--features": 1,
    "-m": 1,
    "--min-contig": 1,
    "-t": 1,
    "--threads": 1,
    "--large": 0,
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
            f"quast.py exited {code}: {shlex.join(argv)}\n{''.join(stderr)}"
        )
    return "".join(stdout), "".join(stderr)


def _nonempty(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"quast.py did not create a nonempty {label}: {path}")
    return str(path)


def _validate_report_tsv(path: Path) -> str:
    _nonempty(path, "report.tsv")
    content = path.read_text(errors="replace")
    if not any("contig" in line.lower() or "length" in line.lower() for line in content.splitlines()):
        raise RuntimeError(f"report.tsv does not contain valid QUAST metrics: {path}")
    return str(path)


def _validate_report_html(path: Path) -> str:
    _nonempty(path, "report.html")
    content = path.read_text(errors="replace")
    if "<html" not in content.lower() and "<!doctype html" not in content.lower():
        raise RuntimeError(f"report.html is not a valid HTML document: {path}")
    return str(path)


class QuastNode:
    CATEGORY = "ComfyBIO/Genome Assembly"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("report_tsv", "report_html")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "assembly_fasta": ("STRING", {"default": ""}),
                "output_dir": ("STRING", {"default": "quast_out"}),
            },
            "optional": {
                "reference": ("STRING", {"default": ""}),
                "features": ("STRING", {"default": ""}),
                "min_contig": ("INT", {"default": 500, "min": 0, "max": 1000000}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 256}),
                "large": ("BOOLEAN", {"default": False}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        assembly_fasta,
        output_dir,
        reference="",
        features="",
        min_contig=500,
        threads=4,
        large=False,
        extra_command="",
    ):
        assembly_path = _file(assembly_fasta, "Assembly FASTA")
        ref_path = _file(reference, "Reference") if reference else None
        feat_path = _file(features, "Features") if features else None

        executable = shutil.which("quast.py") or shutil.which("quast")
        if not executable:
            raise RuntimeError(
                "quast.py executable not found on PATH; install conda package quast=5.3.0"
            )

        out = Path(output_dir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)

        kept, ignored = _filter_extra(extra_command, _MANAGED_OPTIONS)
        if ignored:
            print(
                f"[quast] ignored node-managed extra options: {' '.join(ignored)}",
                file=sys.stderr,
            )

        argv = [
            executable,
            str(assembly_path),
            "-o",
            str(out),
            "-m",
            str(min_contig),
            "-t",
            str(threads),
        ]
        if ref_path:
            argv += ["-r", str(ref_path)]
        if feat_path:
            argv += ["-g", str(feat_path)]
        if large:
            argv.append("--large")

        _run(argv + kept, out)

        report_tsv = out / "report.tsv"
        report_html = out / "report.html"

        return (
            _validate_report_tsv(report_tsv),
            _validate_report_html(report_html),
        )


NODE_CLASS_MAPPINGS = {"QuastNode": QuastNode}
NODE_DISPLAY_NAME_MAPPINGS = {"QuastNode": "QUAST: Evaluate Assembly"}

__all__ = ["QuastNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
