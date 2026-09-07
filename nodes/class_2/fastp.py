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
        raise RuntimeError(f"fastp exited with code {e.returncode}: {shlex.join(argv)}") from e


def _nonempty(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"fastp did not create a nonempty {label}: {path}")
    return str(path)


class Fastp:
    CATEGORY = "ComfyBIO/Read Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = ("read1", "read2", "json_report", "html_report")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "read1": ("STRING", {"default": ""}),
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
        output_dir: str = "",
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

        out = _output_dir("Fastp", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        out1 = out / "R1.fastq.gz"
        out2 = out / "R2.fastq.gz"
        report_json = out / "fastp.json"
        report_html = out / "fastp.html"
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
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        json.loads(report_json.read_text())
        return (
            _nonempty(out1, "read 1 output"),
            _nonempty(out2, "read 2 output") if read2_path else "",
            _nonempty(report_json, "JSON report"),
            _nonempty(report_html, "HTML report"),
        )


NODE_CLASS_MAPPINGS = {"Fastp": Fastp}
NODE_DISPLAY_NAME_MAPPINGS = {"Fastp": "fastp: Trim and QC Reads"}


# Backward compatibility aliases
FastpNode = Fastp

__all__ = ["Fastp",
    "FastpNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
