"""Trim Galore! quality and adapter trimming node.

Python packages: none
External binaries: trim_galore, cutadapt, (optional: fastqc)
Galaxy wrapper: galaxyproject/tools-iuc tools/trim_galore/trim_galore.xml
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
        raise RuntimeError(f"Trim Galore! exited with code {e.returncode}: {shlex.join(argv)}") from e


class TrimGalore:
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("trimmed_reads_fwd", "trimmed_reads_rev", "trimming_report")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "quality": ("INT", {"default": 20, "min": 0, "max": 50}),
                "adapter": ("STRING", {"default": ""}),
                "adapter2": ("STRING", {"default": ""}),
                "stringency": ("INT", {"default": 1, "min": 1, "max": 15}),
                "min_length": ("INT", {"default": 20, "min": 0, "max": 1000}),
                "cores": ("INT", {"default": 2, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        quality: int = 20,
        adapter: str = "",
        adapter2: str = "",
        stringency: int = 1,
        min_length: int = 20,
        cores: int = 2,
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None

        executable = shutil.which("trim_galore")
        if not executable:
            raise RuntimeError(
                "Trim Galore! executable not found on PATH; install bioconda package trim-galore"
            )

        out = _output_dir("TrimGalore", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[TrimGalore] ignored node-managed extra options: {' '.join(ignored)}", file=sys.stderr)

        argv = [
            executable,
            "--output_dir", str(out),
            "--quality", str(quality),
            "--stringency", str(stringency),
            "--length", str(min_length),
            "--cores", str(cores),
        ]
        if adapter.strip():
            argv.extend(["--adapter", adapter.strip()])
        if adapter2.strip():
            argv.extend(["--adapter2", adapter2.strip()])

        if rev_path:
            argv.extend(["--paired", str(fwd_path), str(rev_path)])
        else:
            argv.append(str(fwd_path))

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        report = next(out.glob("*trimming_report.txt"), None)
        report_file = str(report) if report else ""

        def _get_output(*patterns: str, err_msg: str) -> str:
            matched = [p for pat in patterns for p in out.glob(pat)]
            if not matched:
                raise RuntimeError(err_msg)
            if matched[0].stat().st_size == 0:
                raise RuntimeError("Trim Galore! produced empty trimmed fastq output")
            return str(matched[0])

        if rev_path:
            err = "Trim Galore! did not generate paired output fastq files"
            fwd = _get_output("*_val_1.fq*", "*_val_1.fastq*", err_msg=err)
            rev = _get_output("*_val_2.fq*", "*_val_2.fastq*", err_msg=err)
            return fwd, rev, report_file

        fwd = _get_output(
            "*_trimmed.fq*", "*_trimmed.fastq*", "*_val_1.fq*",
            err_msg="Trim Galore! did not generate trimmed fastq output",
        )
        return fwd, "", report_file


NODE_CLASS_MAPPINGS = {"TrimGalore": TrimGalore}
NODE_DISPLAY_NAME_MAPPINGS = {"TrimGalore": "Trim Galore!: Quality & Adapter Trimming"}


# Backward compatibility aliases
TrimGaloreNode = TrimGalore

__all__ = ["TrimGalore",
    "TrimGaloreNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
