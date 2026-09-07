"""MultiQC (Aggregate bioinformatics analysis reports across many samples) node.

Python packages: multiqc
External binaries: multiqc
Galaxy wrapper: galaxyproject/tools-iuc tools/multiqc/multiqc.xml
"""

import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def _dir(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"{label} is not a directory: {path}")
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
        raise RuntimeError(f"MultiQC exited with code {e.returncode}: {shlex.join(argv)}") from e


class MultiQC:
    CATEGORY = "ComfyBIO/Reporting"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("multiqc_html", "multiqc_data_dir")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "analysis_dir": ("STRING", {"default": ""}),
            },
            "optional": {
                "report_filename": ("STRING", {"default": "multiqc_report.html"}),
                "report_title": ("STRING", {"default": "RNA-Seq Analysis QC Summary"}),
                "config_yaml": ("STRING", {"default": ""}),
                "extra_scan_dirs": ("STRING", {"default": ""}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        analysis_dir: str,
        output_dir: str = "",
        report_filename: str = "multiqc_report.html",
        report_title: str = "RNA-Seq Analysis QC Summary",
        config_yaml: str = "",
        extra_scan_dirs: str = "",
        extra_command: str = "",
    ):
        scan_paths = []
        if analysis_dir.strip():
            for p in analysis_dir.split(","):
                clean = p.strip()
                if clean:
                    scan_paths.append(_dir(clean, "Analysis directory"))

        if extra_scan_dirs.strip():
            for p in extra_scan_dirs.split(","):
                clean = p.strip()
                if clean:
                    scan_paths.append(_dir(clean, "Extra scan directory"))

        if not scan_paths:
            raise ValueError("At least one valid scan directory must be provided to MultiQC")

        executable = shutil.which("multiqc")
        if not executable:
            raise RuntimeError("MultiQC executable not found on PATH; install package multiqc")

        out = _output_dir("MultiQC", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[MultiQC] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        argv = [
            executable,
            "--outdir",
            str(out),
            "--filename",
            report_filename,
            "--force",
            "--dirs",
        ]
        if report_title.strip():
            argv += ["--title", report_title.strip()]
        if config_yaml.strip():
            cfg_path = Path(config_yaml).expanduser().resolve()
            if cfg_path.is_file():
                argv += ["--config", str(cfg_path)]

        for p in scan_paths:
            argv.append(str(p))

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        html_out = out / report_filename
        # MultiQC data dir name default: <report_name_without_ext>_data
        stem = Path(report_filename).stem
        data_dir = out / f"{stem}_data"

        if not html_out.is_file() or html_out.stat().st_size == 0:
            raise RuntimeError(f"MultiQC failed to generate HTML report in {out}")

        return (str(html_out), str(data_dir) if data_dir.is_dir() else str(out))


NODE_CLASS_MAPPINGS = {"MultiQC": MultiQC}
NODE_DISPLAY_NAME_MAPPINGS = {"MultiQC": "MultiQC: Comprehensive Analysis QC Report"}


# Backward compatibility aliases
MultiQCNode = MultiQC

__all__ = ["MultiQC",
    "MultiQCNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
