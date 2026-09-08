"""dupRadar duplication rate vs expression level analysis node.

Python packages: none
External binaries: Rscript (with R package dupRadar installed)
"""

import shutil
import subprocess
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



class DupRadar:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("scatter_plot_png", "dup_matrix_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
                "gtf_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "strandedness": ("STRING", {"default": "auto"}),
                "threads": ("INT", {"default": 2, "min": 1, "max": 64}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        gtf_file: str,
        output_dir: str = "",
        strandedness: str = "auto",
        threads: int = 2,
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM")
        gtf_path = _file(gtf_file, "Annotation GTF")

        executable = shutil.which("Rscript")
        if not executable:
            raise RuntimeError("Rscript executable not found on PATH; install R runtime")

        out = _output_dir("DupRadar", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        strand_num = 0
        s = strandedness.strip().lower()
        if "rev" in s or "rf" in s or "isr" in s:
            strand_num = 2
        elif "fwd" in s or "fr" in s or "isf" in s:
            strand_num = 1

        plot_png = out / f"{bam_path.stem}_dupradar.png"
        dup_txt = out / f"{bam_path.stem}_dupradar_matrix.txt"

        r_script = f"""
        suppressPackageStartupMessages({{
            if (!requireNamespace("dupRadar", quietly = TRUE)) {{
                stop("dupRadar R package not installed")
            }}
            library(dupRadar)
        }})
        dm <- analyzeDuprates(
            bam = "{bam_path}",
            gtf = "{gtf_path}",
            stranded = {strand_num},
            threads = {threads}
        )
        write.table(dm, file = "{dup_txt}", sep = "\t", quote = FALSE, row.names = FALSE)
        png("{plot_png}", width = 800, height = 800)
        dupradarPlot(dm)
        dev.off()
        """

        proc = subprocess.run([executable, "-e", r_script], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"dupRadar execution failed: {proc.stderr}")

        if not plot_png.is_file() or plot_png.stat().st_size == 0:
            raise RuntimeError(f"dupRadar failed to create plot at {plot_png}")

        return (str(plot_png), str(dup_txt) if dup_txt.is_file() else "")


NODE_CLASS_MAPPINGS = {"DupRadar": DupRadar}
NODE_DISPLAY_NAME_MAPPINGS = {"DupRadar": "dupRadar: Expression vs Duplication QC"}


# Backward compatibility aliases
DupRadarNode = DupRadar

__all__ = ["DupRadar",
    "DupRadarNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
