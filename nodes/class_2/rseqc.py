"""RSeQC (RNA-seq Quality Control package) nodes.

Python packages: RSeQC
External binaries: geneBody_coverage.py, junction_saturation.py, infer_experiment.py
Galaxy wrapper: galaxyproject/tools-iuc tools/rseqc/
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
        raise RuntimeError(f"RSeQC tool exited with code {e.returncode}: {shlex.join(argv)}") from e


class RSeQCGeneBodyCoverage:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("coverage_curves_pdf", "coverage_curves_txt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
                "ref_gene_bed": ("STRING", {"default": ""}),
            },
            "optional": {
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, bam_file: str, ref_gene_bed: str, output_dir: str = "", extra_command: str = ""):
        bam_path = _file(bam_file, "Input BAM")
        bed_path = _file(ref_gene_bed, "Reference Gene BED12")

        executable = shutil.which("geneBody_coverage.py")
        if not executable:
            raise RuntimeError(
                "geneBody_coverage.py not found on PATH; install bioconda package rseqc"
            )

        out = _output_dir("RSeQCGeneBodyCoverage", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        prefix = str(out / f"{bam_path.stem}_gbc")

        argv = [executable, "-i", str(bam_path), "-r", str(bed_path), "-o", prefix]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        pdf_file = Path(f"{prefix}.geneBodyCoverage.curves.pdf")
        txt_file = Path(f"{prefix}.geneBodyCoverage.txt")

        return (
            str(pdf_file) if pdf_file.is_file() else "",
            str(txt_file) if txt_file.is_file() else "",
        )


class RSeQCJunctionSaturation:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("saturation_plot_pdf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
                "ref_gene_bed": ("STRING", {"default": ""}),
            },
            "optional": {
                "min_coverage": ("INT", {"default": 1, "min": 1, "max": 100}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, bam_file: str, ref_gene_bed: str, output_dir: str = "", min_coverage: int = 1, extra_command: str = ""):
        bam_path = _file(bam_file, "Input BAM")
        bed_path = _file(ref_gene_bed, "Reference Gene BED12")

        executable = shutil.which("junction_saturation.py")
        if not executable:
            raise RuntimeError(
                "junction_saturation.py not found on PATH; install bioconda package rseqc"
            )

        out = _output_dir("RSeQCJunctionSaturation", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        prefix = str(out / f"{bam_path.stem}_junction")

        argv = [executable, "-i", str(bam_path), "-r", str(bed_path), "-o", prefix, "-m", str(min_coverage)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        plot_file = Path(f"{prefix}.junctionSaturation_plot.pdf")
        return (str(plot_file) if plot_file.is_file() else "",)


class RSeQCInferExperiment:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("infer_experiment_txt",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
                "ref_gene_bed": ("STRING", {"default": ""}),
            },
            "optional": {
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, bam_file: str, ref_gene_bed: str, output_dir: str = "", extra_command: str = ""):
        bam_path = _file(bam_file, "Input BAM")
        bed_path = _file(ref_gene_bed, "Reference Gene BED12")

        executable = shutil.which("infer_experiment.py")
        if not executable:
            raise RuntimeError(
                "infer_experiment.py not found on PATH; install bioconda package rseqc"
            )

        out = _output_dir("RSeQCInferExperiment", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        out_txt = out / f"{bam_path.stem}_infer_experiment.txt"
        with open(out_txt, "w", encoding="utf-8") as wf:
            proc = subprocess.run([executable, "-i", str(bam_path), "-r", str(bed_path)] + kept, stdout=wf, stderr=subprocess.PIPE, text=True)

        if proc.returncode != 0:
            raise RuntimeError(f"infer_experiment.py failed: {proc.stderr}")

        return (str(out_txt),)


NODE_CLASS_MAPPINGS = {
    "RSeQCGeneBodyCoverage": RSeQCGeneBodyCoverage,
    "RSeQCJunctionSaturation": RSeQCJunctionSaturation,
    "RSeQCInferExperiment": RSeQCInferExperiment,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "RSeQCGeneBodyCoverage": "RSeQC: Gene Body Coverage",
    "RSeQCJunctionSaturation": "RSeQC: Junction Saturation",
    "RSeQCInferExperiment": "RSeQC: Infer Experiment Strandedness",
}


# Backward compatibility aliases
RSeQCGeneBodyCoverageNode = RSeQCGeneBodyCoverage
RSeQCJunctionSaturationNode = RSeQCJunctionSaturation
RSeQCInferExperimentNode = RSeQCInferExperiment

__all__ = [
    "RSeQCGeneBodyCoverage",
    "RSeQCGeneBodyCoverageNode",
    "RSeQCJunctionSaturation",
    "RSeQCJunctionSaturationNode",
    "RSeQCInferExperiment",
    "RSeQCInferExperimentNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
