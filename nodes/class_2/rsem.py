"""RSEM (RNA-Seq by Expectation-Maximization) node.

Python packages: none
External binaries: rsem-prepare-reference, rsem-calculate-expression
Galaxy wrapper: galaxyproject/tools-iuc tools/rsem/
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
        raise RuntimeError(f"RSEM exited with code {e.returncode}: {shlex.join(argv)}") from e


class RSEMCalculateExpression:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("genes_results_tsv", "isoforms_results_tsv")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "transcriptome_bam": ("STRING", {"default": ""}),
                "rsem_reference_prefix": ("STRING", {"default": ""}),
            },
            "optional": {
                "strandedness": ("STRING", {"default": "auto"}),
                "paired_end": ("BOOLEAN", {"default": True}),
                "threads": ("INT", {"default": 4, "min": 1, "max": 128}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        transcriptome_bam: str,
        rsem_reference_prefix: str,
        output_dir: str = "",
        strandedness: str = "auto",
        paired_end: bool = True,
        threads: int = 4,
        extra_command: str = "",
    ):
        bam_path = _file(transcriptome_bam, "Transcriptome BAM")

        executable = shutil.which("rsem-calculate-expression")
        if not executable:
            raise RuntimeError("rsem-calculate-expression executable not found on PATH; install bioconda package rsem")

        out = _output_dir("RSEMCalculateExpression", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        sample_prefix = str(out / f"{bam_path.stem}_rsem")

        argv = [
            executable,
            "--bam",
            "--num-threads",
            str(threads),
        ]
        if paired_end:
            argv.append("--paired-end")

        s = strandedness.strip().lower()
        if "rev" in s or "rf" in s or "isr" in s:
            argv += ["--strandedness", "reverse"]
        elif "fwd" in s or "fr" in s or "isf" in s:
            argv += ["--strandedness", "forward"]
        else:
            argv += ["--strandedness", "none"]

        argv += [str(bam_path), rsem_reference_prefix, sample_prefix]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        genes_out = Path(f"{sample_prefix}.genes.results")
        isoforms_out = Path(f"{sample_prefix}.isoforms.results")

        if not genes_out.is_file() or genes_out.stat().st_size == 0:
            raise RuntimeError(f"RSEM failed to produce genes results at {genes_out}")

        return (str(genes_out), str(isoforms_out) if isoforms_out.is_file() else "")


NODE_CLASS_MAPPINGS = {"RSEMCalculateExpression": RSEMCalculateExpression}
NODE_DISPLAY_NAME_MAPPINGS = {"RSEMCalculateExpression": "RSEM: Calculate Expression from BAM"}


# Backward compatibility aliases
RSEMCalculateExpressionNode = RSEMCalculateExpression

__all__ = ["RSEMCalculateExpression",
    "RSEMCalculateExpressionNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
