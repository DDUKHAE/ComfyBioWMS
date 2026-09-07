"""BEDTools genome coverage node.

Python packages: none
External binaries: bedtools
Galaxy wrapper: galaxyproject/tools-iuc tools/bedtools/genomeCoverageBed.xml
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



def _run_to_file(argv: list[str], out_file: Path) -> None:
    with open(out_file, "w", encoding="utf-8") as out_fp:
        proc = subprocess.run(argv, stdout=out_fp, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bedtools genomecov failed ({proc.returncode}): {proc.stderr}")


class BedtoolsGenomeCoverage:
    CATEGORY = "ComfyBIO/Visualization"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("fwd_bedgraph", "rev_bedgraph")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "split_strands": ("BOOLEAN", {"default": True}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        bam_file: str,
        output_dir: str = "",
        split_strands: bool = True,
        extra_command: str = "",
    ):
        bam_path = _file(bam_file, "Input BAM")

        executable = shutil.which("bedtools")
        if not executable:
            raise RuntimeError(
                "bedtools executable not found on PATH; install bioconda package bedtools"
            )

        out = _output_dir("BedtoolsGenomeCoverage", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        if ignored:
            print(f"[BEDTools] ignored managed extra options: {" ".join(ignored)}", file=sys.stderr)

        base_argv = [executable, "genomecov", "-ibam", str(bam_path), "-bg", "-split"] + kept

        if split_strands:
            fwd_bg = out / f"{bam_path.stem}_forward.bedgraph"
            rev_bg = out / f"{bam_path.stem}_reverse.bedgraph"

            _run_to_file(base_argv + ["-strand", "+"], fwd_bg)
            _run_to_file(base_argv + ["-strand", "-"], rev_bg)

            if not fwd_bg.is_file() or fwd_bg.stat().st_size == 0:
                raise RuntimeError(f"BEDTools failed to generate forward bedgraph at {fwd_bg}")

            return (str(fwd_bg), str(rev_bg))
        else:
            cov_bg = out / f"{bam_path.stem}_coverage.bedgraph"
            _run_to_file(base_argv, cov_bg)

            if not cov_bg.is_file() or cov_bg.stat().st_size == 0:
                raise RuntimeError(f"BEDTools failed to generate bedgraph at {cov_bg}")

            return (str(cov_bg), "")


NODE_CLASS_MAPPINGS = {"BedtoolsGenomeCoverage": BedtoolsGenomeCoverage}
NODE_DISPLAY_NAME_MAPPINGS = {"BedtoolsGenomeCoverage": "BEDTools: BAM to bedGraph Coverage"}


# Backward compatibility aliases
BedtoolsGenomeCoverageNode = BedtoolsGenomeCoverage

__all__ = ["BedtoolsGenomeCoverage",
    "BedtoolsGenomeCoverageNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
