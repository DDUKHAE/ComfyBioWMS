"""UMI-tools (Tools for handling Unique Molecular Identifiers in NGS data) node.

Python packages: umi_tools
External binaries: umi_tools
Galaxy wrapper: galaxyproject/tools-iuc tools/umi_tools/
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
        raise RuntimeError(f"UMI-tools exited with code {e.returncode}: {shlex.join(argv)}") from e


class UmiToolsExtract:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("extracted_fwd", "extracted_rev", "log_file")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reads_fwd": ("STRING", {"default": ""}),
                "bc_pattern": ("STRING", {"default": "NNNNNNNN"}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "bc_pattern2": ("STRING", {"default": ""}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        reads_fwd: str,
        bc_pattern: str,
        output_dir: str = "",
        reads_rev: str = "",
        bc_pattern2: str = "",
        extra_command: str = "",
    ):
        fwd_path = _file(reads_fwd, "Forward reads FASTQ")
        rev_path = _file(reads_rev, "Reverse reads FASTQ") if reads_rev.strip() else None

        executable = shutil.which("umi_tools")
        if not executable:
            raise RuntimeError(
                "umi_tools executable not found on PATH; install bioconda package umi_tools"
            )

        out = _output_dir("UmiToolsExtract", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        fwd_out = out / f"umi_extracted_{fwd_path.name}"
        log_file = out / "umi_extract.log"

        argv = [
            executable,
            "extract",
            "--bc-pattern",
            bc_pattern.strip(),
            "-I",
            str(fwd_path),
            "-S",
            str(fwd_out),
            "-L",
            str(log_file),
        ]
        rev_out = None
        if rev_path:
            rev_out = out / f"umi_extracted_{rev_path.name}"
            argv += ["--read2-in", str(rev_path), "--read2-out", str(rev_out)]
            if bc_pattern2.strip():
                argv += ["--bc-pattern2", bc_pattern2.strip()]

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not fwd_out.is_file() or fwd_out.stat().st_size == 0:
            raise RuntimeError(f"UMI-tools extract failed to produce output at {fwd_out}")

        return (str(fwd_out), str(rev_out) if rev_out and rev_out.is_file() else "", str(log_file))


class UmiToolsDedup:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("dedup_bam", "log_file")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "indexed_bam": ("STRING", {"default": ""}),
            },
            "optional": {
                "paired": ("BOOLEAN", {"default": False}),
                "method": (["directional", "unique", "percentile", "cluster"], {"default": "directional"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        indexed_bam: str,
        output_dir: str = "",
        paired: bool = False,
        method: str = "directional",
        extra_command: str = "",
    ):
        bam_path = _file(indexed_bam, "Indexed BAM")

        executable = shutil.which("umi_tools")
        if not executable:
            raise RuntimeError(
                "umi_tools executable not found on PATH; install bioconda package umi_tools"
            )

        out = _output_dir("UmiToolsDedup", output_dir)
        out.mkdir(parents=True, exist_ok=True)


        dedup_bam = out / f"dedup_{bam_path.name}"
        log_file = out / "umi_dedup.log"
        stats_prefix = str(out / "dedup_stats")

        argv = [
            executable,
            "dedup",
            "-I",
            str(bam_path),
            "-S",
            str(dedup_bam),
            "-L",
            str(log_file),
            "--method",
            method,
            "--output-stats",
            stats_prefix,
        ]
        if paired:
            argv.append("--paired")

        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        _run(argv, out)

        if not dedup_bam.is_file() or dedup_bam.stat().st_size == 0:
            raise RuntimeError(f"UMI-tools dedup failed to produce output at {dedup_bam}")

        return (str(dedup_bam), str(log_file))


NODE_CLASS_MAPPINGS = {
    "UmiToolsExtract": UmiToolsExtract,
    "UmiToolsDedup": UmiToolsDedup,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "UmiToolsExtract": "UMI-tools: Extract UMI to Header",
    "UmiToolsDedup": "UMI-tools: Deduplicate BAM",
}


# Backward compatibility aliases
UmiToolsExtractNode = UmiToolsExtract
UmiToolsDedupNode = UmiToolsDedup

__all__ = [
    "UmiToolsExtract",
    "UmiToolsExtractNode",
    "UmiToolsDedup",
    "UmiToolsDedupNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
