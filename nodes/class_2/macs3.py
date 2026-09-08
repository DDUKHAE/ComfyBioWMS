"""MACS3 callpeak peak calling node.

Python packages: none
External binary: macs3==3.0.4 (Conda: bioconda::macs3=3.0.4; PyPI: macs3)
Galaxy wrapper: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/macs2/macs2_callpeak.xml
"""

from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Optional, Tuple

# Exact MACS3 spellings:
# -1 indicates variable arity (consumes non-option arguments until next option starting with '-')
FORMAT_CHOICES: List[str] = [
    "AUTO",
    "BAM",
    "SAM",
    "BED",
    "ELAND",
    "ELANDMULTI",
    "ELANDEXPORT",
    "BOWTIE",
    "BAMPE",
    "BEDPE",
]


def _file(value: str, label: str) -> Path:
    if not value or not str(value).strip():
        raise ValueError(f"{label} cannot be empty")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} file not found: {path}")
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



def build_macs3_callpeak_argv(
    executable: str = "macs3",
    treatment: str = "",
    output_dir: str = "",
    sample_name: str = "NA",
    control: str = "",
    format: str = "AUTO",
    genome_size: str = "hs",
    tag_size: int = 0,
    cutoff_mode: str = "qvalue",
    cutoff_value: float = 0.05,
    keep_dup: str = "1",
    nomodel: bool = False,
    extsize: int = 200,
    shift: int = 0,
    mfold_lower: int = 5,
    mfold_upper: int = 50,
    bandwidth: int = 300,
    broad: bool = False,
    broad_cutoff: float = 0.1,
    call_summits: bool = False,
    bedgraph: bool = False,
    spmr: bool = False,
    nolambda: bool = False,
    slocal: int = 1000,
    llocal: int = 10000,
    scale_to: str = "small",
    extra_args: Optional[List[str]] = None,
) -> List[str]:
    """Pure helper to construct macs3 callpeak command-line arguments."""
    argv = [
        executable,
        "callpeak",
        "-t",
        str(treatment),
        "--outdir",
        str(output_dir),
        "-n",
        str(sample_name),
        "-f",
        str(format),
        "-g",
        str(genome_size),
    ]

    if control:
        argv.extend(["-c", str(control)])

    if tag_size > 0:
        argv.extend(["-s", str(tag_size)])

    if str(cutoff_mode).lower() == "pvalue":
        val_str = f"{cutoff_value:g}" if isinstance(cutoff_value, (int, float)) else str(cutoff_value)
        argv.extend(["-p", val_str])
    else:
        val_str = f"{cutoff_value:g}" if isinstance(cutoff_value, (int, float)) else str(cutoff_value)
        argv.extend(["-q", val_str])

    if keep_dup:
        argv.extend(["--keep-dup", str(keep_dup)])

    is_paired_end = str(format).upper() in ("BAMPE", "BEDPE")

    if not is_paired_end:
        # Single-end: use nomodel plus shift/extsize or model plus mfold/bw
        if nomodel:
            argv.append("--nomodel")
            argv.extend(["--extsize", str(extsize)])
            if shift != 0:
                argv.extend(["--shift", str(shift)])
        else:
            argv.extend(["-m", str(int(mfold_lower)), str(int(mfold_upper))])
            if bandwidth:
                argv.extend(["--bw", str(bandwidth)])

    if broad:
        argv.append("--broad")
        b_val = f"{broad_cutoff:g}" if isinstance(broad_cutoff, (int, float)) else str(broad_cutoff)
        argv.extend(["--broad-cutoff", b_val])
    elif call_summits:
        argv.append("--call-summits")

    if bedgraph:
        argv.append("-B")
        if spmr:
            argv.append("--SPMR")

    if nolambda:
        argv.append("--nolambda")
    elif control:
        # slocal/llocal only emitted with a control and without nolambda
        if slocal is not None:
            argv.extend(["--slocal", str(slocal)])
        if llocal is not None:
            argv.extend(["--llocal", str(llocal)])

    st = str(scale_to).lower()
    selected_scale = st if st in ("small", "large") else "small"
    argv.extend(["--scale-to", selected_scale])

    if extra_args:
        argv.extend(extra_args)

    return argv


def resolve_output_paths(
    out_dir: Path,
    sample_name: str,
    broad: bool = False,
    bedgraph: bool = False,
) -> Dict[str, Optional[Path]]:
    """Determine expected output paths based on callpeak configuration."""
    primary_peaks = (
        out_dir / f"{sample_name}_peaks.broadPeak"
        if broad
        else out_dir / f"{sample_name}_peaks.narrowPeak"
    )
    peaks_xls = out_dir / f"{sample_name}_peaks.xls"
    summits = None if broad else (out_dir / f"{sample_name}_summits.bed")
    gapped_peaks = (out_dir / f"{sample_name}_peaks.gappedPeak") if broad else None
    treat_bedgraph = (out_dir / f"{sample_name}_treat_pileup.bdg") if bedgraph else None
    control_bedgraph = (out_dir / f"{sample_name}_control_lambda.bdg") if bedgraph else None

    return {
        "primary_peaks": primary_peaks,
        "peaks_xls": peaks_xls,
        "summits": summits,
        "gapped_peaks": gapped_peaks,
        "treat_bedgraph": treat_bedgraph,
        "control_bedgraph": control_bedgraph,
        "output_dir": out_dir,
    }


def _validate_bed_peaks(path: Path) -> str:
    """Validate that a BED-like peak output exists and has valid interval rows."""
    if not path.is_file():
        raise RuntimeError(f"Peak file does not exist: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"Peak file is empty: {path}")

    row_count = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if (
                not line
                or line.startswith("#")
                or line.startswith("track")
                or line.startswith("browser")
            ):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise RuntimeError(
                    f"Peak file {path} line {line_num} has fewer than 3 tab-separated fields: {line}"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as e:
                raise RuntimeError(
                    f"Peak file {path} line {line_num} has non-integer coordinates ({fields[1]}, {fields[2]}): {line}"
                ) from e
            if start < 0 or start >= end:
                raise RuntimeError(
                    f"Peak file {path} line {line_num} has invalid interval [start={start}, end={end}): {line}"
                )
            row_count += 1

    if row_count == 0:
        raise RuntimeError(f"Peak file has no interval rows: {path}")

    return str(path)


def _validate_nonempty_file(path: Path, label: str) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} file missing or empty: {path}")
    return str(path)


def _run(argv: List[str], cwd: Optional[Path] = None, max_tail_lines: int = 100) -> Tuple[str, str]:
    try:
        subprocess.run(argv, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"macs3 callpeak exited {e.returncode}: {shlex.join(argv)}") from e
    return "", ""


class Macs3Callpeak:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Epigenomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = (
        "primary_peaks",
        "peaks_xls",
        "summits",
        "gapped_peaks",
        "treat_bedgraph",
        "control_bedgraph",
        "output_dir",
    )

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "treatment": ("STRING", {"default": ""}),
            },
            "optional": {
                "control": ("STRING", {"default": ""}),
                "sample_name": ("STRING", {"default": "NA"}),
                "format": (
                    FORMAT_CHOICES,
                    {"default": "AUTO"},
                ),
                "genome_size": ("STRING", {"default": "hs"}),
                "tag_size": ("INT", {"default": 0, "min": 0, "max": 10000}),
                "cutoff_mode": (["qvalue", "pvalue"], {"default": "qvalue"}),
                "cutoff_value": (
                    "FLOAT",
                    {"default": 0.05, "min": 0.0, "max": 1.0, "step": 0.001},
                ),
                "keep_dup": ("STRING", {"default": "1"}),
                "nomodel": ("BOOLEAN", {"default": False}),
                "extsize": ("INT", {"default": 200, "min": 1, "max": 10000}),
                "shift": ("INT", {"default": 0, "min": -10000, "max": 10000}),
                "mfold_lower": (
                    "INT",
                    {"default": 5, "min": 1, "max": 1000},
                ),
                "mfold_upper": (
                    "INT",
                    {"default": 50, "min": 1, "max": 1000},
                ),
                "bandwidth": ("INT", {"default": 300, "min": 1, "max": 10000}),
                "broad": ("BOOLEAN", {"default": False}),
                "broad_cutoff": (
                    "FLOAT",
                    {"default": 0.1, "min": 0.0, "max": 1.0, "step": 0.001},
                ),
                "call_summits": ("BOOLEAN", {"default": False}),
                "bedgraph": ("BOOLEAN", {"default": False}),
                "spmr": ("BOOLEAN", {"default": False}),
                "nolambda": ("BOOLEAN", {"default": False}),
                "slocal": ("INT", {"default": 1000, "min": 0, "max": 1000000}),
                "llocal": (
                    "INT",
                    {"default": 10000, "min": 0, "max": 1000000},
                ),
                "scale_to": (
                    ["small", "large"],
                    {"default": "small"},
                ),
                "extra_command": (
                    "STRING",
                    {"default": "", "multiline": True},
                ),
            },
        }

    def run(
        self,
        treatment: str,
        output_dir: str = "",
        control: str = "",
        sample_name: str = "NA",
        format: str = "AUTO",
        genome_size: str = "hs",
        tag_size: int = 0,
        cutoff_mode: str = "qvalue",
        cutoff_value: float = 0.05,
        keep_dup: str = "1",
        nomodel: bool = False,
        extsize: int = 200,
        shift: int = 0,
        mfold_lower: int = 5,
        mfold_upper: int = 50,
        bandwidth: int = 300,
        broad: bool = False,
        broad_cutoff: float = 0.1,
        call_summits: bool = False,
        bedgraph: bool = False,
        spmr: bool = False,
        nolambda: bool = False,
        slocal: int = 1000,
        llocal: int = 10000,
        scale_to: str = "small",
        extra_command: str = "",
        **kwargs: Any,
    ) -> Tuple[str, str, str, str, str, str, str]:
        # 1. Validate output_dir before resolving input files or checking PATH
        out = _output_dir("Macs3Callpeak", output_dir)
        if out.exists() and not out.is_dir():
            raise ValueError(f"Output path exists and is not a directory: {out}")

        # 2. Validate sample_name before resolving input files or checking PATH
        clean_name = str(sample_name).strip() if sample_name is not None else ""
        if not clean_name:
            raise ValueError("Sample name cannot be empty")
        if "/" in clean_name or "\\" in clean_name:
            raise ValueError(f"Sample name cannot contain path separators: {clean_name}")

        # 3. Resolve and validate input files
        treatment_path = _file(treatment, "Treatment")
        control_path = _file(control, "Control") if control and control.strip() else None

        # 4. Check executable on PATH before directory creation
        executable = shutil.which("macs3")
        if not executable:
            raise RuntimeError(
                "macs3 executable not found on PATH; install bioconda package macs3=3.0.4 or pip install macs3"
            )

        # 5. Output directory created only after all input and sample-name validations pass
        out.mkdir(parents=True, exist_ok=True)

        # Filter extra command
        
        argv = build_macs3_callpeak_argv(
            executable=executable,
            treatment=str(treatment_path),
            output_dir=str(out),
            sample_name=clean_name,
            control=str(control_path) if control_path else "",
            format=format,
            genome_size=genome_size,
            tag_size=tag_size,
            cutoff_mode=cutoff_mode,
            cutoff_value=cutoff_value,
            keep_dup=keep_dup,
            nomodel=nomodel,
            extsize=extsize,
            shift=shift,
            mfold_lower=mfold_lower,
            mfold_upper=mfold_upper,
            bandwidth=bandwidth,
            broad=broad,
            broad_cutoff=broad_cutoff,
            call_summits=call_summits,
            bedgraph=bedgraph,
            spmr=spmr,
            nolambda=nolambda,
            slocal=slocal,
            llocal=llocal,
            scale_to=scale_to,
            extra_args=shlex.split(extra_command) if extra_command.strip() else [],
        )

        _run(argv, cwd=out)

        # Resolve output artifacts and validate
        artifacts = resolve_output_paths(
            out_dir=out,
            sample_name=clean_name,
            broad=broad,
            bedgraph=bedgraph,
        )

        primary_peaks_path = artifacts["primary_peaks"]
        peaks_xls_path = artifacts["peaks_xls"]
        summits_path = artifacts["summits"]
        gapped_peaks_path = artifacts["gapped_peaks"]
        treat_bdg_path = artifacts["treat_bedgraph"]
        control_bdg_path = artifacts["control_bedgraph"]

        # Validate required artifacts (raise if missing or empty)
        validated_primary = _validate_bed_peaks(primary_peaks_path)
        validated_xls = _validate_nonempty_file(peaks_xls_path, "peaks.xls")

        if broad:
            validated_gapped = _validate_bed_peaks(gapped_peaks_path)
            validated_summits = ""
        else:
            validated_summits = _validate_bed_peaks(summits_path)
            validated_gapped = ""

        if bedgraph:
            validated_treat_bdg = _validate_nonempty_file(treat_bdg_path, "treat_pileup.bdg")
            validated_ctrl_bdg = _validate_nonempty_file(control_bdg_path, "control_lambda.bdg")
        else:
            validated_treat_bdg = ""
            validated_ctrl_bdg = ""

        return (
            validated_primary,
            validated_xls,
            validated_summits,
            validated_gapped,
            validated_treat_bdg,
            validated_ctrl_bdg,
            str(out),
        )


NODE_CLASS_MAPPINGS = {"Macs3Callpeak": Macs3Callpeak}
NODE_DISPLAY_NAME_MAPPINGS = {"Macs3Callpeak": "MACS3: Call Peaks"}


# Backward compatibility aliases
Macs3CallpeakNode = Macs3Callpeak

__all__ = [
    "Macs3Callpeak",
    "Macs3CallpeakNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "build_macs3_callpeak_argv",
    "resolve_output_paths",
    "FORMAT_CHOICES",
]
