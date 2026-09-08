"""RNA-seq Library Strandedness Inference node (based on nf-core check_strandedness).

Python packages: none
External binaries: salmon
"""

import gzip
import json
import shlex
import shutil
import subprocess
from pathlib import Path

try:
    from .common import discover_samples
except Exception:
    try:
        from nodes.class_2.common import discover_samples
    except Exception:
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, Path(rev_input).expanduser().resolve() if rev_input.strip() else None)]


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



def _dir(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"{label} is not a directory: {path}")
    return path


def _subsample_fastq(in_path: Path, out_path: Path, max_reads: int) -> None:
    is_gz = in_path.name.endswith(".gz")
    open_in = gzip.open(in_path, "rt", encoding="utf-8") if is_gz else open(in_path, "r", encoding="utf-8")
    open_out = gzip.open(out_path, "wt", encoding="utf-8") if is_gz else open(out_path, "w", encoding="utf-8")

    reads_written = 0
    with open_in, open_out:
        while reads_written < max_reads:
            h = open_in.readline()
            if not h:
                break
            s = open_in.readline()
            p = open_in.readline()
            q = open_in.readline()
            open_out.write(h + s + p + q)
            reads_written += 1


class InferStrandedness:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("strandedness", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "salmon_index_dir": ("STRING", {"default": ""}),
                "reads_fwd": ("STRING", {"default": ""}),
            },
            "optional": {
                "reads_rev": ("STRING", {"default": ""}),
                "subsample_reads": ("INT", {"default": 1000000, "min": 10000, "max": 50000000}),
                "threads": ("INT", {"default": 2, "min": 1, "max": 64}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        salmon_index_dir: str,
        reads_fwd: str,
        output_dir: str = "",
        reads_rev: str = "",
        subsample_reads: int = 1000000,
        threads: int = 2,
        extra_command: str = "",
    ):
        idx_dir = _dir(salmon_index_dir, "Salmon Index")
        samples = discover_samples(reads_fwd, reads_rev)
        _, fwd_path, rev_path = samples[0]

        executable = shutil.which("salmon")
        if not executable:
            raise RuntimeError("Salmon executable not found on PATH; install bioconda package salmon")

        out = _output_dir("InferStrandedness", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # Subsample FASTQ files
        sub_fwd = out / f"subsample_fwd_{fwd_path.name}"
        _subsample_fastq(fwd_path, sub_fwd, subsample_reads)

        sub_rev = None
        if rev_path:
            sub_rev = out / f"subsample_rev_{rev_path.name}"
            _subsample_fastq(rev_path, sub_rev, subsample_reads)

        # Run salmon with automatic library detection (-l A)
        quant_out = out / "salmon_run"
        argv = [
            executable,
            "quant",
            "-i",
            str(idx_dir),
            "-l",
            "A",
            "-o",
            str(quant_out),
            "-p",
            str(threads),
        ]
        if sub_rev:
            argv += ["-1", str(sub_fwd), "-2", str(sub_rev)]
        else:
            argv += ["-r", str(sub_fwd)]

        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Salmon strandedness check failed: {proc.stderr}")

        # Parse lib_format_counts.json
        lib_json = quant_out / "lib_format_counts.json"
        if not lib_json.is_file():
            raise RuntimeError(f"Salmon did not produce lib_format_counts.json in {quant_out}")

        with open(lib_json, "r", encoding="utf-8") as f:
            counts = json.load(f)

        paired = bool(rev_path)
        if paired:
            fwd_count = counts.get("ISF", 0) + counts.get("OSF", 0) + counts.get("MSF", 0)
            rev_count = counts.get("ISR", 0) + counts.get("OSR", 0) + counts.get("MSR", 0)
            unstr_count = counts.get("IU", 0) + counts.get("OU", 0) + counts.get("MU", 0)
        else:
            fwd_count = counts.get("SF", 0)
            rev_count = counts.get("SR", 0)
            unstr_count = counts.get("U", 0)

        total = fwd_count + rev_count + unstr_count
        if total == 0:
            decision = "unstranded"
            ratio = 0.0
        else:
            fwd_ratio = fwd_count / total
            rev_ratio = rev_count / total

            if rev_ratio >= 0.8:
                decision = "reverse"
                ratio = rev_ratio
            elif fwd_ratio >= 0.8:
                decision = "forward"
                ratio = fwd_ratio
            else:
                decision = "unstranded"
                ratio = max(fwd_ratio, rev_ratio)

        summary = {
            "decision": decision,
            "confidence_ratio": round(ratio, 4),
            "paired": paired,
            "counts": counts,
            "fwd_count": fwd_count,
            "rev_count": rev_count,
            "unstranded_count": unstr_count,
            "total_mapped": total,
        }

        summary_path = out / "strandedness_summary.json"
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

        return (decision, str(summary_path))


NODE_CLASS_MAPPINGS = {"InferStrandedness": InferStrandedness}
NODE_DISPLAY_NAME_MAPPINGS = {"InferStrandedness": "RNA-Seq: Infer Library Strandedness"}


# Backward compatibility aliases
InferStrandednessNode = InferStrandedness

__all__ = ["InferStrandedness",
    "InferStrandednessNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
