"""fastp read preprocessing node.

Supports single FASTQ file inputs as well as directory inputs (or comma-separated paths)
to automatically discover and process multiple sample pairs/singles iteratively in a single node.

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

try:
    from .common import discover_samples
except Exception:
    try:
        from nodes.class_2.common import discover_samples
    except Exception:
        def discover_samples(fwd_input: str, rev_input: str = ""):
            p = Path(fwd_input).expanduser().resolve()
            return [(p.stem, p, Path(rev_input).expanduser().resolve() if rev_input.strip() else None)]


try:
    from bioflow.runtime.command_runner import BioCommandRunner, resolve_tool_environment, validate_extra_command_conflicts
    from bioflow.runtime.artifacts import compute_input_fingerprint, get_run_output_dir
except Exception:
    try:
        import sys
        sys.path.append(str(Path(__file__).resolve().parents[2] / "engine" / "src"))
        from bioflow.runtime.command_runner import BioCommandRunner, resolve_tool_environment, validate_extra_command_conflicts
        from bioflow.runtime.artifacts import compute_input_fingerprint, get_run_output_dir
    except Exception:
        BioCommandRunner = None
        resolve_tool_environment = None
        compute_input_fingerprint = None


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


def _run(argv, cwd, inputs=None, outputs=None, extra_command="", managed_flags=None):
    if BioCommandRunner is not None:
        BioCommandRunner.run(
            argv=argv,
            cwd=cwd,
            node_type="Fastp",
            inputs=inputs,
            outputs=outputs,
            extra_command=extra_command,
            managed_flags=managed_flags,
        )
    else:
        try:
            subprocess.run(argv, cwd=str(cwd), check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"fastp exited with code {e.returncode}: {shlex.join(argv)}") from e


def _nonempty(path: Path, label: str, report_json: Path | None = None) -> str:
    if not path.is_file():
        raise RuntimeError(f"fastp did not create output file: {path}")
    if path.stat().st_size == 0:
        # Check if 0-read output is legitimate due to filtering
        if report_json and report_json.is_file():
            try:
                data = json.loads(report_json.read_text(encoding="utf-8"))
                passed_reads = data.get("summary", {}).get("after_filtering", {}).get("total_reads", None)
                if passed_reads == 0:
                    return str(path)
            except Exception:
                pass
        raise RuntimeError(f"fastp did not create a nonempty {label}: {path}")
    return str(path)


class Fastp:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Read Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = ("read1", "read2", "json_report", "html_report")

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        if compute_input_fingerprint is not None:
            return compute_input_fingerprint(**kwargs)
        return float("NaN")


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
                "n_base_limit": ("INT", {"default": 5, "min": 0, "max": 100}),
                "length_required": (
                    "INT",
                    {"default": 15, "min": 0, "max": 1000},
                ),
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
        samples = discover_samples(read1, read2)

        executable = None
        if resolve_tool_environment is not None:
            try:
                _, executable = resolve_tool_environment("fastp")
            except Exception:
                executable = shutil.which("fastp")
        else:
            executable = shutil.which("fastp")

        if not executable:
            raise RuntimeError(
                "fastp executable not found in any conda environment or PATH; install conda package fastp=1.3.6"
            )

        base_out = _output_dir("Fastp", output_dir)
        base_out.mkdir(parents=True, exist_ok=True)

        is_single = len(samples) == 1 and Path(str(read1).strip()).is_file()
        out1_list, out2_list, json_list, html_list = [], [], [], []

        managed_flags = {
            "-i", "-I", "-o", "-O", "-w", "--thread", "-q", "-u", "-n", "-l",
            "-j", "-h", "--detect_adapter_for_pe", "--correction"
        }

        for sample_id, fwd_path, rev_path in samples:
            sample_out = base_out if is_single else (base_out / sample_id)
            sample_out.mkdir(parents=True, exist_ok=True)

            out1 = sample_out / (f"{sample_id}_R1.fastq.gz" if not is_single else "R1.fastq.gz")
            out2 = sample_out / (f"{sample_id}_R2.fastq.gz" if not is_single else "R2.fastq.gz")
            report_json = sample_out / "fastp.json"
            report_html = sample_out / "fastp.html"

            argv = [
                executable,
                "-i", str(fwd_path),
                "-o", str(out1),
                "-w", str(threads),
                "-q", str(qualified_quality_phred),
                "-u", str(unqualified_percent_limit),
                "-n", str(n_base_limit),
                "-l", str(length_required),
                "-j", str(report_json),
                "-h", str(report_html),
            ]
            if rev_path:
                argv += ["-I", str(rev_path), "-O", str(out2)]
            if detect_adapter_for_pe:
                argv.append("--detect_adapter_for_pe")
            if correction:
                argv.append("--correction")

            inputs = [fwd_path]
            if rev_path:
                inputs.append(rev_path)
            outputs = [out1, report_json, report_html]
            if rev_path:
                outputs.append(out2)

            try:
                _run(
                    argv=argv,
                    cwd=sample_out,
                    inputs=inputs,
                    outputs=outputs,
                    extra_command=extra_command,
                    managed_flags=managed_flags,
                )
            except TypeError:
                _run(argv, sample_out)

            out1_list.append(_nonempty(out1, f"read 1 output for {sample_id}", report_json))
            if rev_path:
                out2_list.append(_nonempty(out2, f"read 2 output for {sample_id}", report_json))
            json_list.append(_nonempty(report_json, f"JSON report for {sample_id}"))
            html_list.append(_nonempty(report_html, f"HTML report for {sample_id}"))

        if is_single:
            return out1_list[0], (out2_list[0] if out2_list else ""), json_list[0], html_list[0]
        return ",".join(out1_list), ",".join(out2_list), ",".join(json_list), ",".join(html_list)


NODE_CLASS_MAPPINGS = {"Fastp": Fastp}
NODE_DISPLAY_NAME_MAPPINGS = {"Fastp": "fastp: Trim and QC Reads"}


# Backward compatibility aliases
FastpNode = Fastp

__all__ = ["Fastp",
    "FastpNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
