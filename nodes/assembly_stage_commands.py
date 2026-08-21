import sys
from pathlib import Path

from bioflow.runtime.command_runner import conda_command, parse_extra_command_tokens
from bioflow.runtime.command_runner import fastp_trim_argv as _fastp_trim_argv

ENV_NAME = "genome_assembly"
_REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = _REPO_ROOT / "engine" / "scripts"
REPORT_SCRIPT = _REPO_ROOT / "engine" / "src" / "bioflow" / "runtime" / "assembly_report.py"


def fastp_trim_argv(sample, sample_output_dir, threads, extra_command="") -> list[str]:
    return _fastp_trim_argv(ENV_NAME, sample, sample_output_dir, threads, extra_command)


def spades_assemble_argv(read1, read2, outdir, threads, memory_gb, extra_command="") -> list[str]:
    args = []
    if read2 is not None:
        args += ["--pe1-1", str(read1), "--pe1-2", str(read2)]
    else:
        args += ["--s1", str(read1)]
    args += ["-o", str(outdir), "--threads", str(threads), "--memory", str(memory_gb), "--isolate"]
    return conda_command(ENV_NAME, "spades.py", *args, *parse_extra_command_tokens(extra_command))


def quast_qc_argv(contigs_fasta, outdir, extra_command="") -> list[str]:
    return conda_command(ENV_NAME, "quast.py", str(contigs_fasta), "-o", str(outdir), *parse_extra_command_tokens(extra_command))


def assembly_visualization_argv(qc_dir, plot_dir, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "python", str(SCRIPT_DIR / "assembly_visualization.py"),
        "--qc-dir", str(qc_dir),
        "--output", str(Path(plot_dir) / "assembly_summary.png"),
        *parse_extra_command_tokens(extra_command),
    )


def assembly_report_argv(qc_dir, plot_dir, report_path) -> list[str]:
    return [
        sys.executable, str(REPORT_SCRIPT),
        "--qc-dir", str(qc_dir),
        "--plot-dir", str(plot_dir),
        "--output", str(report_path),
    ]
