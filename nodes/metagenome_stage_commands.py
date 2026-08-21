import sys
from pathlib import Path

from bioflow.runtime.command_runner import conda_command, parse_extra_command_tokens
from bioflow.runtime.command_runner import fastp_trim_argv as _fastp_trim_argv

ENV_NAME = "metagenome"
_REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = _REPO_ROOT / "engine" / "scripts"
REPORT_SCRIPT = _REPO_ROOT / "engine" / "src" / "bioflow" / "runtime" / "metagenome_report.py"


def fastp_trim_argv(sample, sample_output_dir, threads, extra_command="") -> list[str]:
    return _fastp_trim_argv(ENV_NAME, sample, sample_output_dir, threads, extra_command)


def kraken2_classify_argv(db_dir, read1, read2, report_path, output_path, threads, confidence, extra_command="") -> list[str]:
    args = ["--db", str(db_dir), "--threads", str(threads), "--confidence", str(confidence)]
    if read2 is not None:
        args.append("--paired")
    args += ["--report", str(report_path), "--output", str(output_path), str(read1)]
    if read2 is not None:
        args.append(str(read2))
    return conda_command(ENV_NAME, "kraken2", *args, *parse_extra_command_tokens(extra_command))


def bracken_abundance_argv(db_dir, kraken2_report, output_path, report_path, read_length, level, threshold, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "bracken",
        "-d", str(db_dir), "-i", str(kraken2_report),
        "-o", str(output_path), "-w", str(report_path),
        "-r", str(read_length), "-l", str(level), "-t", str(threshold),
        *parse_extra_command_tokens(extra_command),
    )


def metagenome_visualization_argv(reports_dir, plot_dir, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "python", str(SCRIPT_DIR / "metagenome_visualization.py"),
        "--reports-dir", str(reports_dir),
        "--output", str(Path(plot_dir) / "metagenome_summary.png"),
        *parse_extra_command_tokens(extra_command),
    )


def metagenome_report_argv(bracken_dir, plot_dir, report_path) -> list[str]:
    return [
        sys.executable, str(REPORT_SCRIPT),
        "--bracken-dir", str(bracken_dir),
        "--plot-dir", str(plot_dir),
        "--output", str(report_path),
    ]
