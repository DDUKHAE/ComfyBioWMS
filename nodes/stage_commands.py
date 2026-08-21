import sys
from pathlib import Path

from bioflow.runtime.command_runner import conda_command, parse_extra_command_tokens
from bioflow.runtime.command_runner import fastp_trim_argv as _fastp_trim_argv

ENV_NAME = "bulk_rna_seq"
_REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = _REPO_ROOT / "engine" / "scripts"
REPORT_SCRIPT = _REPO_ROOT / "engine" / "src" / "bioflow" / "runtime" / "report.py"


def fastp_qc_argv(sample, output_dir, threads, extra_command="") -> list[str]:
    out = Path(output_dir)
    args = ["-i", str(sample.fastq_1)]
    if sample.fastq_2 is not None:
        args += ["-I", str(sample.fastq_2)]
    args += [
        "-w", str(threads),
        "--json", str(out / f"{sample.sample_id}.fastp.json"),
        "--html", str(out / f"{sample.sample_id}.fastp.html"),
    ]
    return conda_command(ENV_NAME, "fastp", *args, *parse_extra_command_tokens(extra_command))


def fastp_trim_argv(sample, sample_output_dir, threads, extra_command="") -> list[str]:
    return _fastp_trim_argv(ENV_NAME, sample, sample_output_dir, threads, extra_command)


def fastqc_qc_argv(sample, output_dir, threads, extra_command="") -> list[str]:
    out = Path(output_dir)
    args = ["-o", str(out), "-t", str(threads), str(sample.fastq_1)]
    if sample.fastq_2 is not None:
        args.append(str(sample.fastq_2))
    return conda_command(ENV_NAME, "fastqc", *args, *parse_extra_command_tokens(extra_command))


_TRIMMOMATIC_DEFAULT_ADAPTER = "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10"


def trimmomatic_trim_argv(sample, sample_output_dir, threads, extra_command="") -> list[str]:
    out = Path(sample_output_dir)
    extra_tokens = parse_extra_command_tokens(extra_command)
    has_custom_adapter = any(token.startswith("ILLUMINACLIP:") for token in extra_tokens)
    if sample.fastq_2 is not None:
        args = [
            "PE", "-threads", str(threads),
            str(sample.fastq_1), str(sample.fastq_2),
            str(out / "R1.fastq"), str(out / "R1.unpaired.fastq"),
            str(out / "R2.fastq"), str(out / "R2.unpaired.fastq"),
        ]
    else:
        args = ["SE", "-threads", str(threads), str(sample.fastq_1), str(out / "R1.fastq")]
    if not has_custom_adapter:
        args.append(_TRIMMOMATIC_DEFAULT_ADAPTER)
    return conda_command(ENV_NAME, "trimmomatic", *args, *extra_tokens)


def salmon_index_argv(transcriptome_fasta, index_dir, threads, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "salmon", "index",
        "-t", str(transcriptome_fasta),
        "-i", str(index_dir),
        "-p", str(threads),
        *parse_extra_command_tokens(extra_command),
    )


def salmon_quant_argv(index_dir, read1, read2, output_dir, read_layout, threads, extra_command="") -> list[str]:
    args = ["-i", str(index_dir), "-l", str(read_layout), "-1", str(read1)]
    if read2 is not None:
        args += ["-2", str(read2)]
    args += ["-p", str(threads), "-o", str(output_dir)]
    return conda_command(ENV_NAME, "salmon", "quant", *args, *parse_extra_command_tokens(extra_command))


def tximport_argv(salmon_quant_dir, count_matrix, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "Rscript", str(SCRIPT_DIR / "tximport_import.R"),
        str(salmon_quant_dir), str(count_matrix),
        *parse_extra_command_tokens(extra_command),
    )


def deseq2_argv(count_matrix, sample_metadata, results_csv, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "Rscript", str(SCRIPT_DIR / "deseq2_analysis.R"),
        str(count_matrix), str(sample_metadata), str(results_csv),
        *parse_extra_command_tokens(extra_command),
    )


def deseq2_viz_argv(count_matrix, results_csv, plot_dir, extra_command="") -> list[str]:
    return conda_command(
        ENV_NAME, "Rscript", str(SCRIPT_DIR / "deseq2_visualization.R"),
        str(count_matrix), str(results_csv), str(plot_dir),
        *parse_extra_command_tokens(extra_command),
    )


def report_argv(results_csv, plot_dir, report_path) -> list[str]:
    return [
        sys.executable, str(REPORT_SCRIPT),
        "--results", str(results_csv),
        "--plot-dir", str(plot_dir),
        "--output", str(report_path),
    ]
