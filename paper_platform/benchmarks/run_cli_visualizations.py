"""
Execute standard visualizations directly for Native CLI pipeline results across all 5 biological domains.
"""

import os
import subprocess
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
CLI_DIR = _REPO_ROOT / "paper_platform" / "benchmarks" / "results_cli"
SCRIPTS_DIR = _REPO_ROOT / "engine" / "scripts"


def run_conda_cmd(env_name: str, cmd_args: list[str], cwd: Path) -> subprocess.CompletedProcess:
    full_cmd = ["conda", "run", "-n", env_name] + cmd_args
    res = subprocess.run(
        full_cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"Visualization command failed ({res.returncode}): {' '.join(full_cmd)}\nStderr: {res.stderr}\nStdout: {res.stdout}")
    return res


def visualize_all_cli_results():
    print("=" * 80)
    print("▶ [CLI VISUALIZATION] Generating Publication Plots for All 5 Native CLI Pipelines")
    print("=" * 80)

    generated_plots = []

    # 1. PhiX174 Assembly (QUAST -> assembly_summary.png)
    print("\n[1/5] Visualizing PhiX174 Assembly (QUAST)...")
    qc_dir = CLI_DIR / "phix174_assembly_cli" / "quast"
    out_plot = CLI_DIR / "phix174_assembly_cli" / "plots" / "assembly_summary.png"
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    run_conda_cmd("genome_assembly", [
        "python", str(SCRIPTS_DIR / "assembly_visualization.py"),
        "--qc-dir", str(qc_dir),
        "--output", str(out_plot)
    ], CLI_DIR)
    assert out_plot.exists()
    print(f"  ✓ Saved PhiX174 Assembly Plot: {out_plot} ({out_plot.stat().st_size:,} bytes)")
    generated_plots.append(("PhiX174 Assembly", str(out_plot)))

    # 2. RNA-Seq SEQC (Counts + DEG -> pca.png, ma.png)
    print("\n[2/5] Visualizing FDA SEQC RNA-Seq (DESeq2)...")
    counts_csv = CLI_DIR / "rnaseq_seqc_cli" / "counts_matrix.csv"
    deg_csv = CLI_DIR / "rnaseq_seqc_cli" / "deseq2_deg_results.csv"
    plot_dir_rna = CLI_DIR / "rnaseq_seqc_cli" / "plots"
    plot_dir_rna.mkdir(parents=True, exist_ok=True)
    run_conda_cmd("bulk_rna_seq", [
        "Rscript", str(SCRIPTS_DIR / "deseq2_visualization.R"),
        str(counts_csv), str(deg_csv), str(plot_dir_rna)
    ], CLI_DIR)
    pca_plot = plot_dir_rna / "pca.png"
    ma_plot = plot_dir_rna / "ma.png"
    assert pca_plot.exists() and ma_plot.exists()
    print(f"  ✓ Saved RNA-Seq Plots: {pca_plot} & {ma_plot}")
    generated_plots.append(("RNA-Seq SEQC", f"{pca_plot}, {ma_plot}"))

    # 3. Zymo Metagenomics (Kraken2/Bracken -> taxa_bar_chart.png)
    print("\n[3/5] Visualizing Zymo Metagenomics (Kraken2/Bracken)...")
    k_dir = CLI_DIR / "metagenome_zymo_cli" / "kraken2"
    out_plot_meta = CLI_DIR / "metagenome_zymo_cli" / "plots" / "taxa_bar_chart.png"
    out_plot_meta.parent.mkdir(parents=True, exist_ok=True)
    run_conda_cmd("metagenome", [
        "python", str(SCRIPTS_DIR / "metagenome_visualization.py"),
        "--reports-dir", str(k_dir),
        "--output", str(out_plot_meta),
        "--top-n", "8"
    ], CLI_DIR)
    assert out_plot_meta.exists()
    print(f"  ✓ Saved Metagenome Taxa Plot: {out_plot_meta} ({out_plot_meta.stat().st_size:,} bytes)")
    generated_plots.append(("Zymo Metagenomics", str(out_plot_meta)))

    # 4. ENCODE ATAC-Seq (MACS3 NarrowPeak -> atac_summary.png)
    print("\n[4/5] Visualizing ENCODE ATAC-Seq (MACS3 Peaks)...")
    peaks_dir = CLI_DIR / "atacseq_encode_cli" / "macs3_peaks"
    out_plot_atac = CLI_DIR / "atacseq_encode_cli" / "plots" / "atac_summary.png"
    out_plot_atac.parent.mkdir(parents=True, exist_ok=True)
    run_conda_cmd("epigenomics", [
        "python", str(SCRIPTS_DIR / "atac_peak_visualization.py"),
        "--peaks-dir", str(peaks_dir),
        "--output", str(out_plot_atac)
    ], CLI_DIR)
    assert out_plot_atac.exists()
    print(f"  ✓ Saved ATAC-Seq Peak Plot: {out_plot_atac} ({out_plot_atac.stat().st_size:,} bytes)")
    generated_plots.append(("ENCODE ATAC-Seq", str(out_plot_atac)))

    # 5. GIAB Variant Calling (BCFtools stats -> variant_summary.png)
    print("\n[5/5] Visualizing GIAB Variant Calling (BCFtools Stats)...")
    vcf_dir = CLI_DIR / "variant_giab_cli" / "filtered_vcf"
    stats_dir = CLI_DIR / "variant_giab_cli" / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    for sample_dir in vcf_dir.iterdir():
        if sample_dir.is_dir():
            filt_vcf = sample_dir / "filtered.vcf"
            stat_file = stats_dir / f"{sample_dir.name}.stats.txt"
            res = run_conda_cmd("variant_analysis", [
                "bcftools", "stats", str(filt_vcf)
            ], CLI_DIR)
            stat_file.write_text(res.stdout, encoding="utf-8")

    out_plot_var = CLI_DIR / "variant_giab_cli" / "plots" / "variant_summary.png"
    out_plot_var.parent.mkdir(parents=True, exist_ok=True)
    run_conda_cmd("variant_analysis", [
        "python", str(SCRIPTS_DIR / "variant_visualization.py"),
        "--stats-dir", str(stats_dir),
        "--output", str(out_plot_var)
    ], CLI_DIR)
    assert out_plot_var.exists()
    print(f"  ✓ Saved Variant Summary Plot: {out_plot_var} ({out_plot_var.stat().st_size:,} bytes)")
    generated_plots.append(("GIAB Variant Calling", str(out_plot_var)))

    print("\n" + "=" * 80)
    print("✅ All 5 Native CLI Pipeline Visualizations Successfully Generated!")
    print("=" * 80)
    return generated_plots


if __name__ == "__main__":
    visualize_all_cli_results()
