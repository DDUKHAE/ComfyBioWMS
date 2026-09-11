"""Concordance test suite comparing ComfyBIOWMS node outputs directly with native tool executions.

Marked with @pytest.mark.e2e as it requires real biological tools in conda environments.
"""

import gzip
import json
import os
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
import pytest

from bioflow.runtime.command_runner import resolve_tool_environment
from nodes.class_2.fastp import Fastp
from nodes.class_2.deseq2 import DESeq2
from nodes.class_1.tximport import Tximport


@pytest.mark.e2e
def test_fastp_native_cli_concordance(tmp_path):
    """Verify that Fastp node produces identical results to direct native fastp CLI execution."""
    r1_input = Path("data/nf_core_rnaseq/SRR6357070_1.fastq.gz").resolve()
    r2_input = Path("data/nf_core_rnaseq/SRR6357070_2.fastq.gz").resolve()

    if not r1_input.is_file() or not r2_input.is_file():
        pytest.skip("Official RNA-seq test data missing from data/nf_core_rnaseq/")

    _, fastp_bin = resolve_tool_environment("fastp")

    # 1. Direct Native CLI Execution
    native_dir = tmp_path / "native_fastp"
    native_dir.mkdir()
    native_r1 = native_dir / "native_R1.fastq.gz"
    native_r2 = native_dir / "native_R2.fastq.gz"
    native_json = native_dir / "native.json"
    native_html = native_dir / "native.html"

    native_cmd = [
        fastp_bin,
        "-i", str(r1_input),
        "-I", str(r2_input),
        "-o", str(native_r1),
        "-O", str(native_r2),
        "-w", "2",
        "-q", "15",
        "-u", "40",
        "-n", "5",
        "-l", "20",
        "-j", str(native_json),
        "-h", str(native_html),
    ]
    subprocess.run(native_cmd, check=True, cwd=str(native_dir))

    # 2. ComfyBIOWMS Fastp Node Execution
    node_dir = tmp_path / "node_fastp"
    node = Fastp()
    node_r1, node_r2, node_json, node_html = node.run(
        read1=str(r1_input),
        read2=str(r2_input),
        output_dir=str(node_dir),
        threads=2,
        qualified_quality_phred=15,
        unqualified_percent_limit=40,
        n_base_limit=5,
        length_required=20,
    )

    # 3. Concordance Verification
    # Compare summary statistics in JSON reports
    with open(native_json, "r") as f:
        nat_data = json.load(f)
    with open(node_json, "r") as f:
        nod_data = json.load(f)

    nat_after = nat_data["summary"]["after_filtering"]
    nod_after = nod_data["summary"]["after_filtering"]

    assert nat_after["total_reads"] == nod_after["total_reads"], "Total reads must match exactly"
    assert nat_after["total_bases"] == nod_after["total_bases"], "Total bases must match exactly"
    assert nat_after["q30_rate"] == nod_after["q30_rate"], "Q30 rate must match exactly"

    # Compare uncompressed FASTQ sequence records (ignoring gzip header timestamp differences)
    with gzip.open(native_r1, "rt") as f_nat, gzip.open(node_r1, "rt") as f_nod:
        nat_lines = [line.strip() for line in f_nat]
        nod_lines = [line.strip() for line in f_nod]
    assert nat_lines == nod_lines, "Trimmed R1 reads must match line-by-line"

    with gzip.open(native_r2, "rt") as f_nat, gzip.open(node_r2, "rt") as f_nod:
        nat_lines_r2 = [line.strip() for line in f_nat]
        nod_lines_r2 = [line.strip() for line in f_nod]
    assert nat_lines_r2 == nod_lines_r2, "Trimmed R2 reads must match line-by-line"


@pytest.mark.e2e
def test_tximport_length_scaled_tpm_concordance(tmp_path):
    """Verify that Tximport matches the exact Soneson et al. 2015 lengthScaledTPM calculation."""
    # Create two mock Salmon quant.sf outputs with known numeric values
    s1_dir = tmp_path / "sample_1"
    s2_dir = tmp_path / "sample_2"
    s1_dir.mkdir()
    s2_dir.mkdir()

    # gene1: 2 transcripts (tx1, tx2), gene2: 1 transcript (tx3)
    tx2gene_file = tmp_path / "tx2gene.tsv"
    tx2gene_file.write_text("tx1\tgene1\ntx2\tgene1\ntx3\tgene2\n")

    header = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
    # s1: tx1 (len=1000, tpm=20, reads=100), tx2 (len=2000, tpm=40, reads=400), tx3 (len=1500, tpm=60, reads=300)
    (s1_dir / "quant.sf").write_text(
        header +
        "tx1\t1000\t1000\t20.0\t100.0\n" +
        "tx2\t2000\t2000\t40.0\t400.0\n" +
        "tx3\t1500\t1500\t60.0\t300.0\n"
    )
    # s2: tx1 (len=1000, tpm=10, reads=50), tx2 (len=2000, tpm=80, reads=800), tx3 (len=1500, tpm=30, reads=150)
    (s2_dir / "quant.sf").write_text(
        header +
        "tx1\t1000\t1000\t10.0\t50.0\n" +
        "tx2\t2000\t2000\t80.0\t800.0\n" +
        "tx3\t1500\t1500\t30.0\t150.0\n"
    )

    node = Tximport()
    counts_tsv, tpm_tsv, summary_json = node.run(
        quant_files=f"{s1_dir / 'quant.sf'},{s2_dir / 'quant.sf'}",
        output_dir=str(tmp_path / "tximport_out"),
        tx2gene_tsv=str(tx2gene_file),
        sample_names="S1,S2",
        counts_from_abundance="lengthScaledTPM",
    )

    counts_df = pd.read_csv(counts_tsv, sep="\t", index_col=0)
    tpm_df = pd.read_csv(tpm_tsv, sep="\t", index_col=0)

    # In S1: gene1 TPM = 20.0 + 40.0 = 60.0; gene2 TPM = 60.0
    assert np.isclose(tpm_df.loc["gene1", "S1"], 60.0)
    assert np.isclose(tpm_df.loc["gene2", "S1"], 60.0)

    # Total reads in S1 = 100 + 400 + 300 = 800.
    # Total reads in S2 = 50 + 800 + 150 = 1000.
    # lengthScaledTPM scales gene count such that sum(counts) == total original reads
    assert np.isclose(counts_df["S1"].sum(), 800.0, rtol=1e-5)
    assert np.isclose(counts_df["S2"].sum(), 1000.0, rtol=1e-5)


@pytest.mark.e2e
def test_deseq2_native_rscript_concordance(tmp_path):
    """Verify that DESeq2 node outputs match a direct independent native Rscript invocation."""
    _, rscript_bin = resolve_tool_environment("Rscript")

    # Generate synthetic gene matrix and metadata
    counts_file = tmp_path / "counts.csv"
    meta_file = tmp_path / "meta.csv"

    genes = [f"gene_{i}" for i in range(1, 21)]
    # 4 samples: 2 control, 2 treated
    samples = ["ctrl_1", "ctrl_2", "trt_1", "trt_2"]
    meta_df = pd.DataFrame({
        "sample": samples,
        "condition": ["control", "control", "treated", "treated"],
    })
    meta_df.to_csv(meta_file, index=False)

    np.random.seed(123)
    counts_dict = {"gene_id": genes}
    for s in samples:
        base = 50 + np.random.poisson(100, len(genes))
        if s.startswith("trt"):
            base[:5] *= 4  # First 5 genes are upregulated in treatment
        counts_dict[s] = base.tolist()

    counts_df = pd.DataFrame(counts_dict)
    counts_df.to_csv(counts_file, index=False)

    # 1. Run direct independent Rscript with DESeq2
    native_out = tmp_path / "native_r"
    native_out.mkdir()
    native_script = native_out / "native_deseq2.R"
    native_csv = native_out / "native_results.csv"

    r_code = f"""
suppressPackageStartupMessages(library(DESeq2))
counts <- read.csv('{counts_file}', row.names=1)
meta <- read.csv('{meta_file}', row.names=1)
meta$condition <- factor(meta$condition)
meta$condition <- relevel(meta$condition, ref='control')
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)), colData=meta, design=~condition)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds, contrast=c('condition', 'treated', 'control'))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, '{native_csv}', row.names=FALSE)
"""
    native_script.write_text(r_code)
    subprocess.run([rscript_bin, str(native_script)], check=True, cwd=str(native_out))

    # 2. Run ComfyBIOWMS DESeq2 Node
    node_out = tmp_path / "node_deseq2"
    node = DESeq2()
    res_csv, norm_csv, summary_json = node.run(
        count_matrix_csv=str(counts_file),
        sample_metadata_csv=str(meta_file),
        output_dir=str(node_out),
        condition_col="condition",
        contrast_reference="control",
        contrast_target="treated",
    )

    # 3. Compare Results
    df_nat = pd.read_csv(native_csv).set_index("gene_id")
    df_nod = pd.read_csv(res_csv).set_index("gene_id")

    common_genes = df_nat.index.intersection(df_nod.index)
    assert len(common_genes) == len(genes)

    # Verify log2FoldChange matches
    nat_lfc = df_nat.loc[common_genes, "log2FoldChange"].values
    nod_lfc = df_nod.loc[common_genes, "log2FoldChange"].values
    np.testing.assert_allclose(nod_lfc, nat_lfc, rtol=1e-5, atol=1e-5)

    # Verify p-values match
    nat_p = df_nat.loc[common_genes, "pvalue"].dropna().values
    nod_p = df_nod.loc[common_genes, "pvalue"].dropna().values
    np.testing.assert_allclose(nod_p, nat_p, rtol=1e-5, atol=1e-5)
