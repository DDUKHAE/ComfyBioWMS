import pytest
from pathlib import Path
import sys
import pandas as pd

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.ref_nodes import (
    SampleMetadataValidatorNode,
    FastpTrimNode,
    SalmonIndexNode,
    SalmonQuantNode,
    TximportNode,
    DESeq2AnalysisNode,
    DESeq2VisualizationNode,
    ComfyBIOReportNode,
)

SEQC_DATA_DIR = root_dir / "paper_platform" / "benchmarks" / "data" / "rnaseq_seqc"
WORKFLOW_OUT = root_dir / "results" / "workflow_seqc_benchmark"

@pytest.fixture(scope="module")
def setup_seqc():
    WORKFLOW_OUT.mkdir(parents=True, exist_ok=True)
    assert SEQC_DATA_DIR.exists(), f"SEQC dataset directory {SEQC_DATA_DIR} must exist"
    
    # Create SEQC sample metadata CSV (Condition A vs Condition B)
    metadata_path = SEQC_DATA_DIR / "sample_metadata.csv"
    metadata_csv = """sample_id,fastq_1,fastq_2,condition
SAMPLE_A1,SAMPLE_A1_R1.fastq.gz,SAMPLE_A1_R2.fastq.gz,UHRR
SAMPLE_B1,SAMPLE_B1_R1.fastq.gz,SAMPLE_B1_R2.fastq.gz,HBRR
"""
    metadata_path.write_text(metadata_csv, encoding="utf-8")
    return metadata_path


def test_seqc_end_to_end_gold_standard_workflow(setup_seqc):
    metadata_path = setup_seqc
    print("\n" + "=" * 70)
    print("🚀 Running SEQC Bulk RNA-Seq Gold-Standard Workflow (End-to-End)")
    print("=" * 70)

    # 1. Validation Node
    val_node = SampleMetadataValidatorNode()
    val_res = val_node.run(fastq_dir=str(SEQC_DATA_DIR), metadata_csv=str(metadata_path))
    assert Path(val_res[0]).exists()
    print("  [Step 1/7] Metadata validated.")

    # 2. Fastp Trimming Node
    trim_out = WORKFLOW_OUT / "trimmed"
    trim_node = FastpTrimNode()
    trim_res = trim_node.run(
        fastp_qc_json="",
        fastq_dir=str(SEQC_DATA_DIR),
        metadata_csv=str(metadata_path),
        output_dir=str(trim_out),
        threads=2,
    )
    assert Path(trim_res[0]).exists()
    print("  [Step 2/7] Fastp adapter & quality trimming completed.")

    # 3. Salmon Index Node
    idx_out = WORKFLOW_OUT / "salmon_index"
    idx_node = SalmonIndexNode()
    idx_res = idx_node.run(
        transcriptome_fasta_path="",
        transcriptome_fasta=str(SEQC_DATA_DIR / "reference_transcriptome.fasta"),
        index_dir=str(idx_out),
        threads=2,
        extra_command="-k 15",
    )
    assert Path(idx_res[0]).exists()
    print("  [Step 3/7] Salmon transcriptome indexing completed.")

    # 4. Salmon Quant Node
    quant_out = WORKFLOW_OUT / "salmon_quant"
    quant_node = SalmonQuantNode()
    quant_res = quant_node.run(
        salmon_index_dir="",
        index_dir=str(idx_out),
        fastq_dir=str(SEQC_DATA_DIR),
        metadata_csv=str(metadata_path),
        trimmed_dir=str(trim_out),
        output_dir=str(quant_out),
        read_layout="A",
        threads=2,
    )
    assert Path(quant_res[0]).exists()
    assert (quant_out / "SAMPLE_A1" / "quant.sf").exists()
    assert (quant_out / "SAMPLE_B1" / "quant.sf").exists()
    print("  [Step 4/7] Salmon quantification completed for all replicates.")

    # 5. Tximport Node
    counts_out = WORKFLOW_OUT / "deseq2" / "count_matrix.csv"
    tximport_node = TximportNode()
    tximport_res = tximport_node.run(
        salmon_quant_dir_path="",
        salmon_quant_dir=str(quant_out),
        metadata_csv=str(metadata_path),
        output_count_matrix=str(counts_out),
    )
    assert Path(tximport_res[0]).exists()
    df_counts = pd.read_csv(counts_out)
    assert len(df_counts) > 0
    print(f"  [Step 5/7] Tximport generated gene count matrix ({len(df_counts)} genes).")

    # 6. DESeq2 Analysis Node
    deseq2_results_out = WORKFLOW_OUT / "deseq2" / "deseq2_results.csv"
    deseq2_node = DESeq2AnalysisNode()
    deseq2_res = deseq2_node.run(
        deseq2_count_matrix="",
        count_matrix=str(counts_out),
        sample_metadata=str(metadata_path),
        results_csv=str(deseq2_results_out),
    )
    assert Path(deseq2_res[0]).exists()
    df_deg = pd.read_csv(deseq2_results_out)
    assert len(df_deg) > 0
    print(f"  [Step 6/7] DESeq2 DEG analysis completed ({len(df_deg)} tested genes).")

    # 7. DESeq2 Visualization & Reporting Node
    plots_out = WORKFLOW_OUT / "plots"
    viz_node = DESeq2VisualizationNode()
    viz_res = viz_node.run(
        deseq2_results_table="",
        count_matrix=str(counts_out),
        results_csv=str(deseq2_results_out),
        plot_dir=str(plots_out),
    )
    assert Path(viz_res[0]).exists()
    print(f"  [Step 7/7] Visualization plots generated in {plots_out}.")

    # 8. Report Node
    report_out = WORKFLOW_OUT / "report" / "SEQC_Analysis_Report.md"
    report_node = ComfyBIOReportNode()
    report_res = report_node.run(
        plot_dir_path="",
        results_csv=str(deseq2_results_out),
        plot_dir=str(plots_out),
        report_path=str(report_out),
    )
    assert Path(report_res[0]).exists()
    assert report_out.read_text(encoding="utf-8").startswith("# ComfyBIO")
    print(f"  [PASS] Summary Markdown report written to {report_out}.")
    print("=" * 70)
