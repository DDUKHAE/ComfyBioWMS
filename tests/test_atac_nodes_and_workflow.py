import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.atac_nodes import (
    AtacInputValidatorNode,
    AtacFastpTrimNode,
    AtacBwaMem2IndexNode,
    AtacBwaMem2AlignNode,
    AtacMarkDuplicatesNode,
    AtacQualityFilterNode,
    Macs3PeakCallingNode,
    AtacPeakVisualizationNode,
    AtacReportNode,
)

DATA_DIR = root_dir / "data" / "nf_core_atacseq"
OUTPUT_BASE = root_dir / "results" / "test_atacseq"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    assert DATA_DIR.exists(), f"ATAC data directory {DATA_DIR} must exist"
    assert (DATA_DIR / "sample_metadata.csv").exists()
    assert (DATA_DIR / "genome.fasta").exists()


def test_01_atac_input_validator_node(setup_dirs):
    """Test AtacInputValidatorNode on real dataset."""
    node = AtacInputValidatorNode()
    result = node.run(
        fastq_dir=str(DATA_DIR),
        reference_fasta=str(DATA_DIR / "genome.fasta"),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
    )
    assert len(result) == 1
    assert Path(result[0]).exists()
    print("\n[PASS] AtacInputValidatorNode validated inputs successfully.")


def test_02_atac_fastp_trim_node(setup_dirs):
    """Test AtacFastpTrimNode trimming paired-end ATAC reads."""
    trim_out = OUTPUT_BASE / "trimmed"
    node = AtacFastpTrimNode()
    result = node.run(
        sample_metadata_csv="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(trim_out),
        threads=2,
    )
    assert len(result) == 1
    assert (trim_out / "ATAC_REP1" / "R1.fastq").exists()
    print(f"\n[PASS] AtacFastpTrimNode trimmed reads in {trim_out}.")


def test_03_atac_bwa_mem2_index_and_align_nodes(setup_dirs):
    """Test indexing and alignment for ATAC-seq."""
    trim_out = OUTPUT_BASE / "trimmed"
    align_out = OUTPUT_BASE / "aligned"
    
    # 1. Index
    idx_node = AtacBwaMem2IndexNode()
    idx_res = idx_node.run(
        trimmed_fastq_dir="",
        reference_fasta=str(DATA_DIR / "genome.fasta"),
    )
    assert len(idx_res) == 1
    
    # 2. Align
    align_node = AtacBwaMem2AlignNode()
    align_res = align_node.run(
        reference_fasta_indexed="",
        fastq_dir=str(DATA_DIR),
        reference_fasta=str(DATA_DIR / "genome.fasta"),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        trimmed_dir=str(trim_out),
        output_dir=str(align_out),
        threads=2,
    )
    assert len(align_res) == 1
    assert (align_out / "ATAC_REP1" / "sorted.bam").exists()
    assert (align_out / "ATAC_REP1" / "sorted.bam.bai").exists()
    print(f"\n[PASS] AtacBwaMem2AlignNode generated sorted.bam in {align_out}.")


def test_04_atac_mark_duplicates_and_filter_nodes(setup_dirs):
    """Test duplicate marking and quality filtering."""
    align_out = OUTPUT_BASE / "aligned"
    dedup_out = OUTPUT_BASE / "dedup"
    filtered_out = OUTPUT_BASE / "filtered_bam"
    
    # 1. Mark duplicates
    dedup_node = AtacMarkDuplicatesNode()
    dedup_res = dedup_node.run(
        sorted_bam_dir="",
        input_dir=str(align_out),
        output_dir=str(dedup_out),
        threads=2,
    )
    assert len(dedup_res) == 1
    assert (dedup_out / "ATAC_REP1" / "dedup.bam").exists()
    
    # 2. Quality filter
    qf_node = AtacQualityFilterNode()
    qf_res = qf_node.run(
        dedup_bam_dir="",
        input_dir=str(dedup_out),
        output_dir=str(filtered_out),
        min_mapq=10,
    )
    assert len(qf_res) == 1
    assert (filtered_out / "ATAC_REP1" / "final.bam").exists()
    print(f"\n[PASS] AtacQualityFilterNode generated final.bam in {filtered_out}.")


def test_05_macs3_peak_calling_node(setup_dirs):
    """Test MACS3 peak calling for ATAC-seq chromatin accessibility."""
    filtered_out = OUTPUT_BASE / "filtered_bam"
    peaks_out = OUTPUT_BASE / "peaks"
    
    node = Macs3PeakCallingNode()
    result = node.run(
        filtered_bam_dir="",
        input_dir=str(filtered_out),
        output_dir=str(peaks_out),
        genome_size="40000",
    )
    assert len(result) == 1
    assert (peaks_out / "ATAC_REP1" / "ATAC_REP1_peaks.narrowPeak").exists()
    print(f"\n[PASS] Macs3PeakCallingNode identified narrowPeak regions in {peaks_out}.")


def test_06_atac_visualization_and_report_nodes(setup_dirs):
    """Test peak profile visualization and summary report generation."""
    peaks_out = OUTPUT_BASE / "peaks"
    plots_out = OUTPUT_BASE / "plots"
    report_out = OUTPUT_BASE / "report" / "atac_report.md"
    
    # 1. Visualization
    viz_node = AtacPeakVisualizationNode()
    viz_res = viz_node.run(
        peaks_dir="",
        input_dir=str(peaks_out),
        plot_dir=str(plots_out),
    )
    assert len(viz_res) == 2
    assert (plots_out / "atac_summary.png").exists()
    print(f"\n[PASS] AtacPeakVisualizationNode rendered peak summary plot in {plots_out}.")
    
    # 2. Report
    report_node = AtacReportNode()
    report_res = report_node.run(
        plot_dir_path="",
        peaks_dir=str(peaks_out),
        plot_dir=str(plots_out),
        report_path=str(report_out),
    )
    assert len(report_res) == 1
    assert Path(report_res[0]).exists()
    assert report_out.read_text(encoding="utf-8").startswith("# ComfyBIO")
    print(f"\n[PASS] AtacReportNode generated summary report at {report_out}.")
