import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.variant_nodes import (
    VariantInputValidatorNode,
    BwaMem2IndexNode,
    BwaMem2AlignNode,
    MarkDuplicatesNode,
    BcftoolsCallNode,
    BcftoolsFilterNode,
    VariantVisualizationNode,
    VariantReportNode,
)

DATA_DIR = root_dir / "data" / "nf_core_variant"
OUTPUT_BASE = root_dir / "results" / "test_variant_calling"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    assert DATA_DIR.exists(), f"Variant data directory {DATA_DIR} must exist"
    
    # Ensure sample metadata CSV has condition column
    metadata_content = """sample_id,fastq_1,fastq_2,condition
SAMPLE_VAR1,sample1_R1.fastq.gz,sample1_R2.fastq.gz,tumor
"""
    (DATA_DIR / "sample_metadata.csv").write_text(metadata_content, encoding="utf-8")
    assert (DATA_DIR / "genome.fasta").exists()


def test_01_variant_input_validator_node(setup_dirs):
    """Test VariantInputValidatorNode on real nf-core dataset."""
    node = VariantInputValidatorNode()
    result = node.run(
        fastq_dir=str(DATA_DIR),
        reference_fasta=str(DATA_DIR / "genome.fasta"),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
    )
    assert len(result) == 1
    assert Path(result[0]).exists()
    print("\n[PASS] VariantInputValidatorNode validated inputs successfully.")


def test_02_bwa_mem2_index_node(setup_dirs):
    """Test BwaMem2IndexNode indexing reference genome FASTA."""
    node = BwaMem2IndexNode()
    result = node.run(
        sample_metadata_csv="",
        reference_fasta=str(DATA_DIR / "genome.fasta"),
    )
    assert len(result) == 1
    ref_path = Path(result[0])
    assert ref_path.exists()
    print(f"\n[PASS] BwaMem2IndexNode indexed {ref_path}.")


def test_03_bwa_mem2_align_node(setup_dirs):
    """Test BwaMem2AlignNode aligning paired-end reads to indexed reference."""
    align_out = OUTPUT_BASE / "aligned"
    node = BwaMem2AlignNode()
    result = node.run(
        reference_fasta_indexed="",
        fastq_dir=str(DATA_DIR),
        reference_fasta=str(DATA_DIR / "genome.fasta"),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(align_out),
        threads=2,
    )
    assert len(result) == 1
    assert (align_out / "SAMPLE_VAR1" / "sorted.bam").exists()
    assert (align_out / "SAMPLE_VAR1" / "sorted.bam.bai").exists()
    print(f"\n[PASS] BwaMem2AlignNode aligned reads to BAM in {align_out}.")


def test_04_mark_duplicates_node(setup_dirs):
    """Test MarkDuplicatesNode using samtools markdup."""
    align_out = OUTPUT_BASE / "aligned"
    dedup_out = OUTPUT_BASE / "dedup"
    node = MarkDuplicatesNode()
    result = node.run(
        sorted_bam_dir="",
        input_dir=str(align_out),
        output_dir=str(dedup_out),
        threads=2,
    )
    assert len(result) == 1
    assert (dedup_out / "SAMPLE_VAR1" / "dedup.bam").exists()
    assert (dedup_out / "SAMPLE_VAR1" / "dedup.bam.bai").exists()
    print(f"\n[PASS] MarkDuplicatesNode deduplicated BAM in {dedup_out}.")


def test_05_bcftools_call_node(setup_dirs):
    """Test BcftoolsCallNode generating VCF calls from deduplicated BAM."""
    dedup_out = OUTPUT_BASE / "dedup"
    calls_out = OUTPUT_BASE / "calls"
    node = BcftoolsCallNode()
    result = node.run(
        dedup_bam_dir="",
        input_dir=str(dedup_out),
        reference_fasta=str(DATA_DIR / "genome.fasta"),
        output_dir=str(calls_out),
    )
    assert len(result) == 1
    assert (calls_out / "SAMPLE_VAR1" / "raw.vcf").exists()
    print(f"\n[PASS] BcftoolsCallNode generated raw.vcf in {calls_out}.")


def test_06_bcftools_filter_node(setup_dirs):
    """Test BcftoolsFilterNode filtering VCF calls."""
    calls_out = OUTPUT_BASE / "calls"
    filter_out = OUTPUT_BASE / "filtered"
    node = BcftoolsFilterNode()
    result = node.run(
        raw_vcf_dir="",
        input_dir=str(calls_out),
        output_dir=str(filter_out),
        exclude_expression="QUAL<10",
    )
    assert len(result) == 1
    assert (filter_out / "SAMPLE_VAR1" / "filtered.vcf").exists()
    print(f"\n[PASS] BcftoolsFilterNode filtered VCF in {filter_out}.")


def test_07_variant_visualization_and_report_nodes(setup_dirs):
    """Test VariantVisualizationNode and VariantReportNode."""
    filter_out = OUTPUT_BASE / "filtered"
    plots_out = OUTPUT_BASE / "plots"
    report_out = OUTPUT_BASE / "report" / "variant_report.md"

    # Visualization
    viz_node = VariantVisualizationNode()
    viz_res = viz_node.run(
        filtered_vcf_dir="",
        input_dir=str(filter_out),
        plot_dir=str(plots_out),
    )
    assert len(viz_res) == 2
    assert (plots_out / "variant_summary.png").exists()
    print(f"\n[PASS] VariantVisualizationNode created summary plot in {plots_out}.")

    # Report
    report_node = VariantReportNode()
    report_res = report_node.run(
        plot_dir_path="",
        vcf_dir=str(filter_out),
        plot_dir=str(plots_out),
        report_path=str(report_out),
    )
    assert len(report_res) == 1
    assert Path(report_res[0]).exists()
    assert report_out.read_text(encoding="utf-8").startswith("# ComfyBIO")
    print(f"\n[PASS] VariantReportNode generated summary report at {report_out}.")
