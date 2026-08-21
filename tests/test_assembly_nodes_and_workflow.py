import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.assembly_nodes import (
    AssemblyInputValidatorNode,
    AssemblyFastpTrimNode,
    SpadesAssembleNode,
    QuastQcNode,
    AssemblyVisualizationNode,
    AssemblyReportNode,
)

DATA_DIR = root_dir / "paper_platform" / "benchmarks" / "data" / "phix174"
OUTPUT_BASE = root_dir / "results" / "test_assembly"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    assert DATA_DIR.exists(), f"PhiX174 data directory {DATA_DIR} must exist"
    
    # Ensure sample metadata CSV exists
    metadata_content = """sample_id,fastq_1,fastq_2,condition
PHIX174,phix174_R1.fastq.gz,phix174_R2.fastq.gz,control
"""
    (DATA_DIR / "sample_metadata.csv").write_text(metadata_content, encoding="utf-8")
    assert (DATA_DIR / "phix174_ref.fasta").exists()


def test_01_assembly_input_validator_node(setup_dirs):
    """Test AssemblyInputValidatorNode on real PhiX174 dataset."""
    node = AssemblyInputValidatorNode()
    result = node.run(
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
    )
    assert len(result) == 1
    assert Path(result[0]).exists()
    print("\n[PASS] AssemblyInputValidatorNode validated inputs successfully.")


def test_02_assembly_fastp_trim_node(setup_dirs):
    """Test AssemblyFastpTrimNode trimming paired-end reads."""
    trim_out = OUTPUT_BASE / "trimmed"
    node = AssemblyFastpTrimNode()
    result = node.run(
        sample_metadata_csv="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(trim_out),
        threads=2,
    )
    assert len(result) == 1
    assert (trim_out / "PHIX174" / "R1.fastq").exists()
    print(f"\n[PASS] AssemblyFastpTrimNode trimmed reads in {trim_out}.")


def test_03_spades_assemble_node(setup_dirs):
    """Test SpadesAssembleNode assembling de novo contigs."""
    trim_out = OUTPUT_BASE / "trimmed"
    assembly_out = OUTPUT_BASE / "assembly"
    
    node = SpadesAssembleNode()
    result = node.run(
        trimmed_fastq_dir="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        trimmed_dir=str(trim_out),
        output_dir=str(assembly_out),
        threads=2,
        memory_gb=4,
        extra_command="--sc",  # SPAdes single-cell/small genome mode
    )
    assert len(result) == 1
    contigs_file = assembly_out / "PHIX174" / "contigs.fasta"
    assert contigs_file.exists() and contigs_file.stat().st_size > 0
    print(f"\n[PASS] SpadesAssembleNode assembled contigs: {contigs_file} ({contigs_file.stat().st_size:,} bytes).")


def test_04_quast_qc_node(setup_dirs):
    """Test QuastQcNode assessing assembly quality metrics."""
    assembly_out = OUTPUT_BASE / "assembly"
    quast_out = OUTPUT_BASE / "quast"
    
    node = QuastQcNode()
    result = node.run(
        assembly_dir="",
        input_dir=str(assembly_out),
        output_dir=str(quast_out),
        extra_command=f"-r {str(DATA_DIR / 'phix174_ref.fasta')}",
    )
    assert len(result) == 1
    report_tsv = quast_out / "PHIX174" / "report.tsv"
    assert report_tsv.exists()
    print(f"\n[PASS] QuastQcNode generated report.tsv in {quast_out}.")


def test_05_assembly_visualization_and_report_nodes(setup_dirs):
    """Test assembly summary plot visualization and markdown report."""
    quast_out = OUTPUT_BASE / "quast"
    plots_out = OUTPUT_BASE / "plots"
    report_out = OUTPUT_BASE / "report" / "assembly_report.md"
    
    # 1. Visualization
    viz_node = AssemblyVisualizationNode()
    viz_res = viz_node.run(
        qc_dir="",
        input_dir=str(quast_out),
        plot_dir=str(plots_out),
    )
    assert len(viz_res) == 2
    assert (plots_out / "assembly_summary.png").exists()
    print(f"\n[PASS] AssemblyVisualizationNode rendered assembly plot in {plots_out}.")
    
    # 2. Report
    report_node = AssemblyReportNode()
    report_res = report_node.run(
        plot_dir_path="",
        qc_dir=str(quast_out),
        plot_dir=str(plots_out),
        report_path=str(report_out),
    )
    assert len(report_res) == 1
    assert Path(report_res[0]).exists()
    assert report_out.read_text(encoding="utf-8").startswith("# ComfyBIO")
    print(f"\n[PASS] AssemblyReportNode generated summary report at {report_out}.")
