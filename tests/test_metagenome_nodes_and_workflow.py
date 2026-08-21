import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.metagenome_nodes import (
    MetagenomeInputValidatorNode,
    MetagenomeFastpTrimNode,
    Kraken2ClassifyNode,
    BrackenAbundanceNode,
    MetagenomeVisualizationNode,
    MetagenomeReportNode,
)

DATA_DIR = root_dir / "paper_platform" / "benchmarks" / "data" / "metagenome_zymo"
DB_DIR = DATA_DIR / "kraken2_db" / "testdb-kraken2"
OUTPUT_BASE = root_dir / "results" / "test_metagenome"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    assert DATA_DIR.exists(), f"Metagenome data directory {DATA_DIR} must exist"
    assert DB_DIR.exists(), f"Kraken2 DB directory {DB_DIR} must exist"
    
    # Ensure sample metadata CSV exists
    metadata_content = """sample_id,fastq_1,fastq_2,condition
ZYMO_MOCK1,zymo_mock_R1.fastq.gz,zymo_mock_R2.fastq.gz,mock_community
"""
    (DATA_DIR / "sample_metadata.csv").write_text(metadata_content, encoding="utf-8")


def test_01_metagenome_input_validator_node(setup_dirs):
    """Test MetagenomeInputValidatorNode on real Zymo dataset and Kraken2 DB."""
    node = MetagenomeInputValidatorNode()
    result = node.run(
        fastq_dir=str(DATA_DIR),
        kraken2_db_dir=str(DB_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
    )
    assert len(result) == 1
    assert Path(result[0]).exists()
    print("\n[PASS] MetagenomeInputValidatorNode validated inputs successfully.")


def test_02_metagenome_fastp_trim_node(setup_dirs):
    """Test MetagenomeFastpTrimNode trimming paired-end metagenomic reads."""
    trim_out = OUTPUT_BASE / "trimmed"
    node = MetagenomeFastpTrimNode()
    result = node.run(
        sample_metadata_csv="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(trim_out),
        threads=2,
    )
    assert len(result) == 1
    assert (trim_out / "ZYMO_MOCK1" / "R1.fastq").exists()
    print(f"\n[PASS] MetagenomeFastpTrimNode trimmed reads in {trim_out}.")


def test_03_kraken2_classify_node(setup_dirs):
    """Test Kraken2ClassifyNode assigning taxonomy to metagenomic reads."""
    trim_out = OUTPUT_BASE / "trimmed"
    kraken_out = OUTPUT_BASE / "kraken2"
    
    node = Kraken2ClassifyNode()
    result = node.run(
        trimmed_fastq_dir="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        trimmed_dir=str(trim_out),
        kraken2_db_dir=str(DB_DIR),
        output_dir=str(kraken_out),
        threads=2,
        confidence=0.0,
    )
    assert len(result) == 1
    report_file = kraken_out / "ZYMO_MOCK1" / "kraken2_report.txt"
    assert report_file.exists() and report_file.stat().st_size > 0
    print(f"\n[PASS] Kraken2ClassifyNode generated taxonomic report at {report_file}.")


def test_04_bracken_abundance_node(setup_dirs):
    """Test BrackenAbundanceNode re-estimating Bayesian taxonomic abundances."""
    kraken_out = OUTPUT_BASE / "kraken2"
    bracken_out = OUTPUT_BASE / "bracken"
    
    node = BrackenAbundanceNode()
    result = node.run(
        kraken2_output_dir="",
        input_dir=str(kraken_out),
        kraken2_db_dir=str(DB_DIR),
        output_dir=str(bracken_out),
        read_length=100,
        level="S",
        threshold=1,
    )
    assert len(result) == 1
    bracken_report = bracken_out / "ZYMO_MOCK1" / "bracken_report.txt"
    assert bracken_report.exists() and bracken_report.stat().st_size > 0
    print(f"\n[PASS] BrackenAbundanceNode estimated species abundances in {bracken_out}.")


def test_05_metagenome_visualization_and_report_nodes(setup_dirs):
    """Test metagenomic taxonomy visualization and summary report generation."""
    bracken_out = OUTPUT_BASE / "bracken"
    plots_out = OUTPUT_BASE / "plots"
    report_out = OUTPUT_BASE / "report" / "metagenome_report.md"
    
    # 1. Visualization
    viz_node = MetagenomeVisualizationNode()
    viz_res = viz_node.run(
        bracken_dir="",
        input_dir=str(bracken_out),
        plot_dir=str(plots_out),
    )
    assert len(viz_res) == 2
    assert (plots_out / "metagenome_summary.png").exists()
    print(f"\n[PASS] MetagenomeVisualizationNode rendered taxonomic composition plot in {plots_out}.")
    
    # 2. Report
    report_node = MetagenomeReportNode()
    report_res = report_node.run(
        plot_dir_path="",
        bracken_dir=str(bracken_out),
        plot_dir=str(plots_out),
        report_path=str(report_out),
    )
    assert len(report_res) == 1
    assert Path(report_res[0]).exists()
    assert report_out.read_text(encoding="utf-8").startswith("# ComfyBIO")
    print(f"\n[PASS] MetagenomeReportNode generated summary report at {report_out}.")
