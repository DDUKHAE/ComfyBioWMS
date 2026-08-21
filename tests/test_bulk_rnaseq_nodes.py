import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.ref_nodes import (
    SampleMetadataValidatorNode,
    FastpQCNode,
    FastpTrimNode,
    FastQCNode,
    SalmonIndexNode,
    SalmonQuantNode,
    TximportNode,
    DESeq2AnalysisNode,
    DESeq2VisualizationNode,
    ComfyBIOReportNode,
)

DATA_DIR = root_dir / "data" / "nf_core_rnaseq"
OUTPUT_BASE = root_dir / "results" / "test_bulk_rnaseq"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    assert DATA_DIR.exists(), f"Test data directory {DATA_DIR} must exist"
    assert (DATA_DIR / "sample_metadata.csv").exists(), "sample_metadata.csv must exist"
    assert (DATA_DIR / "test_1.fastq.gz").exists(), "test_1.fastq.gz must exist"


def test_01_sample_metadata_validator_node(setup_dirs):
    """Test SampleMetadataValidatorNode on real nf-core dataset."""
    node = SampleMetadataValidatorNode()
    result = node.run(
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
    )
    assert len(result) == 1
    assert Path(result[0]).exists()
    print("\n[PASS] SampleMetadataValidatorNode validated samples successfully.")


def test_02_fastp_qc_node(setup_dirs):
    """Test FastpQCNode executing real fastp inside bulk_rna_seq conda env."""
    qc_out = OUTPUT_BASE / "qc"
    node = FastpQCNode()
    result = node.run(
        fastq_pair="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(qc_out),
        threads=2,
    )
    assert len(result) == 1
    assert qc_out.exists()
    assert any(qc_out.glob("*.json"))
    print(f"\n[PASS] FastpQCNode generated QC artifacts in {qc_out}.")


def test_03_fastp_trim_node(setup_dirs):
    """Test FastpTrimNode executing real fastp trimming."""
    trim_out = OUTPUT_BASE / "trimmed"
    node = FastpTrimNode()
    result = node.run(
        fastp_qc_json="",
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        output_dir=str(trim_out),
        threads=2,
    )
    assert len(result) == 1
    assert trim_out.exists()
    trimmed_files = list(trim_out.rglob("*.fastq")) + list(trim_out.rglob("*.fastq.gz"))
    assert len(trimmed_files) >= 2, f"Expected trimmed fastq files, found: {trimmed_files}"
    print(f"\n[PASS] FastpTrimNode trimmed fastq files: {[f.name for f in trimmed_files]}.")


def test_04_salmon_index_node(setup_dirs):
    """Test SalmonIndexNode generating a real Salmon index from transcriptome FASTA."""
    index_out = OUTPUT_BASE / "salmon_index"
    node = SalmonIndexNode()
    result = node.run(
        transcriptome_fasta_path="",
        transcriptome_fasta=str(DATA_DIR / "transcriptome.fasta"),
        index_dir=str(index_out),
        threads=2,
        extra_command="-k 15",  # small k-mer for test transcriptome
    )
    assert len(result) == 1
    assert index_out.exists()
    assert (index_out / "versionInfo.json").exists() or (index_out / "info.json").exists()
    print(f"\n[PASS] SalmonIndexNode successfully indexed transcriptome in {index_out}.")


def test_05_salmon_quant_node(setup_dirs):
    """Test SalmonQuantNode quantifying reads against Salmon index."""
    index_out = OUTPUT_BASE / "salmon_index"
    trim_out = OUTPUT_BASE / "trimmed"
    quant_out = OUTPUT_BASE / "salmon_quant"
    node = SalmonQuantNode()
    result = node.run(
        salmon_index_dir="",
        index_dir=str(index_out),
        fastq_dir=str(DATA_DIR),
        metadata_csv=str(DATA_DIR / "sample_metadata.csv"),
        trimmed_dir=str(trim_out),
        output_dir=str(quant_out),
        read_layout="A",
        threads=2,
    )
    assert len(result) == 1
    assert quant_out.exists()
    quant_files = list(quant_out.rglob("quant.sf"))
    assert len(quant_files) >= 1, f"Expected quant.sf in {quant_out}, found {quant_files}"
    print(f"\n[PASS] SalmonQuantNode successfully generated quant.sf: {[f.relative_to(quant_out) for f in quant_files]}.")
