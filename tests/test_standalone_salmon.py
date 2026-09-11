import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/salmon.py").resolve()

@pytest.fixture(scope="module")
def salmon_module():
    spec = importlib.util.spec_from_file_location("standalone_salmon", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_salmon_exposes_inputs_and_mappings(salmon_module):
    idx_inputs = salmon_module.SalmonIndex.INPUT_TYPES()
    assert set(idx_inputs["required"]) == {"transcripts_fasta"}
    assert "kmer_len" in idx_inputs["optional"]

    quant_inputs = salmon_module.SalmonQuantReads.INPUT_TYPES()
    assert set(quant_inputs["required"]) == {"salmon_index_dir", "reads_fwd"}
    assert "strandedness" in quant_inputs["optional"]
    assert "validate_mappings" in quant_inputs["optional"]

    bam_inputs = salmon_module.SalmonQuantAlignment.INPUT_TYPES()
    assert set(bam_inputs["required"]) == {"transcripts_fasta", "transcriptome_bam"}

    assert set(salmon_module.NODE_CLASS_MAPPINGS) == {
        "SalmonIndex",
        "SalmonQuantReads",
        "SalmonQuantAlignment",
    }

def test_salmon_reports_missing_binary(salmon_module, tmp_path, monkeypatch):
    dummy_fa = tmp_path / "tx.fa"
    dummy_fa.write_text(">tx1\nACGTACGT\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="Salmon executable not found"):
        salmon_module.SalmonIndex().run(str(dummy_fa), str(tmp_path / "out"))
