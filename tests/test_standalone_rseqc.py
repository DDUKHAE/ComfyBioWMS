import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/rseqc.py").resolve()

@pytest.fixture(scope="module")
def rseqc_module():
    spec = importlib.util.spec_from_file_location("standalone_rseqc", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_rseqc_reports_missing_binary(rseqc_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM")
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t10\t20\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="geneBody_coverage.py not found"):
        rseqc_module.RSeQCGeneBodyCoverage().run(str(bam), str(bed), str(tmp_path / "out"))
