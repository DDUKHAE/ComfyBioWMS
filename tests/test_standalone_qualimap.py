import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/qualimap.py").resolve()

@pytest.fixture(scope="module")
def qualimap_module():
    spec = importlib.util.spec_from_file_location("standalone_qualimap", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_qualimap_reports_missing_binary(qualimap_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM")
    gtf = tmp_path / "test.gtf"
    gtf.write_text("chr1\tref\texon\t1\t10\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="qualimap executable not found"):
        qualimap_module.QualimapRNASeq().run(str(bam), str(gtf), str(tmp_path / "out"))
