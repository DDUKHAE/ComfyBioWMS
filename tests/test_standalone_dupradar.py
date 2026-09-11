import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/dupradar.py").resolve()

@pytest.fixture(scope="module")
def dupradar_module():
    spec = importlib.util.spec_from_file_location("standalone_dupradar", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_dupradar_reports_missing_binary(dupradar_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM")
    gtf = tmp_path / "test.gtf"
    gtf.write_text("chr1\tref\texon\t1\t10\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="Rscript executable not found"):
        dupradar_module.DupRadar().run(str(bam), str(gtf), str(tmp_path / "out"))
