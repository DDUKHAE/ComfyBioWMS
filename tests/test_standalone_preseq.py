import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/preseq.py").resolve()

@pytest.fixture(scope="module")
def preseq_module():
    spec = importlib.util.spec_from_file_location("standalone_preseq", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_preseq_reports_missing_binary(preseq_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="preseq executable not found"):
        preseq_module.Preseq().run(str(bam), str(tmp_path / "out"))
