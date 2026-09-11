import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/picard.py").resolve()

@pytest.fixture(scope="module")
def picard_module():
    spec = importlib.util.spec_from_file_location("standalone_picard", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_picard_reports_missing_binary(picard_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM_DUMMY")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="picard executable not found"):
        picard_module.PicardMarkDuplicates().run(str(bam), str(tmp_path / "out"))
