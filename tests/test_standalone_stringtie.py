import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/stringtie.py").resolve()

@pytest.fixture(scope="module")
def stringtie_module():
    spec = importlib.util.spec_from_file_location("standalone_stringtie", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_stringtie_reports_missing_binary(stringtie_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM_DUMMY")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="stringtie executable not found"):
        stringtie_module.StringTie().run(str(bam), str(tmp_path / "out"))
