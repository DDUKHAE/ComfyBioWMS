import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/sylph.py").resolve()

@pytest.fixture(scope="module")
def sylph_module():
    spec = importlib.util.spec_from_file_location("standalone_sylph", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_sylph_reports_missing_binary(sylph_module, tmp_path, monkeypatch):
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    db = tmp_path / "test.syldb"
    db.write_bytes(b"SYLPH_DB")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="sylph executable not found"):
        sylph_module.SylphProfile().run(str(fwd), str(db), str(tmp_path / "out"))
