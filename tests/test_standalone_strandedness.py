import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/strandedness.py").resolve()

@pytest.fixture(scope="module")
def stranded_module():
    spec = importlib.util.spec_from_file_location("standalone_strandedness", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_strandedness_reports_missing_binary(stranded_module, tmp_path, monkeypatch):
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir()
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="Salmon executable not found"):
        stranded_module.InferStrandedness().run(str(idx_dir), str(fwd), str(tmp_path / "out"))
