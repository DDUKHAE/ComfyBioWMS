import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/sortmerna.py").resolve()

@pytest.fixture(scope="module")
def sortmerna_module():
    spec = importlib.util.spec_from_file_location("standalone_sortmerna", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_sortmerna_reports_missing_binary(sortmerna_module, tmp_path, monkeypatch):
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    db = tmp_path / "rrna.fa"
    db.write_text(">rrna\nACGT\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="sortmerna executable not found"):
        sortmerna_module.SortMeRNA().run(str(fwd), str(db), str(tmp_path / "out"))
