import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/bbtools.py").resolve()

@pytest.fixture(scope="module")
def bbtools_module():
    spec = importlib.util.spec_from_file_location("standalone_bbtools", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_bbtools_reports_missing_binary(bbtools_module, tmp_path, monkeypatch):
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    ref = tmp_path / "ref.fa"
    ref.write_text(">host\nACGT\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bbsplit.sh executable not found"):
        bbtools_module.BBSplit().run(str(fwd), str(ref), str(tmp_path / "out"))
