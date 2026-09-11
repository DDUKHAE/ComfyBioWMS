import importlib.util
from pathlib import Path
import pytest

M_KRAKEN = Path("nodes/class_2/kraken2.py").resolve()
M_BRACKEN = Path("nodes/class_2/bracken.py").resolve()

@pytest.fixture(scope="module")
def kraken_module():
    spec = importlib.util.spec_from_file_location("standalone_kraken2", M_KRAKEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

@pytest.fixture(scope="module")
def bracken_module():
    spec = importlib.util.spec_from_file_location("standalone_bracken", M_BRACKEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def test_kraken2_reports_missing_binary(kraken_module, tmp_path, monkeypatch):
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    db = tmp_path / "kdb"
    db.mkdir()
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="kraken2 executable not found"):
        kraken_module.Kraken2Classify().run(str(fwd), str(db), str(tmp_path / "out"))

def test_bracken_reports_missing_binary(bracken_module, tmp_path, monkeypatch):
    rep = tmp_path / "rep.txt"
    rep.write_text("10.0\t100\t10\tS\t1234\tTaxon")
    db = tmp_path / "bdb"
    db.mkdir()
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bracken executable not found"):
        bracken_module.Bracken().run(str(rep), str(db), str(tmp_path / "out"))
