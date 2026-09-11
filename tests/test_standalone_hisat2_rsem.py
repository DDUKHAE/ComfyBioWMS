import importlib.util
from pathlib import Path
import pytest

M_HISAT2 = Path("nodes/class_2/hisat2.py").resolve()
M_RSEM = Path("nodes/class_2/rsem.py").resolve()

@pytest.fixture(scope="module")
def hisat2_module():
    spec = importlib.util.spec_from_file_location("standalone_hisat2", M_HISAT2)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

@pytest.fixture(scope="module")
def rsem_module():
    spec = importlib.util.spec_from_file_location("standalone_rsem", M_RSEM)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def test_hisat2_reports_missing_binary(hisat2_module, tmp_path, monkeypatch):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="hisat2-build executable not found"):
        hisat2_module.HISAT2Build().run(str(fa), str(tmp_path / "out"))

def test_rsem_reports_missing_binary(rsem_module, tmp_path, monkeypatch):
    bam = tmp_path / "tx.bam"
    bam.write_bytes(b"BAM")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="rsem-calculate-expression executable not found"):
        rsem_module.RSEMCalculateExpression().run(str(bam), "ref_prefix", str(tmp_path / "out"))
