import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/bedtools.py").resolve()

@pytest.fixture(scope="module")
def bedtools_module():
    spec = importlib.util.spec_from_file_location("standalone_bedtools", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_bedtools_reports_missing_binary(bedtools_module, tmp_path, monkeypatch):
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"BAM_DUMMY")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bedtools executable not found"):
        bedtools_module.BedtoolsGenomeCoverage().run(str(bam), str(tmp_path / "out"))
