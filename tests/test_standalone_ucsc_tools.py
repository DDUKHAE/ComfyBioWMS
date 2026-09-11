import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/ucsc_tools.py").resolve()

@pytest.fixture(scope="module")
def ucsc_module():
    spec = importlib.util.spec_from_file_location("standalone_ucsc_tools", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_ucsc_reports_missing_binary(ucsc_module, tmp_path, monkeypatch):
    bg = tmp_path / "test.bedgraph"
    bg.write_text("chr1\t100\t200\t10\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000000\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bedGraphToBigWig executable not found"):
        ucsc_module.BedGraphToBigWig().run(str(bg), str(sizes), str(tmp_path / "out"))
