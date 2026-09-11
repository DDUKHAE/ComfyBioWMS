import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/umi_tools.py").resolve()

@pytest.fixture(scope="module")
def umi_module():
    spec = importlib.util.spec_from_file_location("standalone_umi_tools", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_umi_tools_reports_missing_binary(umi_module, tmp_path, monkeypatch):
    fwd = tmp_path / "r1.fq"
    fwd.write_text("@r1\nACGT\n+\nIIII\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="umi_tools executable not found"):
        umi_module.UmiToolsExtract().run(str(fwd), "NNNN", str(tmp_path / "out"))
