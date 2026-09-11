import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/trim_galore.py").resolve()

@pytest.fixture(scope="module")
def trim_galore_module():
    spec = importlib.util.spec_from_file_location("standalone_trim_galore", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_trim_galore_exposes_galaxy_options_and_extra_command(trim_galore_module):
    inputs = trim_galore_module.TrimGalore.INPUT_TYPES()
    assert set(inputs["required"]) == {"reads_fwd"}
    assert "quality" in inputs["optional"]
    assert "adapter" in inputs["optional"]
    assert "extra_command" in inputs["optional"]
    assert set(trim_galore_module.NODE_CLASS_MAPPINGS) == {"TrimGalore"}

def test_trim_galore_reports_missing_binary(trim_galore_module, tmp_path, monkeypatch):
    dummy_fq = tmp_path / "dummy.fq"
    dummy_fq.write_text("@r1\nACGT\n+\nIIII\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="Trim Galore! executable not found"):
        trim_galore_module.TrimGalore().run(str(dummy_fq), str(tmp_path / "out"))
