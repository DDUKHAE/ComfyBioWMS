import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/multiqc.py").resolve()

@pytest.fixture(scope="module")
def multiqc_module():
    spec = importlib.util.spec_from_file_location("standalone_multiqc", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_multiqc_exposes_inputs_and_mappings(multiqc_module):
    inputs = multiqc_module.MultiQC.INPUT_TYPES()
    assert set(inputs["required"]) == {"analysis_dir"}
    assert "report_filename" in inputs["optional"]
    assert "config_yaml" in inputs["optional"]
    assert set(multiqc_module.NODE_CLASS_MAPPINGS) == {"MultiQC"}

def test_multiqc_reports_missing_binary(multiqc_module, tmp_path, monkeypatch):
    scan_dir = tmp_path / "scan"
    scan_dir.mkdir()
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="MultiQC executable not found"):
        multiqc_module.MultiQC().run(str(scan_dir), str(tmp_path / "out"))
