import importlib.util
import json
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/fastp.py").resolve()


@pytest.fixture(scope="module")
def fastp_module():
    spec = importlib.util.spec_from_file_location("standalone_fastp", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fastp_galaxy_options_and_extra_command_are_optional(fastp_module):
    inputs = fastp_module.FastpNode.INPUT_TYPES()
    assert set(inputs["required"]) == {"read1", "output_dir"}
    assert inputs["optional"]["qualified_quality_phred"][1]["default"] == 15
    assert inputs["optional"]["unqualified_percent_limit"][1]["default"] == 40
    assert inputs["optional"]["length_required"][1]["default"] == 15
    assert "extra_command" in inputs["optional"]
    assert set(fastp_module.NODE_CLASS_MAPPINGS) == {"FastpNode"}


def test_fastp_removes_managed_options_but_keeps_native_extra(fastp_module):
    kept, ignored = fastp_module._filter_extra(
        "--thread=99 -q 3 --overrepresentation_analysis",
        fastp_module._MANAGED_OPTIONS,
    )
    assert kept == ["--overrepresentation_analysis"]
    assert ignored == ["--thread=99", "-q"]


def test_fastp_reports_missing_binary_before_creating_output(fastp_module, tmp_path, monkeypatch):
    read1 = fetch_official_data("fastp_R1.fq", tmp_path)
    output_dir = tmp_path / "absent"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="fastp executable not found on PATH"):
        fastp_module.FastpNode().run(str(read1), str(output_dir))
    assert not output_dir.exists()


@pytest.mark.e2e
def test_fastp_runs_on_upstream_paired_data(fastp_module, tmp_path):
    r1 = fetch_official_data("fastp_R1.fq", tmp_path)
    r2 = fetch_official_data("fastp_R2.fq", tmp_path)
    out1, out2, report_json, report_html = fastp_module.FastpNode().run(
        str(r1),
        str(tmp_path / "out"),
        str(r2),
        threads=2,
        extra_command="--overrepresentation_analysis --thread 99",
    )
    assert Path(out1).stat().st_size > 0
    assert Path(out2).stat().st_size > 0
    report = json.loads(Path(report_json).read_text())
    assert report["summary"]["before_filtering"]["total_reads"] > 0
    assert Path(report_html).stat().st_size > 0
