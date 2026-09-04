import importlib.util
import zipfile
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/fastqc.py").resolve()


@pytest.fixture(scope="module")
def fastqc_module():
    spec = importlib.util.spec_from_file_location("standalone_fastqc", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fastqc_exposes_galaxy_options_and_optional_extra_command(fastqc_module):
    inputs = fastqc_module.FastQCNode.INPUT_TYPES()
    assert set(inputs["required"]) == {"input_file", "output_dir"}
    assert inputs["optional"]["kmers"][1] == {"default": 7, "min": 2, "max": 10}
    assert "extra_command" in inputs["optional"]
    assert set(fastqc_module.NODE_CLASS_MAPPINGS) == {"FastQCNode"}


def test_fastqc_filters_managed_options_and_keeps_native_extra(fastqc_module):
    kept, ignored = fastqc_module._filter_extra(
        "--kmers=3 -t 99 --noextract", fastqc_module._MANAGED_OPTIONS
    )
    assert kept == ["--noextract"]
    assert ignored == ["--kmers=3", "-t"]


def test_fastqc_reports_missing_binary_before_creating_output(fastqc_module, tmp_path, monkeypatch):
    input_file = fetch_official_data("fastqc_minimal.fastq", tmp_path)
    output_dir = tmp_path / "absent"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="FastQC executable not found on PATH"):
        fastqc_module.FastQCNode().run(str(input_file), str(output_dir))
    assert not output_dir.exists()


@pytest.mark.e2e
def test_fastqc_runs_on_official_minimal_fastq(fastqc_module, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = fetch_official_data("fastqc_minimal.fastq", tmp_path)
    html, archive = fastqc_module.FastQCNode().run(
        str(input_file), "out", threads=2, kmers=7
    )
    assert Path(html).is_absolute()
    assert Path(archive).is_absolute()
    assert Path(html).stat().st_size > 0
    assert Path(archive).stat().st_size > 0
    with zipfile.ZipFile(archive) as zipped:
        data_name = next(name for name in zipped.namelist() if name.endswith("fastqc_data.txt"))
        assert zipped.read(data_name).startswith(b"##FastQC")
