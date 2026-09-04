import importlib.util
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/spades.py").resolve()


@pytest.fixture(scope="module")
def spades_module():
    spec = importlib.util.spec_from_file_location("standalone_spades", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_spades_galaxy_options_and_extra_command_are_optional(spades_module):
    inputs = spades_module.SpadesNode.INPUT_TYPES()
    assert set(inputs["required"]) == {"read1", "output_dir"}
    assert inputs["optional"]["threads"][1]["default"] == 16
    assert inputs["optional"]["memory_gb"][1]["default"] == 16
    assert inputs["optional"]["cov_cutoff"][1]["default"] == "off"
    assert inputs["optional"]["kmers"][1]["default"] == "auto"
    assert inputs["optional"]["phred_offset"][1]["default"] == "auto"
    assert inputs["optional"]["careful"][1]["default"] is False
    assert inputs["optional"]["sc"][1]["default"] is False
    assert inputs["optional"]["meta"][1]["default"] is False
    assert inputs["optional"]["isolate"][1]["default"] is False
    assert inputs["optional"]["only_assembler"][1]["default"] is False
    assert "extra_command" in inputs["optional"]
    assert set(spades_module.NODE_CLASS_MAPPINGS) == {"SpadesNode"}
    assert set(spades_module.NODE_DISPLAY_NAME_MAPPINGS) == {"SpadesNode"}


def test_spades_removes_managed_options_but_keeps_native_extra(spades_module):
    kept, ignored = spades_module._filter_extra(
        "-t8 -m16 -k 21,33,55 --threads=8 --memory=16 --careful --cov-cutoff=auto --phred-offset=33 --disable-rr --only-error-correction --tmp-dir /tmp/other",
        spades_module._MANAGED_OPTIONS,
    )
    assert kept == ["--disable-rr", "--only-error-correction", "--tmp-dir", "/tmp/other"]
    assert "--only-error-correction" not in ignored
    assert "-t8" in ignored
    assert "-m16" in ignored
    assert "-k" in ignored
    assert "--threads=8" in ignored
    assert "--memory=16" in ignored
    assert "--careful" in ignored
    assert "--cov-cutoff=auto" in ignored
    assert "--phred-offset=33" in ignored


def test_spades_reports_missing_binary_before_creating_output(spades_module, tmp_path, monkeypatch):
    r1 = fetch_official_data("spades_1K_1.fq.gz", tmp_path)
    output_dir = tmp_path / "absent"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="spades.py executable not found on PATH"):
        spades_module.SpadesNode().run(str(r1), str(output_dir))
    assert not output_dir.exists()


def test_spades_reports_missing_input_before_creating_output(spades_module, tmp_path):
    r1 = tmp_path / "absent_read1.fastq.gz"
    output_dir = tmp_path / "absent"
    with pytest.raises(FileNotFoundError, match="is not a file"):
        spades_module.SpadesNode().run(str(r1), str(output_dir))
    assert not output_dir.exists()


@pytest.mark.e2e
def test_spades_runs_on_official_paired_dataset_and_produces_contigs(spades_module, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r1 = fetch_official_data("spades_1K_1.fq.gz", tmp_path)
    r2 = fetch_official_data("spades_1K_2.fq.gz", tmp_path)
    contigs, scaffolds, log = spades_module.SpadesNode().run(
        str(r1),
        "spades_out",
        read2=str(r2),
        threads=2,
        memory_gb=4,
        extra_command="--frugal -t 99",
    )
    assert Path(contigs).is_absolute()
    assert Path(scaffolds).is_absolute()
    assert Path(log).is_absolute()
    assert Path(contigs).stat().st_size > 0

    records = []
    current_header = None
    current_seq = []
    with Path(contigs).open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                current_header = stripped[1:]
                current_seq = []
            else:
                current_seq.append(stripped)
        if current_header is not None:
            records.append((current_header, "".join(current_seq)))

    assert len(records) == 1
    header, seq = records[0]
    assert "NODE_1_length_1000" in header
    assert len(seq) == 1000
    assert Path(scaffolds).stat().st_size > 0
    assert Path(log).stat().st_size > 0
