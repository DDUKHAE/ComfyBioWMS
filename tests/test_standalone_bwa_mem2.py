import importlib.util
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/bwa_mem2.py").resolve()


@pytest.fixture(scope="module")
def bwa_module():
    spec = importlib.util.spec_from_file_location("standalone_bwa_mem2", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bwa_mem2_file_owns_both_nodes_and_galaxy_defaults(bwa_module):
    assert set(bwa_module.NODE_CLASS_MAPPINGS) == {
        "BwaMem2Index",
        "BwaMem2Align",
    }
    inputs = bwa_module.BwaMem2Align.INPUT_TYPES()
    assert inputs["optional"]["minimum_seed_length"][1]["default"] == 19
    assert inputs["optional"]["band_width"][1]["default"] == 100
    assert inputs["optional"]["minimum_score"][1]["default"] == 30
    assert "extra_command" in inputs["optional"]


def test_bwa_mem2_missing_binary_does_not_create_output(bwa_module, tmp_path, monkeypatch):
    reference = fetch_official_data("sarek_ref.fasta", tmp_path)
    output_dir = tmp_path / "absent"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bwa-mem2 executable not found on PATH"):
        bwa_module.BwaMem2Index().run(str(reference), str(output_dir))
    assert not output_dir.exists()


@pytest.mark.e2e
def test_bwa_mem2_indexes_and_aligns_official_sarek_data(bwa_module, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    reference = fetch_official_data("sarek_ref.fasta", tmp_path)
    read1 = fetch_official_data("sarek_R1.fastq.gz", tmp_path)
    read2 = fetch_official_data("sarek_R2.fastq.gz", tmp_path)
    (indexed_reference,) = bwa_module.BwaMem2Index().run(
        str(reference), "index"
    )
    assert Path(indexed_reference).is_absolute()
    (sam_path,) = bwa_module.BwaMem2Align().run(
        indexed_reference,
        str(read1),
        "aligned/reads.sam",
        read2=str(read2),
        threads=2,
        extra_command="-t 99 -Y",
    )
    assert Path(sam_path).is_absolute()
    lines = Path(sam_path).read_text().splitlines()
    assert lines[0].startswith("@")
    assert any(not line.startswith("@") for line in lines)
