import importlib.util
import inspect
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/quast.py").resolve()
SPADES_MODULE = Path("nodes/class_2/spades.py").resolve()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def quast_module():
    return load_module("standalone_quast", MODULE)


def test_quast_run_signature_has_no_var_keyword(quast_module):
    sig = inspect.signature(quast_module.Quast.run)
    assert not any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
    assert "gff" not in sig.parameters
    assert "features" in sig.parameters


def test_quast_galaxy_options_and_extra_command_are_optional(quast_module):
    inputs = quast_module.Quast.INPUT_TYPES()
    assert set(inputs["required"]) == {"assembly_fasta"}
    assert inputs["optional"]["min_contig"][1]["default"] == 500
    assert inputs["optional"]["threads"][1]["default"] == 4
    assert inputs["optional"]["large"][1]["default"] is False
    assert "reference" in inputs["optional"]
    assert "features" in inputs["optional"]
    assert "extra_command" in inputs["optional"]
    assert set(quast_module.NODE_CLASS_MAPPINGS) == {"Quast"}
    assert set(quast_module.NODE_DISPLAY_NAME_MAPPINGS) == {"Quast"}


def test_quast_reports_missing_binary_before_creating_output(quast_module, tmp_path, monkeypatch):
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nACGTACGT\n")
    output_dir = tmp_path / "absent"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="quast.py executable not found on PATH"):
        quast_module.Quast().run(str(assembly), str(output_dir))
    assert not output_dir.exists()


def test_quast_reports_missing_input_before_creating_output(quast_module, tmp_path):
    assembly = tmp_path / "absent_assembly.fasta"
    output_dir = tmp_path / "absent"
    with pytest.raises(FileNotFoundError, match="is not a file"):
        quast_module.Quast().run(str(assembly), str(output_dir))
    assert not output_dir.exists()


def test_quast_reports_missing_reference_before_creating_output(quast_module, tmp_path):
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nACGTACGT\n")
    reference = tmp_path / "absent_reference.fasta"
    output_dir = tmp_path / "absent"
    with pytest.raises(FileNotFoundError, match="is not a file"):
        quast_module.Quast().run(str(assembly), str(output_dir), reference=str(reference))
    assert not output_dir.exists()


@pytest.mark.e2e
def test_quast_runs_on_spades_assembled_contigs(quast_module, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    spades = load_module("standalone_spades_for_quast", SPADES_MODULE)
    r1 = fetch_official_data("spades_1K_1.fq.gz", tmp_path)
    r2 = fetch_official_data("spades_1K_2.fq.gz", tmp_path)
    contigs, scaffolds, log = spades.Spades().run(
        str(r1),
        "spades_out",
        read2=str(r2),
        threads=2,
        memory_gb=4,
    )
    assert Path(contigs).is_absolute()
    assert Path(contigs).stat().st_size > 0

    report_tsv, report_html = quast_module.Quast().run(
        contigs,
        "quast_out",
        min_contig=100,
        threads=2,
        extra_command="-m 500 --silent",
    )
    assert Path(report_tsv).is_absolute()
    assert Path(report_html).is_absolute()
    assert Path(report_tsv).stat().st_size > 0
    assert Path(report_html).stat().st_size > 0

    metrics = dict(
        line.split("\t", 1)
        for line in Path(report_tsv).read_text().splitlines()
        if "\t" in line
    )
    assert int(metrics.get("# contigs", 0)) > 0
    assert int(metrics.get("Total length", 0)) > 0
    assert int(metrics.get("Largest contig", 0)) > 0
    assert float(metrics.get("GC (%)", 0.0)) > 0.0
