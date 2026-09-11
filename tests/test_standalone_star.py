import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/star.py").resolve()

@pytest.fixture(scope="module")
def star_module():
    spec = importlib.util.spec_from_file_location("standalone_star", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_star_exposes_galaxy_options_and_extra_command(star_module):
    gen_inputs = star_module.STARGenomeGenerate.INPUT_TYPES()
    assert set(gen_inputs["required"]) == {"genome_fasta"}
    assert "gtf_file" in gen_inputs["optional"]
    assert "threads" in gen_inputs["optional"]

    align_inputs = star_module.STARAlignReads.INPUT_TYPES()
    assert set(align_inputs["required"]) == {"star_index_dir", "reads_fwd"}
    assert "quant_mode" in align_inputs["optional"]
    assert "twopass_mode" in align_inputs["optional"]
    assert set(star_module.NODE_CLASS_MAPPINGS) == {"STARGenomeGenerate", "STARAlignReads"}

def test_star_reports_missing_binary(star_module, tmp_path, monkeypatch):
    dummy_fa = tmp_path / "ref.fa"
    dummy_fa.write_text(">chr1\nACGTACGT\n")
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="STAR executable not found"):
        star_module.STARGenomeGenerate().run(str(dummy_fa), str(tmp_path / "out"))
