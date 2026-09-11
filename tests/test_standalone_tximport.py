import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_1/tximport.py").resolve()

@pytest.fixture(scope="module")
def tximport_module():
    spec = importlib.util.spec_from_file_location("standalone_tximport", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_tximport_aggregates_salmon_quant(tximport_module, tmp_path):
    s1_dir = tmp_path / "sample1"
    s1_dir.mkdir()
    s1_quant = s1_dir / "quant.sf"
    s1_quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTX1\t1000\t800\t10.5\t100.0\nTX2\t2000\t1800\t5.0\t50.0\n")

    s2_dir = tmp_path / "sample2"
    s2_dir.mkdir()
    s2_quant = s2_dir / "quant.sf"
    s2_quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTX1\t1000\t800\t20.0\t200.0\nTX2\t2000\t1800\t15.0\t150.0\n")

    tx2gene = tmp_path / "tx2gene.tsv"
    tx2gene.write_text("TX1\tGENE_A\nTX2\tGENE_A\n")

    counts_tsv, tpm_tsv, summary_json = tximport_module.Tximport().run(
        quant_files=f"{s1_quant},{s2_quant}",
        output_dir=str(tmp_path / "out"),
        tx2gene_tsv=str(tx2gene),
        sample_names="S1,S2",
    )

    assert Path(counts_tsv).is_file()
    assert Path(tpm_tsv).is_file()
    
    counts_content = Path(counts_tsv).read_text()
    assert "gene_id\tS1\tS2" in counts_content
    # GENE_A sum in S1: 100 + 50 = 150.0, in S2: 200 + 150 = 350.0
    assert "GENE_A\t150.0000\t350.0000" in counts_content
