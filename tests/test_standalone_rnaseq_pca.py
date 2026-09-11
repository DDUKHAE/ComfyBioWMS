import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_1/rnaseq_pca.py").resolve()

@pytest.fixture(scope="module")
def pca_module():
    spec = importlib.util.spec_from_file_location("standalone_rnaseq_pca", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_deseq2_sample_qc_runs(pca_module, tmp_path):
    tsv = tmp_path / "counts.tsv"
    tsv.write_text(
        "gene_id\tS1\tS2\tS3\n"
        "G1\t100\t200\t300\n"
        "G2\t500\t400\t600\n"
        "G3\t10\t20\t30\n"
    )
    pca_png, heat_png, summary = pca_module.DESeq2SampleQC().run(str(tsv), str(tmp_path / "out"))
    assert Path(pca_png).is_file()
    assert Path(heat_png).is_file()
    assert "S1" in summary
