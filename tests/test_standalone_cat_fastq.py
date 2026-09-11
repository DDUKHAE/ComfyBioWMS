import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_2/cat_fastq.py").resolve()

@pytest.fixture(scope="module")
def cat_module():
    spec = importlib.util.spec_from_file_location("standalone_cat_fastq", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_cat_fastq_merges_files(cat_module, tmp_path):
    f1 = tmp_path / "lane1.fq"
    f1.write_text("@r1\nAAAA\n+\nIIII\n")
    f2 = tmp_path / "lane2.fq"
    f2.write_text("@r2\nTTTT\n+\nJJJJ\n")

    res = cat_module.CatFastq().run(
        fastq_files=f"{f1},{f2}",
        output_dir=str(tmp_path / "out"),
        output_filename="merged.fq"
    )
    merged_path = Path(res[0])
    assert merged_path.is_file()
    assert merged_path.read_text() == "@r1\nAAAA\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n"
