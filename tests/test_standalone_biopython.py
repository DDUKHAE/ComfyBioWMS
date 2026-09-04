import importlib.util
import json
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_1/biopython.py").resolve()


def load_module():
    spec = importlib.util.spec_from_file_location("standalone_biopython", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_biopython_file_loads_standalone_and_owns_mappings():
    module = load_module()
    assert set(module.NODE_CLASS_MAPPINGS) == {
        "BiopythonSeqIOStatsNode",
        "BiopythonAlignmentStatsNode",
    }
    assert set(module.NODE_DISPLAY_NAME_MAPPINGS) == set(module.NODE_CLASS_MAPPINGS)


@pytest.mark.e2e
def test_seqio_stats_parse_official_fasta(tmp_path):
    module = load_module()
    fasta = fetch_official_data("biopython_dups.fasta", tmp_path)
    summary, count = module.BiopythonSeqIOStatsNode().run(str(fasta), "fasta")
    rows = json.loads(summary)
    assert count == 5
    assert len({row["id"] for row in rows}) == 4


@pytest.mark.e2e
def test_alignment_stats_parse_official_clustal(tmp_path):
    module = load_module()
    alignment = fetch_official_data("biopython_protein.aln", tmp_path)
    summary, rows, columns, identity = module.BiopythonAlignmentStatsNode().run(
        str(alignment), "clustal"
    )
    assert rows == 20
    assert columns == 411
    assert 0 <= identity <= 100
    assert json.loads(summary)["rows"] == 20


def test_biopython_nodes_reject_missing_files():
    module = load_module()
    with pytest.raises(FileNotFoundError, match="Sequence input is not a file"):
        module.BiopythonSeqIOStatsNode().run("missing.fasta", "fasta")
