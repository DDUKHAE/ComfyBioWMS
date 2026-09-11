import importlib.util
import json
from pathlib import Path
import pytest

MODULE = Path("nodes/class_1/biopython.py").resolve()


def load_module():
    spec = importlib.util.spec_from_file_location("standalone_biopython", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_biopython_module_loads_standalone_and_owns_mappings():
    mod = load_module()
    assert len(mod.NODE_CLASS_MAPPINGS) == 8
    assert set(mod.NODE_CLASS_MAPPINGS) == set(mod.NODE_DISPLAY_NAME_MAPPINGS)
    expected_classes = [
        "BiopythonSeqIOStats",
        "BiopythonAlignmentStats",
        "BiopythonSeqTransform",
        "BiopythonGCContent",
        "BiopythonPairwiseAlign",
        "BiopythonProtParam",
        "BiopythonRestrictionDigest",
        "BiopythonSeqFilter",
    ]
    for name in expected_classes:
        assert name in mod.NODE_CLASS_MAPPINGS
        cls = mod.NODE_CLASS_MAPPINGS[name]
        assert hasattr(cls, "INPUT_TYPES")
        assert hasattr(cls, "RETURN_TYPES")
        assert cls.CATEGORY == "ComfyBIO/Biopython"
        assert cls.FUNCTION == "run"


def test_seq_transform_reverse_complement_transcribe_translate():
    mod = load_module()
    node = mod.BiopythonSeqTransform()

    # DNA: ATGCGATCGATCGATCGATAG
    # RevComp: CTATCGATCGATCGATCGCAT
    # Comp: TACGCTAGCTAGCTAGCTATC
    # RNA: AUGCGAUCGAUCGAUCGAUAG
    # Protein: MRSIDR*
    rev_comp, comp, rna, protein = node.run(sequence="ATGCGATCGATCGATCGATAG")
    assert rev_comp == "CTATCGATCGATCGATCGCAT"
    assert comp == "TACGCTAGCTAGCTAGCTATC"
    assert rna == "AUGCGAUCGAUCGAUCGAUAG"
    assert protein.startswith("MRSIDR")


def test_gc_content_and_metrics():
    mod = load_module()
    node = mod.BiopythonGCContent()

    # Sequence: 10 bp, 5 GC (50.0%), 5 AT (50.0%)
    gc_pct, at_pct, length, mw, tm, summary_json = node.run(sequence="ATGCATGCGC")
    assert gc_pct == 60.0
    assert at_pct == 40.0
    assert length == 10
    assert mw > 0
    assert tm > 0
    summary = json.loads(summary_json)
    assert summary["length"] == 10
    assert summary["counts"]["G"] == 3
    assert summary["counts"]["C"] == 3


def test_pairwise_align():
    mod = load_module()
    node = mod.BiopythonPairwiseAlign()

    aln_text, score, pid = node.run(
        seq1="ATGCGATC",
        seq2="ATGCGATC",
        mode="global",
    )
    assert score == 8.0
    assert pid == 100.0
    assert "ATGCGATC" in aln_text


def test_prot_param_physicochemical_properties():
    mod = load_module()
    node = mod.BiopythonProtParam()

    mw, pi, instability, gravy, summary_json = node.run(protein_sequence="MSIVMGRKGAR")
    assert 1100 < mw < 1300
    assert 10.0 < pi <= 14.0
    assert isinstance(instability, float)
    assert isinstance(gravy, float)
    summary = json.loads(summary_json)
    assert summary["length"] == 11
    assert "helix" in summary["secondary_structure_fractions"]


def test_restriction_digest_cut_sites_and_fragments():
    mod = load_module()
    node = mod.BiopythonRestrictionDigest()

    # EcoRI: GAATTC (cut after G at 1, 1-based pos 2)
    # BamHI: GGATCC (cut after G at 1, 1-based pos 8)
    dna = "GAATTCGGATCC"
    cut_sites_json, fragments_json, total_cuts = node.run(
        dna_sequence=dna,
        enzymes="EcoRI,BamHI",
        linear=True,
    )
    cut_sites = json.loads(cut_sites_json)
    fragments = json.loads(fragments_json)
    assert "EcoRI" in cut_sites
    assert "BamHI" in cut_sites
    assert total_cuts == 2
    assert len(fragments) >= 2


def test_seqio_stats_and_filter_on_real_fasta(tmp_path):
    mod = load_module()
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">seq1 High GC\n"
        "GCGCGCGCGCGCGCGCGCGC\n"
        ">seq2 Short\n"
        "ATGC\n"
        ">seq3 Normal\n"
        "ATGCGATCGATCGATCGATCGATCGATC\n"
    )

    # 1. SeqIOStats
    stats_json, count = mod.BiopythonSeqIOStats().run(str(fasta), "fasta")
    assert count == 3
    stats = json.loads(stats_json)
    assert stats[0]["gc_percent"] == 100.0
    assert stats[1]["length"] == 4

    # 2. SeqFilter
    filtered_file, passed, total = mod.BiopythonSeqFilter().run(
        sequence_file=str(fasta),
        min_length=10,
        max_length=100,
        min_gc=20.0,
        max_gc=80.0,
        file_format="fasta",
        output_dir=str(tmp_path / "filtered_out"),
    )
    assert total == 3
    assert passed == 1  # only seq3 (normal) passes
    assert Path(filtered_file).is_file()


def test_alignment_stats_on_real_clustal(tmp_path):
    mod = load_module()
    aln_file = tmp_path / "test.aln"
    aln_file.write_text(
        "CLUSTAL W (1.82) multiple sequence alignment\n\n"
        "seq1      ATGCGATCG\n"
        "seq2      ATGCGATCG\n"
        "seq3      ATGCCATCG\n"
        "          **** ****\n"
    )
    summary_json, rows, cols, identity = mod.BiopythonAlignmentStats().run(str(aln_file), "clustal")
    assert rows == 3
    assert cols == 9
    assert identity > 80.0


def test_biopython_rejects_missing_file():
    mod = load_module()
    with pytest.raises(FileNotFoundError, match="Sequence input is not a file"):
        mod.BiopythonSeqIOStats().run("missing.fasta", "fasta")

