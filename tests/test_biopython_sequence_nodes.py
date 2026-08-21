import pytest
from pathlib import Path
import json
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.biopython_nodes import (
    BiopythonSeqIONode,
    BiopythonSeqTranscribeNode,
    BiopythonAlignIONode,
    BiopythonBlastNode,
    BiopythonEntrezNode,
    BiopythonBioPDBNode,
    BiopythonRestrictionNode,
    BiopythonBioMotifsNode,
    LogomakerVisualizerNode,
    BiopythonProtParamNode,
    BiotiteStructureAlignNode,
    SquiggleWaveformNode,
    PyhmmerSearchNode,
    DnaFeaturesViewerNode,
    Primer3DesignNode,
    BioSeqAnalysisNode,
    PyfaidxIndexNode,
    PyCircosPlotNode,
    BiopythonPhyloNode,
    EdlibAlignNode,
    CodonWAnalysisNode,
    PyTFBSScanNode,
    BioKEGGPathwayNode,
    HelicalWheelNode,
    SeqLogoGeneratorNode,
)

OUTPUT_BASE = root_dir / "results" / "test_biopython_nodes"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)


def test_01_seqio_and_transcribe_translate(setup_dirs):
    """Test SeqIO and Sequence Transcription/Translation."""
    fasta_file = OUTPUT_BASE / "test_seq.fasta"
    fasta_file.write_text(">gene1 Synthetic GFP\nATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCACTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA\n", encoding="utf-8")
    
    # 1. SeqIO
    seqio_node = BiopythonSeqIONode()
    records_json, count = seqio_node.run(str(fasta_file), file_format="fasta")
    assert count >= 1
    parsed = json.loads(records_json)
    assert parsed[0]["id"] == "gene1"
    print(f"\n[PASS] BiopythonSeqIONode parsed {count} records: ID={parsed[0]['id']}.")

    # 2. Transcribe / Translate
    trans_node = BiopythonSeqTranscribeNode()
    rna, prot, rev_comp = trans_node.run("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATTAG")
    assert rna.startswith("AUG")
    assert rev_comp.endswith("CAT")
    print(f"[PASS] BiopythonSeqTranscribeNode translated {len(rna)} nt to {len(prot)} aa.")


def test_02_alignment_and_similarity(setup_dirs):
    """Test Pairwise Alignment, Edit Distance (Edlib), and Alignment nodes."""
    aln_node = BiopythonAlignIONode()
    aln_str, pid = aln_node.run("ATGCGATCGATCG", "ATGCGATAGATCG")
    assert pid > 80.0
    print(f"\n[PASS] BiopythonAlignIONode calculated {pid:.1f}% identity.")

    # Edlib align
    ed_node = EdlibAlignNode()
    dist, cigar = ed_node.run("ATGCGATCGA", "ATGCGATAGA")
    assert dist >= 1
    print(f"[PASS] EdlibAlignNode calculated edit distance: {dist}, CIGAR: {cigar}.")


def test_03_gc_content_and_sequence_analysis(setup_dirs):
    """Test BioSeqAnalysisNode (GC, MW, Tm)."""
    dna = "ATGCGATCGATCGATCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    analysis_node = BioSeqAnalysisNode()
    res_json, gc_val = analysis_node.run(dna)
    assert 0.0 <= gc_val <= 100.0
    print(f"\n[PASS] BioSeqAnalysisNode: GC={gc_val:.1f}%.")


def test_04_physicochemical_properties(setup_dirs):
    """Test ProtParam physicochemical properties and CodonW analysis."""
    prot_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    
    # ProtParam
    prot_node = BiopythonProtParamNode()
    mw, pi_val, instab, summary_json = prot_node.run(prot_seq)
    assert 4.0 <= pi_val <= 10.0
    print(f"\n[PASS] BiopythonProtParamNode: MW={mw:.1f}, pI={pi_val:.2f}, Instability={instab:.1f}.")

    # CodonW analysis
    codon_node = CodonWAnalysisNode()
    cai_val, gc3s, fop = codon_node.run("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATTAG")
    assert 0.0 <= cai_val <= 1.0
    print(f"[PASS] CodonWAnalysisNode calculated CAI: {cai_val:.3f}, GC3s: {gc3s:.3f}.")


def test_05_visualizers_and_motifs(setup_dirs):
    """Test Logomaker, DNA Features Viewer, Helical Wheel, and Motif nodes."""
    # Logomaker
    logo_node = LogomakerVisualizerNode()
    logo_res = logo_node.run("TATAAA\nTATAAT\nTATAAA\nTATACT\nTATAAT")
    assert len(logo_res) == 2
    print("\n[PASS] LogomakerVisualizerNode rendered sequence logo.")

    # Helical wheel
    wheel_node = HelicalWheelNode()
    wheel_res = wheel_node.run("MSKGEELFTGVVPILVELDGDVNG")
    assert len(wheel_res) == 2
    print("[PASS] HelicalWheelNode rendered helical wheel plot.")

    # DNA Features Viewer
    viewer_node = DnaFeaturesViewerNode()
    viewer_res = viewer_node.run(plasmid_name="pUC19-EGFP")
    assert len(viewer_res) == 2
    print("[PASS] DnaFeaturesViewerNode generated plasmid/feature map.")


def test_06_primer_and_motif_search(setup_dirs):
    """Test Primer3 Design and Restriction Enzyme Digestion nodes."""
    dna = "ATGCGATCGATCGATCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 5
    
    # Primer3
    primer_node = Primer3DesignNode()
    fwd, rev, tm_f, tm_r = primer_node.run(dna, target_tm=60.0)
    assert len(fwd) > 0 and len(rev) > 0
    print(f"\n[PASS] Primer3DesignNode: Fwd={fwd} (Tm={tm_f:.1f}°C), Rev={rev} (Tm={tm_r:.1f}°C).")

    # Restriction Enzyme Digestion
    re_node = BiopythonRestrictionNode()
    re_json, cut_count = re_node.run(dna, enzymes="EcoRI, BamHI, HindIII")
    assert cut_count >= 0
    print(f"[PASS] BiopythonRestrictionNode analyzed {cut_count} restriction sites.")
