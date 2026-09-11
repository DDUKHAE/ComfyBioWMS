import io
import json
import sys
from pathlib import Path

import pytest

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from bioflow.runtime.command_runner import CondaCommandRunner
from nodes.class_1 import (
    BioKEGGPathwayNode,
    BioSeqAnalysisNode,
    BiopythonAlignIONode,
    BiopythonBioMotifsNode,
    BiopythonBioPDBNode,
    BiopythonEntrezNode,
    BiopythonPhyloNode,
    BiopythonProtParamNode,
    BiopythonRestrictionNode,
    BiopythonSeqIONode,
    BiopythonSeqTranscribeNode,
    BiotiteStructureAlignNode,
    CodonWAnalysisNode,
    DnaFeaturesViewerNode,
    EdlibAlignNode,
    HelicalWheelNode,
    LogomakerVisualizerNode,
    Primer3DesignNode,
    PyCircosPlotNode,
    PyTFBSScanNode,
    PyfaidxIndexNode,
    SeqLogoGeneratorNode,
    SquiggleWaveformNode,
)
from nodes.class_2 import (
    BiopythonBlastNode,
    PyhmmerSearchNode,
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


def test_07_blast_local_and_manifest(setup_dirs):
    """Test BiopythonBlastNode local BLAST+ CLI dispatch and manifest generation."""
    query_fa = OUTPUT_BASE / "blast_query.fasta"
    query_fa.write_text(">query1\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIA\n", encoding="utf-8")

    subject_fa = OUTPUT_BASE / "blast_subject.fasta"
    subject_fa.write_text(">subj1 Serum albumin precursor\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAEEGPKLVAASQAALA\n", encoding="utf-8")

    runner = CondaCommandRunner()
    blast_node = BiopythonBlastNode()
    blast_out_dir = OUTPUT_BASE / "blast_run"
    res_json, top_hit = blast_node.run(
        query_sequence=str(query_fa),
        program="blastp",
        database="swissprot",
        subject_file=str(subject_fa),
        output_dir=str(blast_out_dir),
        runner=runner,
    )

    data = json.loads(res_json)
    assert data["total_hits"] >= 1
    assert data["hits"][0]["identity_percent"] == 100.0
    assert "Serum albumin" in top_hit
    # Verify execution manifest
    assert (blast_out_dir / "run_manifest.sh").exists()
    assert (blast_out_dir / "run_manifest.json").exists()
    print(f"\n[PASS] BiopythonBlastNode dispatched local blastp with {data['total_hits']} hits and generated manifests.")


def test_08_pdb_structure_parsing_and_rmsd_alignment(setup_dirs):
    """Test BiopythonBioPDBNode and BiotiteStructureAlignNode."""
    pdb1_content = """ATOM      1  N   MET A   1      20.154  13.447  11.585  1.00 15.00           N
ATOM      2  CA  MET A   1      19.643  12.123  11.234  1.00 15.00           C
ATOM      3  C   MET A   1      18.234  12.234  10.745  1.00 15.00           C
ATOM      4  O   MET A   1      17.543  13.234  10.987  1.00 15.00           O
ATOM      5  N   ALA A   2      17.843  11.123  10.123  1.00 15.00           N
ATOM      6  CA  ALA A   2      16.456  11.023   9.654  1.00 15.00           C
ATOM      7  C   ALA A   2      15.543  11.890  10.456  1.00 15.00           C
ATOM      8  O   ALA A   2      15.890  13.023  10.789  1.00 15.00           O
TER       9      ALA A   2
END
"""
    pdb2_content = """ATOM      1  N   MET A   1      20.200  13.500  11.600  1.00 15.00           N
ATOM      2  CA  MET A   1      19.700  12.200  11.200  1.00 15.00           C
ATOM      3  C   MET A   1      18.300  12.300  10.800  1.00 15.00           C
ATOM      4  O   MET A   1      17.600  13.300  11.000  1.00 15.00           O
ATOM      5  N   ALA A   2      17.900  11.150  10.150  1.00 15.00           N
ATOM      6  CA  ALA A   2      16.500  11.050   9.700  1.00 15.00           C
ATOM      7  C   ALA A   2      15.600  11.900  10.500  1.00 15.00           C
ATOM      8  O   ALA A   2      15.900  13.050  10.800  1.00 15.00           O
TER       9      ALA A   2
END
"""
    p1 = OUTPUT_BASE / "target.pdb"
    p2 = OUTPUT_BASE / "mobile.pdb"
    p1.write_text(pdb1_content, encoding="utf-8")
    p2.write_text(pdb2_content, encoding="utf-8")

    # 1. PDB Parser
    pdb_node = BiopythonBioPDBNode()
    summary_json, res_count, chain_count = pdb_node.run(str(p1))
    assert res_count == 2
    assert chain_count == 1
    info = json.loads(summary_json)
    assert info["total_atoms"] == 8
    print(f"\n[PASS] BiopythonBioPDBNode parsed {res_count} residues and {info['total_atoms']} atoms.")

    # 2. Structural Alignment
    align_node = BiotiteStructureAlignNode()
    out_pdb = OUTPUT_BASE / "aligned.pdb"
    aligned_path, rmsd = align_node.run(str(p1), str(p2), output_aligned_pdb=str(out_pdb))
    assert Path(aligned_path).exists()
    assert 0.0 <= rmsd < 1.0
    print(f"[PASS] BiotiteStructureAlignNode computed RMSD={rmsd:.3f} Å.")


def test_09_pyfaidx_indexed_fasta(setup_dirs):
    """Test PyfaidxIndexNode for fast chromosome subsequence retrieval."""
    fa = OUTPUT_BASE / "indexed_ref.fasta"
    fa.write_text(">chr1\nACGTACGTACGTACGT\n>chr2\nGGGGCCCCAAAATTTT\n", encoding="utf-8")

    pyfaidx_node = PyfaidxIndexNode()
    subseq, sub_len = pyfaidx_node.run(str(fa), chrom="chr1", start_bp=2, end_bp=8)
    assert subseq == "GTACGT"
    assert sub_len == 6
    print(f"\n[PASS] PyfaidxIndexNode extracted subseq: {subseq} (length={sub_len}).")


def test_10_phylo_tree_and_kegg_pathway(setup_dirs):
    """Test BiopythonPhyloNode and BioKEGGPathwayNode."""
    # Phylo
    nwk_file = OUTPUT_BASE / "test_tree.nwk"
    nwk_file.write_text("(((Human:0.1,Chimp:0.1):0.05,Gorilla:0.15):0.1,Orangutan:0.25);\n", encoding="utf-8")
    phylo_node = BiopythonPhyloNode()
    nwk_out, leaf_count = phylo_node.run(str(nwk_file))
    assert leaf_count == 4
    assert "Human" in nwk_out
    print(f"\n[PASS] BiopythonPhyloNode parsed tree with {leaf_count} terminal clades.")

    # KEGG
    kegg_node = BioKEGGPathwayNode()
    path_json, p_name = kegg_node.run("hsa04110")
    p_info = json.loads(path_json)
    assert len(p_info["genes"]) > 0
    print(f"[PASS] BioKEGGPathwayNode retrieved pathway: {p_name} ({len(p_info['genes'])} genes).")


def test_11_pyhmmer_search_and_manifest(setup_dirs):
    """Test PyhmmerSearchNode execution and manifest export."""
    from pyhmmer import easel, plan7

    alphabet = easel.Alphabet.amino()
    seq1 = easel.TextSequence(name=b"seq1", sequence="MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIA")
    builder = plan7.Builder(alphabet)
    background = plan7.Background(alphabet)
    hmm, _, _ = builder.build(seq1.digitize(alphabet), background)

    hmm_path = OUTPUT_BASE / "test_model.hmm"
    with open(hmm_path, "wb") as hf:
        hmm.write(hf)

    prot_path = OUTPUT_BASE / "test_proteins.fasta"
    prot_path.write_text(">seq1 albumin-like\nMKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIA\n", encoding="utf-8")

    runner = CondaCommandRunner()
    hmmer_node = PyhmmerSearchNode()
    hmmer_out = OUTPUT_BASE / "hmmer_run"
    hits_json, num_hits = hmmer_node.run(str(prot_path), str(hmm_path), output_dir=str(hmmer_out), runner=runner)
    assert num_hits >= 1
    hits = json.loads(hits_json)
    assert hits[0]["target_name"] == "seq1"
    assert (hmmer_out / "run_manifest.sh").exists()
    assert (hmmer_out / "run_manifest.json").exists()
    print(f"\n[PASS] PyhmmerSearchNode identified {num_hits} domain hits and generated manifest.")


def test_12_motifs_and_tfbs_scanning(setup_dirs):
    """Test BiopythonBioMotifsNode and PyTFBSScanNode."""
    # Motifs PWM
    motif_node = BiopythonBioMotifsNode()
    pwm_json, consensus = motif_node.run("TACGAT\nTATGAT\nTACAAT\nTATAAT")
    pwm_dict = json.loads(pwm_json)
    assert "A" in pwm_dict
    assert len(consensus) == 6
    print(f"\n[PASS] BiopythonBioMotifsNode generated consensus: {consensus}.")

    # TFBS Scanner
    tfbs_node = PyTFBSScanNode()
    tsv_out, site_count = tfbs_node.run("TTTTGTAAACAGGGG", jaspar_id="MA0139.1")
    assert "MA0139.1" in tsv_out
    assert site_count >= 1
    print(f"[PASS] PyTFBSScanNode found {site_count} binding sites for MA0139.1.")


def test_13_file_not_found_validation(setup_dirs):
    """Test that nodes check for input file existence and raise proper FileNotFoundError."""
    missing_file = str(OUTPUT_BASE / "non_existent_file.xyz")

    with pytest.raises(FileNotFoundError):
        BiopythonSeqIONode().run(missing_file)

    with pytest.raises(FileNotFoundError):
        BiopythonBioPDBNode().run(missing_file)

    with pytest.raises(FileNotFoundError):
        PyfaidxIndexNode().run(missing_file)

    with pytest.raises(FileNotFoundError):
        BiopythonPhyloNode().run(missing_file)

    with pytest.raises(FileNotFoundError):
        PyhmmerSearchNode().run(missing_file, str(OUTPUT_BASE / "non_existent.hmm"))

    with pytest.raises(FileNotFoundError):
        BiotiteStructureAlignNode().run(missing_file, missing_file)

    with pytest.raises(FileNotFoundError):
        BiopythonBlastNode().run("MKWV", subject_file=missing_file)

    print("\n[PASS] All nodes properly validated input file existence and raised FileNotFoundError.")
