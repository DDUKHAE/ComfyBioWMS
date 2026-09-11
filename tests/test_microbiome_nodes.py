import json
import os
import sys
from pathlib import Path
import pytest
import pandas as pd
import numpy as np

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from bioflow.runtime.command_runner import DryRunCommandRunner
from nodes.class_1 import (
    BiomFormatNode,
    Ete3TreeParserNode,
    FastUniFracDistanceNode,
    ScikitBioDiversityNode,
)
from nodes.class_2 import (
    AmrFinderPlusNode,
    BaktaAnnotationNode,
    CheckVQualityNode,
    Dada2AmpliconNode,
    GeNomadViromeNode,
    Humann3PathwayNode,
    Metaphlan4ProfileNode,
    NextstrainAugurNode,
    ProkkaAnnotationNode,
    VirSorter2Node,
)

OUTPUT_DIR = root_dir / "results" / "test_microbiome_nodes"


@pytest.fixture(scope="module", autouse=True)
def setup_test_environment():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------------------
# 1. Missing File Error Tests
# ------------------------------------------------------------------------------

def test_missing_input_file_errors():
    """Verify that all nodes raise FileNotFoundError when inputs are missing."""
    bogus = str(OUTPUT_DIR / "non_existent_file.xyz")

    with pytest.raises(FileNotFoundError):
        ScikitBioDiversityNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        BiomFormatNode().run(bogus, str(OUTPUT_DIR / "out.csv"))

    with pytest.raises(FileNotFoundError):
        Ete3TreeParserNode().run(bogus, str(OUTPUT_DIR / "out.nwk"))

    with pytest.raises(FileNotFoundError):
        FastUniFracDistanceNode().run(bogus, bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        Dada2AmpliconNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        Humann3PathwayNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        Metaphlan4ProfileNode().run(bogus, str(OUTPUT_DIR / "out.txt"))

    with pytest.raises(FileNotFoundError):
        GeNomadViromeNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        VirSorter2Node().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        CheckVQualityNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        ProkkaAnnotationNode().run(bogus, str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        BaktaAnnotationNode().run(bogus, str(OUTPUT_DIR), str(OUTPUT_DIR))

    with pytest.raises(FileNotFoundError):
        AmrFinderPlusNode().run(bogus, str(OUTPUT_DIR / "out.tsv"))

    with pytest.raises(FileNotFoundError):
        NextstrainAugurNode().run(bogus, bogus, str(OUTPUT_DIR / "out.json"))


# ------------------------------------------------------------------------------
# 2. In-Memory Ecological & Phylogenetic Tests
# ------------------------------------------------------------------------------

def test_biom_format_node():
    """Test BiomFormatNode on a real BIOM 1.0 JSON table."""
    biom_file = OUTPUT_DIR / "test_table.biom"
    out_csv = OUTPUT_DIR / "converted_table.csv"

    biom_data = {
        "id": "test_biom",
        "format": "1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "shape": [3, 2],
        "rows": [{"id": "OTU_1"}, {"id": "OTU_2"}, {"id": "OTU_3"}],
        "columns": [{"id": "Sample_A"}, {"id": "Sample_B"}],
        "matrix_type": "dense",
        "matrix_element_type": "float",
        "data": [
            [120.0, 45.0],
            [30.0, 80.0],
            [0.0, 15.0],
        ]
    }
    biom_file.write_text(json.dumps(biom_data), encoding="utf-8")

    node = BiomFormatNode()
    res_csv, num_obs, num_samples = node.run(str(biom_file), str(out_csv))

    assert Path(res_csv).is_file()
    assert num_obs == 3
    assert num_samples == 2

    df = pd.read_csv(res_csv, index_col=0)
    assert df.shape == (3, 2)
    assert list(df.columns) == ["Sample_A", "Sample_B"]
    assert df.loc["OTU_1", "Sample_A"] == 120.0


def test_ete3_tree_parser_node():
    """Test Ete3TreeParserNode parsing and leaf counting on a real Newick string."""
    nwk_file = OUTPUT_DIR / "test_tree.nwk"
    annotated_nwk_file = OUTPUT_DIR / "annotated_tree.nwk"

    nwk_content = "((OTU_1:0.12,OTU_2:0.34):0.56,(OTU_3:0.78,OTU_4:0.90):0.21);"
    nwk_file.write_text(nwk_content, encoding="utf-8")

    node = Ete3TreeParserNode()
    res_file, leaf_count = node.run(str(nwk_file), str(annotated_nwk_file))

    assert Path(res_file).is_file()
    assert leaf_count == 4
    annotated_text = Path(res_file).read_text()
    assert "OTU_1" in annotated_text
    assert "OTU_4" in annotated_text


def test_fast_unifrac_distance_node():
    """Test FastUniFracDistanceNode computation of weighted and unweighted UniFrac."""
    biom_file = OUTPUT_DIR / "unifrac_test.biom"
    tree_file = OUTPUT_DIR / "unifrac_tree.nwk"
    out_dir = OUTPUT_DIR / "unifrac_out"

    biom_data = {
        "id": "unifrac_biom",
        "format": "1.0.0",
        "type": "OTU table",
        "shape": [4, 2],
        "rows": [{"id": "OTU_1"}, {"id": "OTU_2"}, {"id": "OTU_3"}, {"id": "OTU_4"}],
        "columns": [{"id": "Sample_1"}, {"id": "Sample_2"}],
        "matrix_type": "dense",
        "matrix_element_type": "float",
        "data": [
            [50.0, 0.0],
            [50.0, 0.0],
            [0.0, 60.0],
            [0.0, 40.0],
        ]
    }
    biom_file.write_text(json.dumps(biom_data), encoding="utf-8")
    tree_file.write_text("((OTU_1:0.2,OTU_2:0.2):0.5,(OTU_3:0.2,OTU_4:0.2):0.5);\n", encoding="utf-8")

    node = FastUniFracDistanceNode()
    unw_csv, w_csv = node.run(str(biom_file), str(tree_file), str(out_dir))

    assert Path(unw_csv).is_file()
    assert Path(w_csv).is_file()

    df_unw = pd.read_csv(unw_csv, index_col=0)
    df_w = pd.read_csv(w_csv, index_col=0)

    assert df_unw.shape == (2, 2)
    assert df_w.shape == (2, 2)
    # Diagonal should be 0.0
    assert df_unw.loc["Sample_1", "Sample_1"] == 0.0
    assert df_w.loc["Sample_1", "Sample_1"] == 0.0
    # Completely disjoint clades should have distance > 0
    assert df_unw.loc["Sample_1", "Sample_2"] > 0.5
    assert df_w.loc["Sample_1", "Sample_2"] > 0.5


def test_scikit_bio_diversity_node():
    """Test ScikitBioDiversityNode calculating alpha diversity (Shannon, Simpson) and beta diversity (PCoA)."""
    abund_csv = OUTPUT_DIR / "abundance.csv"
    out_dir = OUTPUT_DIR / "diversity_out"

    abund_df = pd.DataFrame(
        [
            [100, 20, 5, 0],
            [90, 30, 10, 0],
            [0, 5, 80, 120],
            [0, 2, 90, 110],
        ],
        index=["S1", "S2", "S3", "S4"],
        columns=["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D"],
    )
    abund_df.index.name = "sample_id"
    abund_df.to_csv(abund_csv)

    node = ScikitBioDiversityNode()
    alpha_csv, pcoa_csv = node.run(str(abund_csv), str(out_dir))

    assert Path(alpha_csv).is_file()
    assert Path(pcoa_csv).is_file()

    alpha = pd.read_csv(alpha_csv)
    assert "shannon" in alpha.columns
    assert "simpson" in alpha.columns
    assert len(alpha) == 4
    # Real computed indices must be positive
    assert all(alpha["shannon"] > 0)
    assert all(alpha["simpson"] > 0)

    pcoa = pd.read_csv(pcoa_csv, index_col=0)
    assert "PC1" in pcoa.columns
    assert len(pcoa) == 4
    # S1 and S2 should have similar PC1 signs, opposite of S3 and S4
    s1_pc1 = pcoa.loc["S1", "PC1"]
    s2_pc1 = pcoa.loc["S2", "PC1"]
    s3_pc1 = pcoa.loc["S3", "PC1"]
    s4_pc1 = pcoa.loc["S4", "PC1"]
    assert (s1_pc1 * s2_pc1) > 0
    assert (s3_pc1 * s4_pc1) > 0
    assert (s1_pc1 * s3_pc1) < 0


# ------------------------------------------------------------------------------
# 3. CLI Subprocess Nodes Execution & Manifest Tests
# ------------------------------------------------------------------------------

def test_dada2_amplicon_node():
    """Verify Dada2AmpliconNode dispatches Rscript and generates run manifests."""
    fastq_dir = OUTPUT_DIR / "16s_fastqs"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    (fastq_dir / "sample1_R1.fastq").write_text("@read1\nAGCT\n+\nIIII\n")

    out_dir = OUTPUT_DIR / "dada2_out"
    runner = DryRunCommandRunner()

    node = Dada2AmpliconNode()
    asv_seqtab, asv_fasta = node.run(str(fastq_dir), str(out_dir), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "Rscript"
    assert "dada2_pipeline.R" in cmd.argv[1]
    assert "--fastq_dir" in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()
    assert (out_dir / "run_manifest.json").is_file()


def test_humann3_pathway_node():
    """Verify Humann3PathwayNode dispatches humann -i ... -o ... and generates manifests."""
    fq_file = OUTPUT_DIR / "sample.fastq.gz"
    fq_file.write_text("FAKE_FASTQ_GZ")
    out_dir = OUTPUT_DIR / "humann_out"
    runner = DryRunCommandRunner()

    node = Humann3PathwayNode()
    gf, pa = node.run(str(fq_file), str(out_dir), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "humann"
    assert "-i" in cmd.argv
    assert "-o" in cmd.argv
    assert str(fq_file) in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_metaphlan4_profile_node():
    """Verify Metaphlan4ProfileNode dispatches metaphlan --input_type ... -o ..."""
    fq_file = OUTPUT_DIR / "sample.fq.gz"
    fq_file.write_text("FAKE_FASTQ_GZ")
    out_profile = OUTPUT_DIR / "metaphlan" / "profile.txt"
    runner = DryRunCommandRunner()

    node = Metaphlan4ProfileNode()
    (res_profile,) = node.run(str(fq_file), str(out_profile), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "metaphlan"
    assert "--input_type" in cmd.argv
    assert "fastq" in cmd.argv
    assert "-o" in cmd.argv
    assert (out_profile.parent / "run_manifest.sh").is_file()


def test_genomad_virome_node():
    """Verify GeNomadViromeNode dispatches genomad end-to-end."""
    fa_file = OUTPUT_DIR / "contigs.fa"
    fa_file.write_text(">contig1\nAGCTAGCTAGCT\n")
    out_dir = OUTPUT_DIR / "genomad_out"
    runner = DryRunCommandRunner()

    node = GeNomadViromeNode()
    v_sum, p_sum = node.run(str(fa_file), str(out_dir), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "genomad"
    assert cmd.argv[1] == "end-to-end"
    assert str(fa_file) in cmd.argv
    assert str(out_dir) in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_virsorter2_node():
    """Verify VirSorter2Node dispatches virsorter run -w ... -i ..."""
    fa_file = OUTPUT_DIR / "contigs.fa"
    fa_file.write_text(">contig1\nAGCTAGCTAGCT\n")
    out_dir = OUTPUT_DIR / "virsorter2_out"
    runner = DryRunCommandRunner()

    node = VirSorter2Node()
    (v_boundary,) = node.run(str(fa_file), str(out_dir), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "virsorter"
    assert cmd.argv[1] == "run"
    assert "-w" in cmd.argv
    assert "-i" in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_checkv_quality_node():
    """Verify CheckVQualityNode dispatches checkv end_to_end."""
    fa_file = OUTPUT_DIR / "viruses.fa"
    fa_file.write_text(">virus1\nAGCTAGCTAGCT\n")
    out_dir = OUTPUT_DIR / "checkv_out"
    runner = DryRunCommandRunner()

    node = CheckVQualityNode()
    q_sum, c_fa = node.run(str(fa_file), str(out_dir), runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "checkv"
    assert cmd.argv[1] == "end_to_end"
    assert str(fa_file) in cmd.argv
    assert str(out_dir) in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_prokka_annotation_node():
    """Verify ProkkaAnnotationNode dispatches prokka --outdir ... --prefix ..."""
    fa_file = OUTPUT_DIR / "contigs.fa"
    fa_file.write_text(">contig1\nAGCTAGCTAGCT\n")
    out_dir = OUTPUT_DIR / "prokka_out"
    runner = DryRunCommandRunner()

    node = ProkkaAnnotationNode()
    gff, faa = node.run(str(fa_file), str(out_dir), prefix="MYPROKKA", runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "prokka"
    assert "--outdir" in cmd.argv
    assert "--prefix" in cmd.argv
    assert "MYPROKKA" in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_bakta_annotation_node():
    """Verify BaktaAnnotationNode dispatches bakta --db ... --output ... --prefix ..."""
    fa_file = OUTPUT_DIR / "contigs.fa"
    fa_file.write_text(">contig1\nAGCTAGCTAGCT\n")
    out_dir = OUTPUT_DIR / "bakta_out"
    db_dir = OUTPUT_DIR / "bakta_db"
    runner = DryRunCommandRunner()

    node = BaktaAnnotationNode()
    gff, summary = node.run(str(fa_file), str(out_dir), str(db_dir), prefix="MYBAKTA", runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "bakta"
    assert "--db" in cmd.argv
    assert "--output" in cmd.argv
    assert "--prefix" in cmd.argv
    assert "MYBAKTA" in cmd.argv
    assert (out_dir / "run_manifest.sh").is_file()


def test_amrfinder_plus_node():
    """Verify AmrFinderPlusNode dispatches amrfinder with correct input flags."""
    faa_file = OUTPUT_DIR / "proteins.faa"
    faa_file.write_text(">prot1\nMKVIL\n")
    out_tsv = OUTPUT_DIR / "amr" / "report.tsv"
    runner = DryRunCommandRunner()

    node = AmrFinderPlusNode()
    (res_tsv,) = node.run(str(faa_file), str(out_tsv), organism="Salmonella", runner=runner)

    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "amrfinder"
    assert "-p" in cmd.argv
    assert "-O" in cmd.argv
    assert "Salmonella" in cmd.argv
    assert (out_tsv.parent / "run_manifest.sh").is_file()


def test_nextstrain_augur_node():
    """Verify NextstrainAugurNode dispatches augur tree and augur export."""
    aln_file = OUTPUT_DIR / "aligned.fa"
    aln_file.write_text(">seq1\nACGT\n>seq2\nACGA\n")
    meta_file = OUTPUT_DIR / "metadata.tsv"
    meta_file.write_text("strain\tdate\nseq1\t2026-01-01\nseq2\t2026-01-02\n")
    out_json = OUTPUT_DIR / "auspice" / "tree.json"
    runner = DryRunCommandRunner()

    node = NextstrainAugurNode()
    (res_json,) = node.run(str(aln_file), str(meta_file), str(out_json), runner=runner)

    assert len(runner.commands) == 2
    cmd1 = runner.commands[0]
    assert cmd1.argv[0] == "augur"
    assert cmd1.argv[1] == "tree"
    assert "--alignment" in cmd1.argv

    cmd2 = runner.commands[1]
    assert cmd2.argv[0] == "augur"
    assert cmd2.argv[1] == "export"
    assert (out_json.parent / "run_manifest.sh").is_file()
