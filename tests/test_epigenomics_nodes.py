import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from bioflow.runtime.command_runner import DryRunCommandRunner
from nodes.class_1 import ScikitFusionNode
from nodes.class_2 import (
    ChromosightLoopNode,
    CoolerMatrixNode,
    CooltoolsTadNode,
    Crispresso2ScreenNode,
    DeepToolsProfileNode,
    HomerMotifNode,
    MageckScreenNode,
    MemeSuiteNode,
    MethyldackelNode,
    PyGenrichNode,
    SeacrPeakNode,
    TobiasFootprintNode,
)

TEST_BASE = root_dir / "results" / "test_epigenomics_nodes"


@pytest.fixture(scope="module")
def setup_epigenomics_env():
    TEST_BASE.mkdir(parents=True, exist_ok=True)

    # Create dummy mock inputs for CLI execution verification
    bam_file = TEST_BASE / "test.bam"
    bam_file.touch()

    bed_file = TEST_BASE / "peaks.bed"
    bed_file.write_text("chr1\t1000\t2000\tpeak1\n", encoding="utf-8")

    fasta_file = TEST_BASE / "genome.fa"
    fasta_file.write_text(">chr1\nACGTACGTACGTACGT\n", encoding="utf-8")

    fastq_file = TEST_BASE / "amplicon_R1.fq"
    fastq_file.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")

    counts_file = TEST_BASE / "counts.txt"
    counts_file.write_text("sgRNA\tGene\ttreat_1\tctrl_1\nsg1\tG1\t100\t10\n", encoding="utf-8")

    bedgraph_file = TEST_BASE / "treat.bedgraph"
    bedgraph_file.write_text("chr1\t1000\t2000\t5.5\n", encoding="utf-8")

    ctrl_bedgraph = TEST_BASE / "ctrl.bedgraph"
    ctrl_bedgraph.write_text("chr1\t1000\t2000\t1.0\n", encoding="utf-8")

    chrom_sizes = TEST_BASE / "chrom.sizes"
    chrom_sizes.write_text("chr1\t5000000\n", encoding="utf-8")

    pairs_file = TEST_BASE / "contacts.pairs"
    pairs_file.write_text("read1\tchr1\t1000\tchr1\t5000\n", encoding="utf-8")

    cool_file = TEST_BASE / "matrix.cool"
    cool_file.touch()

    # Multi-omics matrices for ScikitFusionNode
    samples = [f"sample_{i}" for i in range(10)]
    df_exp = pd.DataFrame(
        np.random.RandomState(42).randn(10, 20),
        index=samples,
        columns=[f"gene_{j}" for j in range(20)],
    )
    df_meth = pd.DataFrame(
        np.random.RandomState(43).randn(10, 15),
        index=samples,
        columns=[f"cpg_{k}" for k in range(15)],
    )
    exp_csv = TEST_BASE / "gene_expression.csv"
    meth_csv = TEST_BASE / "methylation.csv"
    df_exp.to_csv(exp_csv)
    df_meth.to_csv(meth_csv)

    return {
        "bam_file": str(bam_file),
        "bed_file": str(bed_file),
        "fasta_file": str(fasta_file),
        "fastq_file": str(fastq_file),
        "counts_file": str(counts_file),
        "bedgraph_file": str(bedgraph_file),
        "ctrl_bedgraph": str(ctrl_bedgraph),
        "chrom_sizes": str(chrom_sizes),
        "pairs_file": str(pairs_file),
        "cool_file": str(cool_file),
        "exp_csv": str(exp_csv),
        "meth_csv": str(meth_csv),
    }


def test_missing_files_raise_filenotfound(setup_epigenomics_env):
    """Verify that all nodes raise FileNotFoundError when inputs do not exist."""
    fake = str(TEST_BASE / "non_existent_file.xyz")

    with pytest.raises(FileNotFoundError):
        DeepToolsProfileNode().run(fake, str(TEST_BASE / "out.bw"))

    with pytest.raises(FileNotFoundError):
        TobiasFootprintNode().run(fake, fake, fake, str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        PyGenrichNode().run(fake, str(TEST_BASE / "out.bed"))

    with pytest.raises(FileNotFoundError):
        Crispresso2ScreenNode().run(fake, "ACGT", "CGAT", str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        MageckScreenNode().run(fake, "t1", "c1", str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        ScikitFusionNode().run(fake, fake)

    with pytest.raises(FileNotFoundError):
        SeacrPeakNode().run(fake, "0.01", "stringent", str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        HomerMotifNode().run(fake, "hg38", str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        MemeSuiteNode().run(fake, str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        MethyldackelNode().run(fake, fake, str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        CoolerMatrixNode().run(fake, fake, str(TEST_BASE / "out.cool"))

    with pytest.raises(FileNotFoundError):
        CooltoolsTadNode().run(fake, 10000, str(TEST_BASE / "out"))

    with pytest.raises(FileNotFoundError):
        ChromosightLoopNode().run(fake, str(TEST_BASE / "out.bedpe"))


def test_deeptools_profile_node(setup_epigenomics_env):
    """Test DeepToolsProfileNode CLI command dispatch and manifest generation."""
    runner = DryRunCommandRunner()
    node = DeepToolsProfileNode()
    out_bw = str(TEST_BASE / "deeptools" / "signal.bw")
    bw, mat = node.run(
        bam_file=setup_epigenomics_env["bam_file"],
        output_bigwig=out_bw,
        regions_bed=setup_epigenomics_env["bed_file"],
        runner=runner,
    )
    assert bw == out_bw
    assert mat.endswith("_matrix.gz")
    assert len(runner.commands) == 2
    assert runner.commands[0].argv[0] == "bamCoverage"
    assert runner.commands[1].argv[0] == "computeMatrix"
    assert (Path(out_bw).parent / "run_manifest.sh").exists()
    assert (Path(out_bw).parent / "run_manifest.json").exists()


def test_tobias_footprint_node(setup_epigenomics_env):
    """Test TobiasFootprintNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = TobiasFootprintNode()
    out_dir = str(TEST_BASE / "tobias")
    bw, tsv = node.run(
        bam_file=setup_epigenomics_env["bam_file"],
        peaks_bed=setup_epigenomics_env["bed_file"],
        ref_fasta=setup_epigenomics_env["fasta_file"],
        output_dir=out_dir,
        runner=runner,
    )
    assert bw.endswith("footprints.bw")
    assert tsv.endswith("bindetect_results.txt")
    assert len(runner.commands) == 2
    assert runner.commands[0].argv[:2] == ["TOBIAS", "ATACorrect"]
    assert runner.commands[1].argv[:2] == ["TOBIAS", "ScoreBigwig"]
    assert (Path(out_dir) / "run_manifest.sh").exists()


def test_pygenrich_node(setup_epigenomics_env):
    """Test PyGenrichNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = PyGenrichNode()
    out_bed = str(TEST_BASE / "genrich" / "peaks.bed")
    bed, bg = node.run(
        bam_file=setup_epigenomics_env["bam_file"],
        output_bed=out_bed,
        runner=runner,
    )
    assert bed == out_bed
    assert bg.endswith(".bedgraph")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[0] == "Genrich"
    assert "-t" in runner.commands[0].argv
    assert "-o" in runner.commands[0].argv
    assert (Path(out_bed).parent / "run_manifest.sh").exists()


def test_crispresso2_screen_node(setup_epigenomics_env):
    """Test Crispresso2ScreenNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = Crispresso2ScreenNode()
    out_dir = str(TEST_BASE / "crispresso")
    stats, report_dir = node.run(
        read1_fastq=setup_epigenomics_env["fastq_file"],
        amplicon_seq="ATGCGATCGATCGATCGATCGATCGATAG",
        sgrna_seq="CGATCGATCGATCGATCGAT",
        output_dir=out_dir,
        runner=runner,
    )
    assert stats.endswith("CRISPResso_quantification_of_editing_frequency.txt")
    assert report_dir == out_dir
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[0] == "CRISPResso"
    assert "--fastq_r1" in runner.commands[0].argv
    assert (Path(out_dir) / "run_manifest.sh").exists()


def test_mageck_screen_node(setup_epigenomics_env):
    """Test MageckScreenNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = MageckScreenNode()
    prefix = str(TEST_BASE / "mageck" / "screen")
    gene_s, sgrna_s = node.run(
        count_table_tsv=setup_epigenomics_env["counts_file"],
        treatment_samples="treat_1",
        control_samples="ctrl_1",
        output_prefix=prefix,
        runner=runner,
    )
    assert gene_s.endswith(".gene_summary.txt")
    assert sgrna_s.endswith(".sgrna_summary.txt")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[:2] == ["mageck", "test"]
    assert "-k" in runner.commands[0].argv
    assert "-t" in runner.commands[0].argv
    assert "-c" in runner.commands[0].argv
    assert "-n" in runner.commands[0].argv
    assert (Path(prefix).parent / "run_manifest.sh").exists()


def test_seacr_peak_node(setup_epigenomics_env):
    """Test SeacrPeakNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = SeacrPeakNode()
    prefix = str(TEST_BASE / "seacr" / "peaks")
    (peaks_bed,) = node.run(
        treatment_bedgraph=setup_epigenomics_env["bedgraph_file"],
        control_bedgraph=setup_epigenomics_env["ctrl_bedgraph"],
        mode="stringent",
        output_prefix=prefix,
        runner=runner,
    )
    assert peaks_bed == f"{prefix}.stringent.bed"
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[0] == "SEACR_1.3.sh"
    assert (Path(prefix).parent / "run_manifest.sh").exists()


def test_homer_motif_node(setup_epigenomics_env):
    """Test HomerMotifNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = HomerMotifNode()
    out_dir = str(TEST_BASE / "homer")
    motifs_dir, html = node.run(
        peaks_bed=setup_epigenomics_env["bed_file"],
        genome="hg38",
        output_dir=out_dir,
        runner=runner,
    )
    assert motifs_dir == out_dir
    assert html.endswith("knownResults.html")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[0] == "findMotifsGenome.pl"
    assert (Path(out_dir) / "run_manifest.sh").exists()


def test_meme_suite_node(setup_epigenomics_env):
    """Test MemeSuiteNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = MemeSuiteNode()
    out_dir = str(TEST_BASE / "meme")
    xml, tsv = node.run(
        sequences_fasta=setup_epigenomics_env["fasta_file"],
        output_dir=out_dir,
        runner=runner,
    )
    assert xml.endswith("meme.xml")
    assert tsv.endswith("fimo.tsv")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[0] == "meme"
    assert (Path(out_dir) / "run_manifest.sh").exists()


def test_methyldackel_node(setup_epigenomics_env):
    """Test MethyldackelNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = MethyldackelNode()
    prefix = str(TEST_BASE / "methyldackel" / "sample")
    cpg_bg, summary = node.run(
        ref_fasta=setup_epigenomics_env["fasta_file"],
        wgbs_bam=setup_epigenomics_env["bam_file"],
        output_prefix=prefix,
        runner=runner,
    )
    assert cpg_bg.endswith("_CpG.bedGraph")
    assert summary.endswith("_summary.txt")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[:2] == ["MethylDackel", "extract"]
    assert "-o" in runner.commands[0].argv
    assert (Path(prefix).parent / "run_manifest.sh").exists()


def test_cooler_matrix_node(setup_epigenomics_env):
    """Test CoolerMatrixNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = CoolerMatrixNode()
    out_cool = str(TEST_BASE / "cooler" / "matrix.cool")
    cool, mcool = node.run(
        pairs_file=setup_epigenomics_env["pairs_file"],
        chrom_sizes=setup_epigenomics_env["chrom_sizes"],
        output_cool=out_cool,
        resolution=50000,
        runner=runner,
    )
    assert cool == out_cool
    assert mcool.endswith(".mcool")
    assert len(runner.commands) == 2
    assert runner.commands[0].argv[:2] == ["cooler", "cload"]
    assert runner.commands[1].argv[:2] == ["cooler", "zoomify"]
    assert (Path(out_cool).parent / "run_manifest.sh").exists()


def test_cooltools_tad_node(setup_epigenomics_env):
    """Test CooltoolsTadNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = CooltoolsTadNode()
    prefix = str(TEST_BASE / "cooltools" / "tad")
    ins_tsv, bnd_bed = node.run(
        cool_matrix=setup_epigenomics_env["cool_file"],
        window_size=100000,
        output_prefix=prefix,
        runner=runner,
    )
    assert ins_tsv.endswith("_insulation.tsv")
    assert bnd_bed.endswith("_boundaries.bed")
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[:2] == ["cooltools", "insulation"]
    assert (Path(prefix).parent / "run_manifest.sh").exists()


def test_chromosight_loop_node(setup_epigenomics_env):
    """Test ChromosightLoopNode CLI command dispatch."""
    runner = DryRunCommandRunner()
    node = ChromosightLoopNode()
    out_bedpe = str(TEST_BASE / "chromosight" / "loops.bedpe")
    (loops,) = node.run(
        cool_matrix=setup_epigenomics_env["cool_file"],
        output_bedpe=out_bedpe,
        runner=runner,
    )
    assert loops == out_bedpe
    assert len(runner.commands) == 1
    assert runner.commands[0].argv[:2] == ["chromosight", "detect"]
    assert (Path(out_bedpe).parent / "run_manifest.sh").exists()


def test_scikit_fusion_node_real_execution(setup_epigenomics_env):
    """Test ScikitFusionNode performing real matrix factorization on multi-omics data."""
    node = ScikitFusionNode()
    out_npy = str(TEST_BASE / "scikit_fusion" / "fused_matrix.npy")
    npy_path, factors_json = node.run(
        matrix1_csv=setup_epigenomics_env["exp_csv"],
        matrix2_csv=setup_epigenomics_env["meth_csv"],
        rank=5,
        output_npy=out_npy,
    )
    assert Path(npy_path).exists()
    fused_arr = np.load(npy_path)
    assert fused_arr.shape == (10, 5)

    meta = json.loads(factors_json)
    assert meta["rank"] == 5
    assert meta["num_samples"] == 10
    assert meta["matrix1_features"] == 20
    assert meta["matrix2_features"] == 15
    assert meta["total_features"] == 35
    assert 0.0 <= meta["reconstruction_error"] <= 1.0
    assert 0.0 <= meta["variance_explained_ratio"] <= 1.0
    assert len(meta["singular_values"]) == 5
