import json
import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from bioflow.runtime.command_runner import DryRunCommandRunner
from nodes.class_1 import (
    Cyvcf2VariantNode,
    MappyAlignNode,
    PybedtoolsIntervalNode,
    PyfastxIndexNode,
    PysamAnalysisNode,
)
from nodes.class_2 import (
    Bowtie2AlignNode,
    CuteSvNode,
    DeepVariantCallNode,
    FlyeAssembleNode,
    HifiasmAssembleNode,
    MedakaConsensusNode,
    MosdepthCoverageNode,
    RaconPolishNode,
    SeqKitToolNode,
    Sniffles2SvNode,
)

TEST_DIR = root_dir / "results" / "test_genomics_longread"


@pytest.fixture(scope="module")
def setup_test_files():
    TEST_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Reference FASTA
    ref_fa = TEST_DIR / "ref.fa"
    ref_seq = "ACGTACGTTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 10
    ref_fa.write_text(f">chr1\n{ref_seq}\n", encoding="utf-8")

    # 2. Reads FASTQ
    reads_fq = TEST_DIR / "reads.fq"
    read_seq = ref_seq[20:120]
    reads_fq.write_text(f"@read_01\n{read_seq}\n+\n{'I'*len(read_seq)}\n", encoding="utf-8")

    # 3. BED files
    bed_a = TEST_DIR / "peaks.bed"
    bed_a.write_text("chr1\t100\t200\tpeak1\nchr1\t300\t400\tpeak2\n", encoding="utf-8")
    bed_b = TEST_DIR / "promoters.bed"
    bed_b.write_text("chr1\t150\t250\tprom1\nchr1\t500\t600\tprom2\n", encoding="utf-8")

    # 4. VCF file
    vcf_path = TEST_DIR / "sample.vcf"
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t50\t.\tA\tG\t30\tPASS\t.\n"  # Transition (A->G)
        "chr1\t60\t.\tC\tT\t30\tPASS\t.\n"  # Transition (C->T)
        "chr1\t70\t.\tA\tC\t30\tPASS\t.\n"  # Transversion (A->C)
        "chr1\t80\t.\tA\tAT\t30\tPASS\t.\n" # Indel
    )
    vcf_path.write_text(vcf_content, encoding="utf-8")

    # 5. Bowtie2 index dummy file
    bt2_idx = TEST_DIR / "bowtie2_idx" / "genome.1.bt2"
    bt2_idx.parent.mkdir(parents=True, exist_ok=True)
    bt2_idx.write_text("dummy bowtie2 index", encoding="utf-8")

    # 6. Dedicated valid BAM file using pysam
    bam_path = TEST_DIR / "sample.sorted.bam"
    import pysam
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'chr1'}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        a = pysam.AlignedSegment()
        a.query_name = "read_1"
        a.query_sequence = "ACGT"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigar = ((0, 4),)
        outf.write(a)

    return {
        "ref_fa": ref_fa,
        "reads_fq": reads_fq,
        "bed_a": bed_a,
        "bed_b": bed_b,
        "vcf_path": vcf_path,
        "bt2_prefix": bt2_idx.parent / "genome",
        "bam_path": bam_path,
    }


# =============================================================================
# Python Library Functional Tests
# =============================================================================

def test_01_pysam_analysis_real(setup_test_files):
    node = PysamAnalysisNode()
    bam_file = str(setup_test_files["bam_path"])
    stats_json, mapped_reads, mean_mapq = node.run(bam_file)

    assert mapped_reads >= 1
    assert mean_mapq > 0.0
    parsed = json.loads(stats_json)
    assert parsed["bam"] == bam_file
    assert parsed["mapped_reads"] == mapped_reads
    assert parsed["mean_mapq"] == mean_mapq

    with pytest.raises(FileNotFoundError):
        node.run("non_existent_file.bam")


def test_02_cyvcf2_variant_real(setup_test_files):
    node = Cyvcf2VariantNode()
    vcf_file = str(setup_test_files["vcf_path"])
    summary_json, snvs, indels = node.run(vcf_file)

    assert snvs == 3
    assert indels == 1
    parsed = json.loads(summary_json)
    assert parsed["transitions"] == 2
    assert parsed["transversions"] == 1
    assert parsed["ti_tv_ratio"] == 2.0

    with pytest.raises(FileNotFoundError):
        node.run("non_existent_file.vcf")


def test_03_pybedtools_interval_real(setup_test_files):
    node = PybedtoolsIntervalNode()
    bed_a = str(setup_test_files["bed_a"])
    bed_b = str(setup_test_files["bed_b"])
    out_bed = str(TEST_DIR / "intersected.bed")

    out_path, count = node.run(bed_a, bed_b, out_bed)
    assert Path(out_path).is_file()
    assert count == 1  # 100-200 overlaps with 150-250 (at 150-200)

    with pytest.raises(FileNotFoundError):
        node.run("non_existent_a.bed", bed_b)
    with pytest.raises(FileNotFoundError):
        node.run(bed_a, "non_existent_b.bed")


def test_04_pyfastx_index_real(setup_test_files):
    node = PyfastxIndexNode()
    fastx_file = str(setup_test_files["reads_fq"])
    summary_json, total_seqs, n50 = node.run(fastx_file)

    assert total_seqs == 1
    assert n50 == 100.0
    parsed = json.loads(summary_json)
    assert parsed["total_seqs"] == 1
    assert parsed["total_bases"] == 100
    assert parsed["n50"] == 100.0
    assert parsed["gc_content"] > 0.0

    with pytest.raises(FileNotFoundError):
        node.run("non_existent_fastx.fq")


def test_05_mappy_align_real(setup_test_files):
    node = MappyAlignNode()
    ref_fa = str(setup_test_files["ref_fa"])
    reads_fq = str(setup_test_files["reads_fq"])
    out_paf = str(TEST_DIR / "minimap.paf")

    paf_path, aligned_reads = node.run(ref_fa, reads_fq, out_paf, preset="sr")
    assert Path(paf_path).is_file()
    assert aligned_reads >= 1

    with pytest.raises(FileNotFoundError):
        node.run("non_existent_ref.fa", reads_fq)
    with pytest.raises(FileNotFoundError):
        node.run(ref_fa, "non_existent_reads.fq")


# =============================================================================
# CLI Dispatch Functional Tests
# =============================================================================

def test_06_seqkit_tool_node_dispatch(setup_test_files):
    node = SeqKitToolNode()
    runner = DryRunCommandRunner()
    in_fq = str(setup_test_files["reads_fq"])
    out_fq = str(TEST_DIR / "seqkit" / "filtered.fq")

    res = node.run(in_fq, out_fq, min_length=50, extra_command="--quiet", runner=runner)
    assert res == (out_fq,)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "seqkit"
    assert cmd.argv[1] == "seq"
    assert "-m" in cmd.argv
    assert str(in_fq) in cmd.argv
    assert (Path(out_fq).parent / "run_manifest.sh").is_file()
    assert (Path(out_fq).parent / "run_manifest.json").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.fq", out_fq, runner=runner)


def test_07_bowtie2_align_node_dispatch(setup_test_files):
    node = Bowtie2AlignNode()
    runner = DryRunCommandRunner()
    idx = str(setup_test_files["bt2_prefix"])
    r1 = str(setup_test_files["reads_fq"])
    out_bam = str(TEST_DIR / "bowtie2" / "aligned.bam")

    res = node.run(idx, r1, out_bam, runner=runner)
    assert res == (out_bam,)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "bowtie2"
    assert "-x" in cmd.argv
    assert idx in cmd.argv
    assert (Path(out_bam).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing_idx", r1, out_bam, runner=runner)
    with pytest.raises(FileNotFoundError):
        node.run(idx, "missing_r1.fq", out_bam, runner=runner)


def test_08_mosdepth_coverage_node_dispatch(setup_test_files):
    node = MosdepthCoverageNode()
    runner = DryRunCommandRunner()
    bam = str(setup_test_files["bam_path"])
    prefix = str(TEST_DIR / "mosdepth" / "sample_cov")

    res = node.run(bam, prefix, runner=runner)
    assert res == (f"{prefix}.regions.bed.gz", f"{prefix}.global.dist.txt")
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "mosdepth"
    assert prefix in cmd.argv
    assert bam in cmd.argv
    assert (Path(prefix).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.bam", prefix, runner=runner)


def test_09_sniffles2_sv_node_dispatch(setup_test_files):
    node = Sniffles2SvNode()
    runner = DryRunCommandRunner()
    bam = str(setup_test_files["bam_path"])
    vcf = str(TEST_DIR / "sniffles" / "sv.vcf")

    res = node.run(bam, vcf, runner=runner)
    assert res == (vcf,)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "sniffles"
    assert "-i" in cmd.argv
    assert "-v" in cmd.argv
    assert (Path(vcf).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.bam", vcf, runner=runner)


def test_10_cutesv_node_dispatch(setup_test_files):
    node = CuteSvNode()
    runner = DryRunCommandRunner()
    bam = str(setup_test_files["bam_path"])
    ref = str(setup_test_files["ref_fa"])
    vcf = str(TEST_DIR / "cutesv" / "sv.vcf")

    res = node.run(bam, ref, vcf, runner=runner)
    assert res == (vcf,)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "cuteSV"
    assert bam in cmd.argv
    assert ref in cmd.argv
    assert (Path(vcf).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.bam", ref, vcf, runner=runner)
    with pytest.raises(FileNotFoundError):
        node.run(bam, "missing_ref.fa", vcf, runner=runner)


def test_11_flye_assemble_node_dispatch(setup_test_files):
    node = FlyeAssembleNode()
    runner = DryRunCommandRunner()
    reads = str(setup_test_files["reads_fq"])
    out_dir = str(TEST_DIR / "flye_out")

    fasta, gfa = node.run(reads, out_dir, genome_size="10m", runner=runner)
    assert fasta.endswith("assembly.fasta")
    assert gfa.endswith("assembly_graph.gfa")
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "flye"
    assert "--nano-raw" in cmd.argv
    assert "--genome-size" in cmd.argv
    assert (Path(out_dir) / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing_reads.fq", out_dir, runner=runner)


def test_12_hifiasm_assemble_node_dispatch(setup_test_files):
    node = HifiasmAssembleNode()
    runner = DryRunCommandRunner()
    reads = str(setup_test_files["reads_fq"])
    prefix = str(TEST_DIR / "hifiasm_out" / "asm")

    p_gfa, a_gfa = node.run(reads, prefix, runner=runner)
    assert p_gfa.endswith(".p_ctg.gfa")
    assert a_gfa.endswith(".a_ctg.gfa")
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "hifiasm"
    assert "-o" in cmd.argv
    assert (Path(prefix).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing_reads.fq", prefix, runner=runner)


def test_13_racon_polish_node_dispatch(setup_test_files):
    node = RaconPolishNode()
    runner = DryRunCommandRunner()
    draft = str(setup_test_files["ref_fa"])
    paf = str(setup_test_files["bed_a"])  # any existing file as mock PAF
    reads = str(setup_test_files["reads_fq"])
    out_fa = str(TEST_DIR / "racon" / "polished.fa")

    res = node.run(draft, paf, reads, out_fa, runner=runner)
    assert res == (out_fa,)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "racon"
    assert reads in cmd.argv
    assert paf in cmd.argv
    assert draft in cmd.argv
    assert (Path(out_fa).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.fa", paf, reads, out_fa, runner=runner)


def test_14_medaka_consensus_node_dispatch(setup_test_files):
    node = MedakaConsensusNode()
    runner = DryRunCommandRunner()
    draft = str(setup_test_files["ref_fa"])
    reads = str(setup_test_files["reads_fq"])
    out_dir = str(TEST_DIR / "medaka_out")

    res = node.run(draft, reads, out_dir, model="r1041_e82_400bps_sup_v4.2.0", runner=runner)
    assert res == (str(Path(out_dir) / "consensus.fasta"),)
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "medaka_consensus"
    assert "-i" in cmd.argv
    assert "-d" in cmd.argv
    assert "-o" in cmd.argv
    assert "-m" in cmd.argv
    assert (Path(out_dir) / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.fa", reads, out_dir, runner=runner)


def test_15_deepvariant_call_node_dispatch(setup_test_files):
    node = DeepVariantCallNode()
    runner = DryRunCommandRunner()
    bam = str(setup_test_files["bam_path"])
    ref = str(setup_test_files["ref_fa"])
    vcf = str(TEST_DIR / "deepvariant" / "calls.vcf.gz")

    vcf_out, gvcf_out = node.run(bam, ref, vcf, model_type="WGS", runner=runner)
    assert vcf_out == vcf
    assert gvcf_out == str(TEST_DIR / "deepvariant" / "calls.g.vcf.gz")
    assert len(runner.commands) == 1
    cmd = runner.commands[0]
    assert cmd.argv[0] == "run_deepvariant"
    assert "--model_type=WGS" in cmd.argv
    assert f"--ref={ref}" in cmd.argv
    assert f"--reads={bam}" in cmd.argv
    assert (Path(vcf).parent / "run_manifest.sh").is_file()

    with pytest.raises(FileNotFoundError):
        node.run("missing.bam", ref, vcf, runner=runner)
    with pytest.raises(FileNotFoundError):
        node.run(bam, "missing_ref.fa", vcf, runner=runner)
