#!/usr/bin/env python
"""Generate ComfyUI workflow JSON for legacy and reduced case-study sets.

Widget order and port names are read from the live registry, so a node signature
change here is a regenerate away instead of a hand-edited JSON drift.

Usage:
  python engine/scripts/generate_case_study_workflows.py
  python engine/scripts/generate_case_study_workflows.py paper
  python engine/scripts/generate_case_study_workflows.py final
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from nodes.registry import NODE_CLASS_MAPPINGS  # noqa: E402

# Terminal driver: ComfyUI core PreviewAny / PreviewImage nodes, OUTPUT_NODE=True.
# Note: All ComfyBIO nodes now also declare OUTPUT_NODE = True for standalone execution.
MODE = sys.argv[1] if len(sys.argv) > 1 else "nfcore"
SOURCE = "paper" if MODE == "paper" else "nfcore"

# Data dirs per --source. The paper sets carry ground truth (mock composition,
# GIAB truth regions, the exact PhiX genome) that the nf-core test data cannot
# provide; CS-1/CS-5 have no paper variant here (different graph shape / same file).
_DIRS = {
    "nfcore": ("data/nf_core_taxprofiler", "data/nf_core_sarek", "data/nf_core_bacass"),
    "paper": ("data/paper_zymo_d6300", "data/paper_giab_hg001", "data/paper_phix174"),
}[SOURCE]

PREVIEW_ANY = "PreviewAny"
PREVIEW_IMAGE = "PreviewImage"


def _widgets(node_type):
    """(names, defaults) for every widget of a node, in ComfyUI's order."""
    spec = NODE_CLASS_MAPPINGS[node_type].INPUT_TYPES()
    names, defaults = [], []
    for section in ("required", "optional"):
        for name, decl in spec.get(section, {}).items():
            kind = decl[0]
            cfg = decl[1] if len(decl) > 1 else {}
            names.append(name)
            defaults.append(cfg.get("default", kind[0] if isinstance(kind, list) else ""))
    return names, defaults


def _outputs(node_type):
    cls = NODE_CLASS_MAPPINGS[node_type]
    return list(getattr(cls, "RETURN_NAMES", ())), list(getattr(cls, "RETURN_TYPES", ()))


def build(steps, title, benchmark=None):
    """steps: list of dicts {id, type, w={widget: value}, i={input: (src_id, out_name)}}."""
    steps_by_id = {s["id"]: s for s in steps}
    depth, order = {}, []
    for s in steps:
        d = max((depth[src] + 1 for src, _ in s.get("i", {}).values()), default=0)
        depth[s["id"]] = d
        order.append(s["id"])
    rows = {}
    nodes, links, link_id = [], [], 0
    outs = {}  # node_id -> {slot_index: [link_ids]}

    for s in steps:
        nid, ntype = s["id"], s["type"]
        d = depth[nid]
        row = rows.get(d, 0)
        rows[d] = row + 1

        if ntype in (PREVIEW_ANY, PREVIEW_IMAGE):
            names, defaults, out_names, out_types = [], [], [], []
        else:
            names, defaults = _widgets(ntype)
            out_names, out_types = _outputs(ntype)

        values = list(defaults)
        for k, v in s.get("w", {}).items():
            values[names.index(k)] = v

        inputs = []
        for in_name, (src_id, out_name) in s.get("i", {}).items():
            src_names, src_types = _outputs(steps_by_id[src_id]["type"])
            slot = src_names.index(out_name)
            link_id += 1
            links.append([link_id, src_id, slot, nid, len(inputs), src_types[slot]])
            outs.setdefault(src_id, {}).setdefault(slot, []).append(link_id)
            entry = {"name": in_name, "type": src_types[slot], "link": link_id}
            if in_name in names:  # widget-backed input
                entry["widget"] = {"name": in_name}
            inputs.append(entry)

        nodes.append({
            "id": nid,
            "type": ntype,
            "pos": [80 + d * 400, 80 + row * 260],
            "size": [340, 26 + 24 * max(len(values), 1)],
            "flags": {},
            "order": order.index(nid),
            "mode": 0,
            "inputs": inputs,
            "outputs": [
                {"name": n, "type": t, "slot_index": i, "links": []}
                for i, (n, t) in enumerate(zip(out_names, out_types))
            ],
            "properties": {"Node name for S&R": ntype},
            "widgets_values": values,
            "title": s.get("title", ntype),
        })

    by_id = {n["id"]: n for n in nodes}
    for src_id, slots in outs.items():
        for slot, lids in slots.items():
            by_id[src_id]["outputs"][slot]["links"] = lids

    extra = {"comfybio_case_study": title}
    if benchmark:
        extra["benchmark"] = benchmark

    return {
        "last_node_id": max(depth),
        "last_link_id": link_id,
        "version": 0.4,
        "config": {},
        "extra": extra,
        "groups": [],
        "nodes": nodes,
        "links": links,
    }


def validate(graph):
    """Fail loudly on the mistakes a hand-written workflow JSON actually makes."""
    by_id = {n["id"]: n for n in graph["nodes"]}
    seen_output_node = False
    for n in graph["nodes"]:
        if n["type"] in (PREVIEW_ANY, PREVIEW_IMAGE):
            seen_output_node = True
            continue
        cls = NODE_CLASS_MAPPINGS[n["type"]]  # KeyError => unregistered node
        names, _ = _widgets(n["type"])
        assert len(n["widgets_values"]) == len(names), f"{n['type']}: widget count drift"
        required = set(cls.INPUT_TYPES().get("required", {}))
        supplied = {i["name"] for i in n["inputs"]} | {
            w for w, v in zip(names, n["widgets_values"]) if v not in ("", None)
        }
        missing = required - supplied
        assert not missing, f"{n['type']}: required input unset: {sorted(missing)}"
    assert seen_output_node, "graph has no OUTPUT_NODE; ComfyUI would execute nothing"

    for lid, src, slot, dst, dslot, ltype in graph["links"]:
        assert src in by_id and dst in by_id, f"link {lid}: dangling node"
        assert by_id[src]["outputs"][slot]["type"] == ltype, f"link {lid}: type drift"
        assert lid in by_id[src]["outputs"][slot]["links"], f"link {lid}: not registered on source"
        assert by_id[dst]["inputs"][dslot]["link"] == lid, f"link {lid}: input slot mismatch"
        assert by_id[dst]["order"] > by_id[src]["order"], f"link {lid}: not topologically ordered"


# --------------------------------------------------------------------------
# CS-1  nf-core/rnaseq — quantification-route concordance (Salmon / STAR+Salmon / HISAT2+StringTie)
# --------------------------------------------------------------------------
RNA = "data/nf_core_rnaseq"
CS1 = [
    {"id": 1, "type": "FastQC", "title": "FastQC (raw reads)", "w": {"input_file": RNA}},
    {"id": 2, "type": "TrimGalore", "title": "TrimGalore (batch trimming)",
     "w": {"reads_fwd": RNA}},
    {"id": 3, "type": "SalmonIndex", "w": {"transcripts_fasta": f"{RNA}/transcriptome.fasta"}},
    {"id": 4, "type": "InferStrandedness",
     "i": {"salmon_index_dir": (3, "salmon_index_dir"),
           "reads_fwd": (2, "trimmed_reads_fwd"), "reads_rev": (2, "trimmed_reads_rev")}},
    # route A - Salmon selective alignment (nf-core --pseudo_aligner salmon)
    {"id": 5, "type": "SalmonQuantReads", "title": "Route A: Salmon (selective alignment)",
     "i": {"salmon_index_dir": (3, "salmon_index_dir"), "reads_fwd": (2, "trimmed_reads_fwd"),
           "reads_rev": (2, "trimmed_reads_rev"), "strandedness": (4, "strandedness")}},
    # route B - STAR + Salmon alignment mode (nf-core default: star_salmon)
    {"id": 6, "type": "STARGenomeGenerate", "title": "STAR index",
     "w": {"genome_fasta": f"{RNA}/genome.fasta", "gtf_file": f"{RNA}/genes.gtf",
           "genome_sa_index_nbases": 9}},
    {"id": 7, "type": "STARAlignReads",
     "i": {"star_index_dir": (6, "star_index_dir"), "reads_fwd": (2, "trimmed_reads_fwd"),
           "reads_rev": (2, "trimmed_reads_rev")}},
    {"id": 8, "type": "SalmonQuantAlignment", "title": "Route B: STAR + Salmon (alignment mode)",
     "w": {"transcripts_fasta": f"{RNA}/transcriptome.fasta"},
     "i": {"transcriptome_bam": (7, "transcriptome_bam"), "strandedness": (4, "strandedness")}},
    # route C - HISAT2 + StringTie
    {"id": 9, "type": "HISAT2Build", "w": {"reference_fasta": f"{RNA}/genome.fasta"}},
    {"id": 10, "type": "HISAT2Align",
     "i": {"hisat2_index_prefix": (9, "hisat2_index_prefix"),
           "reads_fwd": (2, "trimmed_reads_fwd"), "reads_rev": (2, "trimmed_reads_rev"),
           "strandedness": (4, "strandedness")}},
    {"id": 11, "type": "SamtoolsSort", "i": {"input_alignment": (10, "aligned_sam")}},
    {"id": 12, "type": "StringTie", "title": "Route C: HISAT2 + StringTie",
     "w": {"guide_gtf": f"{RNA}/genes.gtf"}, "i": {"bam_file": (11, "sorted_bam")}},
    # aggregate + QC
    # tximport directly consumes batch quant.sf from Salmon route A
    {"id": 13, "type": "Tximport", "title": "tximport (4-sample gene matrix)",
     "i": {"quant_files": (5, "quant_sf")}},
    {"id": 14, "type": "DESeq2SampleQC", "w": {"metadata_tsv": f"{RNA}/sample_metadata.csv"},
     "i": {"counts_tsv": (13, "gene_counts_tsv")}},
    {"id": 15, "type": "MultiQC", "w": {"analysis_dir": "results/cs1_rnaseq",
                                        "report_title": "CS-1 nf-core/rnaseq route concordance"},
     "i": {"extra_scan_dirs": (5, "salmon_out_dir")}},
    # DEG analysis + Volcano visualization
    {"id": 16, "type": "DESeq2", "title": "DESeq2 (Differential Expression)",
     "w": {"sample_metadata_csv": f"{RNA}/sample_metadata.csv"},
     "i": {"count_matrix_csv": (13, "gene_counts_tsv")}},
    {"id": 17, "type": "VolcanoPlot", "title": "Volcano Plot (DEGs)",
     "i": {"deg_results_csv": (16, "deseq2_results_csv")}},
    {"id": 18, "type": PREVIEW_ANY, "title": "Route A gene counts",
     "i": {"source": (13, "gene_counts_tsv")}},
    {"id": 19, "type": PREVIEW_ANY, "title": "Route B quant.sf", "i": {"source": (8, "quant_sf")}},
    {"id": 20, "type": PREVIEW_ANY, "title": "Route C abundances",
     "i": {"source": (12, "gene_abundances_tsv")}},
    {"id": 21, "type": PREVIEW_ANY, "title": "Sample QC summary", "i": {"source": (14, "summary_json")}},
    {"id": 22, "type": PREVIEW_ANY, "title": "DESeq2 DEG summary", "i": {"source": (16, "deg_summary_json")}},
    {"id": 23, "type": PREVIEW_ANY, "title": "Volcano Plot PNG", "i": {"source": (17, "plot_path")}},
    {"id": 24, "type": PREVIEW_ANY, "title": "MultiQC report", "i": {"source": (15, "multiqc_html")}},
    {"id": 25, "type": PREVIEW_ANY, "title": "FastQC report", "i": {"source": (1, "html_report")}},
]

# --------------------------------------------------------------------------
# CS-2  nf-core/taxprofiler — Zymo D6300 mock, 3 profilers
# --------------------------------------------------------------------------
TAXP = _DIRS[0]
# Reference DBs ship with the nf-core test set and are reused by both sources.
TAXDB = "data/nf_core_taxprofiler"
CS2 = [
    {"id": 1, "type": "Fastp", "w": {"read1": f"{TAXP}/reads_R1.fastq.gz",
                                     "read2": f"{TAXP}/reads_R2.fastq.gz",
                                     "detect_adapter_for_pe": True}},
    {"id": 4, "type": "Kraken2Classify", "title": "Profiler 1: Kraken2",
     "w": {"kraken2_db": f"{TAXDB}/kraken2_db", "confidence": 0.05},
     "i": {"reads_fwd": (1, "read1"), "reads_rev": (1, "read2")}},
    {"id": 5, "type": "Bracken", "w": {"bracken_db": f"{TAXDB}/bracken_db", "read_length": 150},
     "i": {"kraken_report": (4, "kraken_report_txt")}},
    {"id": 6, "type": "MetaPhlAn", "title": "Profiler 2: MetaPhlAn",
     "i": {"reads_fwd": (1, "read1"), "reads_rev": (1, "read2")}},
    {"id": 7, "type": "SylphProfile", "title": "Profiler 3: sylph",
     "w": {"sylph_db": f"{TAXDB}/sylph_db.syldb"},
     "i": {"reads_fwd": (1, "read1"), "reads_rev": (1, "read2")}},
    {"id": 8, "type": "MicrobiomeStackedBar",
     "w": {"figure_title": "Bracken species composition", "top_taxa": 10},
     "i": {"abundance_csv": (5, "bracken_abundance_tsv")}},
    {"id": 9, "type": PREVIEW_IMAGE, "i": {"images": (8, "preview_image")}},
    {"id": 10, "type": PREVIEW_ANY, "title": "Bracken abundance",
     "i": {"source": (5, "bracken_abundance_tsv")}},
    {"id": 11, "type": PREVIEW_ANY, "title": "MetaPhlAn profile",
     "i": {"source": (6, "profiled_metagenome_tsv")}},
    {"id": 12, "type": PREVIEW_ANY, "title": "sylph profile", "i": {"source": (7, "sylph_profile_tsv")}},
    {"id": 14, "type": PREVIEW_ANY, "title": "fastp report", "i": {"source": (1, "json_report")}},
]

# --------------------------------------------------------------------------
# CS-3  nf-core/sarek (bcftools branch) — GIAB NA12878 / HG001
# --------------------------------------------------------------------------
GIAB = _DIRS[1]
CS3 = [
    {"id": 1, "type": "Fastp", "w": {"read1": f"{GIAB}/sample1_R1.fastq.gz",
                                     "read2": f"{GIAB}/sample1_R2.fastq.gz"}},
    {"id": 2, "type": "BwaMem2Index", "w": {"reference_fasta": f"{GIAB}/genome.fasta"}},
    {"id": 3, "type": "BwaMem2Align",
     "w": {"threads": 4, "read_group": r"@RG\tID:HG001\tSM:NA12878\tPL:ILLUMINA"},
     "i": {"indexed_reference": (2, "indexed_reference"),
           "read1": (1, "read1"), "read2": (1, "read2")}},
    {"id": 4, "type": "SamtoolsSort", "w": {"threads": 4}, "i": {"input_alignment": (3, "sam_file")}},
    {"id": 5, "type": "PicardMarkDuplicates", "i": {"input_bam": (4, "sorted_bam")}},
    {"id": 6, "type": "SamtoolsIndex", "i": {"input_bam": (5, "marked_bam")}},
    {"id": 7, "type": "Mosdepth", "title": "Coverage QC", "i": {"bam_file": (5, "marked_bam")}},
    {"id": 8, "type": "BcftoolsMpileup", "w": {"reference_fasta": f"{GIAB}/genome.fasta",
                                               "max_depth": 250, "min_base_quality": 20},
     "i": {"input_bam": (5, "marked_bam")}},
    {"id": 9, "type": "BcftoolsCall", "i": {"input_bcf": (8, "pileup_bcf")}},
    {"id": 10, "type": "BcftoolsFilter", "w": {"exclude": "QUAL<30 || DP<10"},
     "i": {"input_vcf": (9, "called_vcf")}},
    {"id": 11, "type": "Cyvcf2Stats", "i": {"vcf_file": (10, "filtered_vcf")}},
    {"id": 12, "type": "PybedtoolsIntersect", "title": "Calls inside GIAB high-conf regions",
     "w": {"bed_b": f"{GIAB}/giab_truth.bed"},
     "i": {"bed_a": (10, "filtered_vcf")}},
    {"id": 13, "type": PREVIEW_ANY, "title": "Variant stats", "i": {"source": (11, "summary_json")}},
    {"id": 14, "type": PREVIEW_ANY, "title": "High-conf call count",
     "i": {"source": (12, "feature_count")}},
    {"id": 15, "type": PREVIEW_ANY, "title": "Coverage summary", "i": {"source": (7, "summary_txt")}},
    {"id": 16, "type": PREVIEW_ANY, "title": "Duplicate metrics", "i": {"source": (5, "metrics_txt")}},
    {"id": 17, "type": PREVIEW_ANY, "title": "BAM index", "i": {"source": (6, "alignment_index")}},
]

# The GIAB truth-region intersect only means anything against the paper truth
# set; the nf-core test reference has no benchmark BED.
if SOURCE != "paper":
    CS3 = [step for step in CS3 if step["id"] not in (12, 14)]


# --------------------------------------------------------------------------
# CS-4  nf-core/bacass — PhiX174 (NC_001422.1) assembly + annotation
# --------------------------------------------------------------------------
BACASS = _DIRS[2]
CS4 = [
    {"id": 1, "type": "Fastp", "w": {"read1": f"{BACASS}/reads_R1.fastq.gz",
                                     "read2": f"{BACASS}/reads_R2.fastq.gz",
                                     "detect_adapter_for_pe": True, "correction": True}},
    {"id": 2, "type": "FastQC", "i": {"input_file": (1, "read1")}},
    {"id": 3, "type": "SeqKitStats", "i": {"sequence_file": (1, "read1")}},
    {"id": 4, "type": "Spades", "w": {"careful": True, "threads": 8, "memory_gb": 8},
     "i": {"read1": (1, "read1"), "read2": (1, "read2")}},
    {"id": 5, "type": "Quast",
     "title": "QUAST vs NC_001422.1 (5,386 bp)" if SOURCE == "paper" else "QUAST (reference-free)",
     "w": {"min_contig": 200,
           **({"reference": f"{BACASS}/phix174_ref.fasta"} if SOURCE == "paper" else {})},
     "i": {"assembly_fasta": (4, "contigs")}},
    {"id": 6, "type": "Prokka", "w": {"kingdom": "Bacteria"},
     "i": {"genome_fasta": (4, "contigs")}},
    {"id": 7, "type": "SeqKitStats", "title": "Assembly stats",
     "i": {"sequence_file": (4, "contigs")}},
    {"id": 8, "type": PREVIEW_ANY, "title": "QUAST report.tsv", "i": {"source": (5, "report_tsv")}},
    {"id": 9, "type": PREVIEW_ANY, "title": "Prokka GFF", "i": {"source": (6, "annotation_gff")}},
    {"id": 10, "type": PREVIEW_ANY, "title": "Assembly stats", "i": {"source": (7, "stats_tsv")}},
    {"id": 11, "type": PREVIEW_ANY, "title": "Read stats", "i": {"source": (3, "stats_tsv")}},
    {"id": 12, "type": PREVIEW_ANY, "title": "FastQC report", "i": {"source": (2, "html_report")}},
]

# --------------------------------------------------------------------------
# CS-5  Scanpy best-practices — pancreas endocrinogenesis (E15.5)
# --------------------------------------------------------------------------
PANC = "data/Pancreas/endocrinogenesis_day15.h5ad"
CS5 = [
    {"id": 1, "type": "AnnDataInspect", "w": {"h5ad_file": PANC}},
    {"id": 2, "type": "ScanpyQC", "w": {"input_h5ad": PANC, "min_genes": 200, "min_cells": 3,
                                        "max_percent_mito": 20.0}},
    {"id": 3, "type": "ScanpyNormalize", "w": {"target_sum": 10000.0, "n_top_genes": 2000},
     "i": {"input_h5ad": (2, "filtered_h5ad")}},
    {"id": 4, "type": "ScanpyCluster", "w": {"n_pcs": 30, "resolution": 0.5},
     "i": {"input_h5ad": (3, "normalized_h5ad")}},
    {"id": 5, "type": PREVIEW_ANY, "title": "Dataset summary", "i": {"source": (1, "summary_json")}},
    {"id": 6, "type": PREVIEW_ANY, "title": "QC metrics", "i": {"source": (2, "qc_metrics_json")}},
    {"id": 7, "type": PREVIEW_ANY, "title": "UMAP PNG path", "i": {"source": (4, "umap_plot_png")}},
    {"id": 8, "type": PREVIEW_ANY, "title": "Clustered h5ad", "i": {"source": (4, "clustered_h5ad")}},
]

CASES = [
    ("cs1_nfcore_rnaseq_route_concordance", CS1, "CS-1 nf-core/rnaseq quantification-route concordance"),
    ("cs2_nfcore_taxprofiler_multiprofiler", CS2, "CS-2 nf-core/taxprofiler multi-profiler comparison"),
    ("cs3_nfcore_sarek_variant_calling", CS3, "CS-3 nf-core/sarek bcftools variant-calling branch"),
    ("cs4_nfcore_bacass_assembly", CS4, "CS-4 nf-core/bacass assembly + annotation"),
    ("cs5_scanpy_pancreas_endocrinogenesis", CS5, "CS-5 Scanpy best-practices scRNA-seq"),
]

# Public-facing reduced set. These labels intentionally distinguish a protocol
# smoke test from a biological reproduction claim.
ZYMO_PAPER_DIR = "data/paper_zymo_d6300"
ZYMO_METAPHLAN = [
    {"id": 1, "type": "Fastp", "w": {"read1": f"{ZYMO_PAPER_DIR}/reads_R1.fastq.gz",
                                     "read2": f"{ZYMO_PAPER_DIR}/reads_R2.fastq.gz",
                                     "detect_adapter_for_pe": True}},
    {"id": 6, "type": "MetaPhlAn", "title": "Profiler 2: MetaPhlAn",
     "i": {"reads_fwd": (1, "read1"), "reads_rev": (1, "read2")}},
    {"id": 11, "type": PREVIEW_ANY, "title": "MetaPhlAn profile",
     "i": {"source": (6, "profiled_metagenome_tsv")}},
    {"id": 14, "type": PREVIEW_ANY, "title": "fastp report", "i": {"source": (1, "json_report")}},
]
FINAL_CASES = [
    (
        "final_01_pancreas_scanpy_e2e",
        CS5,
        "Pancreas E15.5 dataset-backed Scanpy technical E2E",
        {
            "case_id": "CS-R1",
            "claim_level": "technical_e2e",
            "dataset_status": "verified_exact_file",
            "dataset_sha256": "9e3e459eca00ba06b496ec80def32941b5b2889918720e3e7aa6ffb811fbe7c6",
            "dimensions": {"cells": 3696, "genes": 27998},
            "limitation": "This graph validates Scanpy execution; it is not an scVelo result reproduction.",
        },
    ),
    (
        "final_02_zymo_metaphlan_e2e",
        ZYMO_METAPHLAN,
        "ZymoBIOMICS D6300 read-backed MetaPhlAn E2E",
        {
            "case_id": "CS-E1",
            "claim_level": "e2e",
            "dataset_status": "verified_accession_first_100000_pairs",
            "source_accession": "SRR12324253",
            "read_pairs": 100000,
            "dataset_sha256": {
                "reads_R1.fastq.gz": "ee8597b4ad203b66521fb3fb88960520bfdb87d2609107477c28a7f3fe58b043",
                "reads_R2.fastq.gz": "0efc31fb5835eafb1b7a657ef6be6ca8222ab34c217c213ffef6b60fc72e7d84",
            },
            "runtime_database_required": True,
            "limitation": "A compatible full MetaPhlAn database is not bundled with the downloaded reads.",
        },
    ),
    (
        "final_03_gse110004_chri_smoke",
        CS1,
        "GSE110004 chromosome-I RNA-seq workflow smoke test",
        {
            "case_id": "CS-S1",
            "claim_level": "smoke",
            "dataset_status": "verified_nfcore_fixture",
            "source_revision": "nf-core/test-datasets rnaseq@e07c1b1",
            "runs": ["SRR6357070", "SRR6357072", "SRR6357076", "SRR6357077"],
            "read_pairs_per_run": 50000,
            "reference_scope": "Saccharomyces cerevisiae chromosome I only",
            "limitation": "Subsampled reads and a chromosome-I reference cannot reproduce the full paper result.",
        },
    ),
]

FINAL_STATUS = [
    {"case_id": "CS-R1", "status": "ready_technical_e2e", "workflow": "final_01_pancreas_scanpy_e2e.json"},
    {"case_id": "CS-E1", "status": "requires_runtime_database", "workflow": "final_02_zymo_metaphlan_e2e.json"},
    {
        "case_id": "CS-E2",
        "status": "blocked",
        "workflow": None,
        "reason": "Downloaded SRR5458066 is a human-skin 16S run, not a PhiX174 assembly read set.",
    },
    {"case_id": "CS-S1", "status": "ready_smoke", "workflow": "final_03_gse110004_chri_smoke.json"},
]

if __name__ == "__main__":
    out_dir = ROOT / "workflows"
    if MODE == "final":
        cases = FINAL_CASES
    else:
        cases = [c for c in CASES if c[0][:3] in ("cs2", "cs3", "cs4")] if SOURCE == "paper" else CASES
    for case in cases:
        name, steps, title, *metadata = case
        if SOURCE == "paper":
            name, title = f"{name}_paper", f"{title} (paper ground truth)"
        graph = build(steps, title, metadata[0] if metadata else None)
        validate(graph)
        dest = out_dir / f"{name}.json"
        dest.write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")
        print(f"[OK] {dest.relative_to(ROOT)}  ({len(graph['nodes'])} nodes, {len(graph['links'])} links)")
    if MODE == "final":
        dest = out_dir / "final_case_study_status.json"
        dest.write_text(json.dumps(FINAL_STATUS, indent=2) + "\n", encoding="utf-8")
        print(f"[OK] {dest.relative_to(ROOT)}")
