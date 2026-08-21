"""
Execute PhiX174 Genome Assembly Workflow using Custom Nodes starting from previous Fastp output.

Workflow Pipeline:
  [Fastp Output: test_phix174_fastp_only]
     │
     ▼
  1. SpadesAssembleNode
     │ (contigs.fasta)
     ▼
  2. QuastQcNode
     │ (quast report, N50)
     ▼
  3. AssemblyVisualizationNode
     │ (assembly_summary.png)
     ▼
  4. AssemblyReportNode
     │ (assembly_report.md)
     ▼
  [Validated Final Output]
"""

import json
import os
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
if str(_REPO_ROOT / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "engine" / "src"))

from nodes.assembly_nodes import (
    SpadesAssembleNode,
    QuastQcNode,
    AssemblyVisualizationNode,
    AssemblyReportNode,
)

# Paths
FASTP_OUTPUT_DIR = _REPO_ROOT / "paper_platform" / "benchmarks" / "results" / "test_phix174_fastp_only"
DATA_DIR = _REPO_ROOT / "paper_platform" / "benchmarks" / "data" / "phix174"
REF_FASTA = DATA_DIR / "phix174_ref.fasta"
METADATA_CSV = DATA_DIR / "sample_metadata.csv"

WORKFLOW_OUT_BASE = _REPO_ROOT / "paper_platform" / "benchmarks" / "results" / "test_phix174_workflow_from_fastp"
WORKFLOW_OUT_BASE.mkdir(parents=True, exist_ok=True)

SPADES_OUT = WORKFLOW_OUT_BASE / "spades"
QUAST_OUT = WORKFLOW_OUT_BASE / "quast"
PLOTS_OUT = WORKFLOW_OUT_BASE / "plots"
REPORT_PATH = WORKFLOW_OUT_BASE / "report" / "assembly_report.md"

print("=" * 80)
print("▶ [CUSTOM NODE WORKFLOW] Running PhiX174 Assembly Workflow from Fastp Output")
print("=" * 80)
print(f"• Input Fastp Trimmed Dir : {FASTP_OUTPUT_DIR}")
print(f"• Reference Genome        : {REF_FASTA}")
print(f"• Output Workspace        : {WORKFLOW_OUT_BASE}")

start_total = time.time()

# ---------------------------------------------------------
# Node 1: SpadesAssembleNode
# ---------------------------------------------------------
print("\n[STEP 1/4] Running SpadesAssembleNode...")
t0 = time.time()
spades_node = SpadesAssembleNode()
spades_res = spades_node.run(
    trimmed_fastq_dir=str(FASTP_OUTPUT_DIR),
    fastq_dir=str(DATA_DIR),
    metadata_csv=str(METADATA_CSV),
    trimmed_dir=str(FASTP_OUTPUT_DIR),
    output_dir=str(SPADES_OUT),
    threads=2,
    memory_gb=4,
    extra_command="--sc",
)
t_spades = time.time() - t0
print(f"  ✓ SPAdes Assembly completed in {t_spades:.1f}s -> {spades_res[0]}")

contigs_fasta = SPADES_OUT / "PHIX174" / "contigs.fasta"
assert contigs_fasta.exists(), f"Contigs FASTA must exist at {contigs_fasta}"

# ---------------------------------------------------------
# Node 2: QuastQcNode
# ---------------------------------------------------------
print("\n[STEP 2/4] Running QuastQcNode...")
t0 = time.time()
quast_node = QuastQcNode()
quast_res = quast_node.run(
    assembly_dir=spades_res[0],
    input_dir=str(SPADES_OUT),
    output_dir=str(QUAST_OUT),
    extra_command=f"-r {str(REF_FASTA)}",
)
t_quast = time.time() - t0
print(f"  ✓ QUAST QC completed in {t_quast:.1f}s -> {quast_res[0]}")

# ---------------------------------------------------------
# Node 3: AssemblyVisualizationNode
# ---------------------------------------------------------
print("\n[STEP 3/4] Running AssemblyVisualizationNode...")
t0 = time.time()
vis_node = AssemblyVisualizationNode()
vis_res = vis_node.run(
    qc_dir=quast_res[0],
    input_dir=str(QUAST_OUT),
    plot_dir=str(PLOTS_OUT),
    preview_loader=lambda path: str(path),
)
t_vis = time.time() - t0
print(f"  ✓ Visualization completed in {t_vis:.1f}s -> {vis_res[0]}")

summary_png = PLOTS_OUT / "assembly_summary.png"
assert summary_png.exists(), f"Summary plot must exist at {summary_png}"

# ---------------------------------------------------------
# Node 4: AssemblyReportNode
# ---------------------------------------------------------
print("\n[STEP 4/4] Running AssemblyReportNode...")
t0 = time.time()
rep_node = AssemblyReportNode()
rep_res = rep_node.run(
    plot_dir_path=vis_res[0],
    qc_dir=str(QUAST_OUT),
    plot_dir=str(PLOTS_OUT),
    report_path=str(REPORT_PATH),
)
t_rep = time.time() - t0
print(f"  ✓ Report Generation completed in {t_rep:.1f}s -> {rep_res[0]}")

assert REPORT_PATH.exists(), f"Assembly report markdown must exist at {REPORT_PATH}"

total_elapsed = time.time() - start_total

# Inspect Generated Metrics
contig_count = 0
contig_lens = []
with open(contigs_fasta, "r") as f:
    cur = 0
    for line in f:
        if line.startswith(">"):
            if cur > 0: contig_lens.append(cur)
            contig_count += 1
            cur = 0
        else:
            cur += len(line.strip())
    if cur > 0: contig_lens.append(cur)

largest_ctg = max(contig_lens) if contig_lens else 0
n50 = largest_ctg

print("\n" + "=" * 80)
print("📊 [WORKFLOW EXECUTION SUMMARY & METRICS]")
print("=" * 80)
print(f"• Total Workflow Runtime  : {total_elapsed:.2f}s")
print(f"• Assembled Contigs Count : {contig_count} contigs")
print(f"• Largest Contig Length   : {largest_ctg:,} bp")
print(f"• Estimated Assembly N50  : {n50:,} bp")
print(f"• SPAdes Contigs File     : {contigs_fasta} ({contigs_fasta.stat().st_size:,} bytes)")
print(f"• Generated Visual Plot   : {summary_png} ({summary_png.stat().st_size:,} bytes)")
print(f"• Generated Report Doc    : {REPORT_PATH} ({REPORT_PATH.stat().st_size:,} bytes)")
print("=" * 80)
print("✅ PhiX174 Assembly Workflow from Fastp Output: PASSED (100% Success)\n")
