"""
Direct Test for AssemblyFastpTrimNode on PhiX174 Biological FASTQ Dataset.
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

from nodes.assembly_nodes import AssemblyFastpTrimNode

PHIX_DATA_DIR = _REPO_ROOT / "paper_platform" / "benchmarks" / "data" / "phix174"
OUTPUT_DIR = _REPO_ROOT / "paper_platform" / "benchmarks" / "results" / "test_phix174_fastp_only"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("▶ [CUSTOM NODE TEST] Testing AssemblyFastpTrimNode on PhiX174 Dataset")
print("=" * 80)
print(f"• Input Data Dir : {PHIX_DATA_DIR}")
print(f"• Target Output  : {OUTPUT_DIR}")

start_time = time.time()

# Invoke AssemblyFastpTrimNode
node = AssemblyFastpTrimNode()
result = node.run(
    sample_metadata_csv="",
    fastq_dir=str(PHIX_DATA_DIR),
    metadata_csv=str(PHIX_DATA_DIR / "sample_metadata.csv"),
    output_dir=str(OUTPUT_DIR),
    threads=2,
    extra_command="--detect_adapter_for_pe",
)

elapsed = time.time() - start_time

print(f"\n[EXECUTION COMPLETED in {elapsed:.2f}s]")
print(f"• Returned OutDir Tuple : {result}")

# Inspect Generated Files
sample_dir = OUTPUT_DIR / "PHIX174"
r1 = sample_dir / "R1.fastq"
r2 = sample_dir / "R2.fastq"
fastp_json = sample_dir / "fastp.json"
fastp_html = sample_dir / "fastp.html"

assert sample_dir.exists(), f"Sample directory {sample_dir} must exist"
assert r1.exists(), f"Trimmed R1 file {r1} must exist"
assert r2.exists(), f"Trimmed R2 file {r2} must exist"
assert fastp_json.exists(), f"fastp.json report {fastp_json} must exist"
assert fastp_html.exists(), f"fastp.html report {fastp_html} must exist"

# Parse fastp QC stats from JSON
with open(fastp_json, "r") as f:
    qc = json.load(f)

summary = qc.get("summary", {})
before = summary.get("before_filtering", {})
after = summary.get("after_filtering", {})

print("\n" + "-" * 80)
print("📊 [FASTP QC SUMMARY STATISTICS]")
print("-" * 80)
print(f"• Total Reads Before Filtering : {before.get('total_reads', 0):,} reads")
print(f"• Total Reads After Filtering  : {after.get('total_reads', 0):,} reads")
print(f"• Passed Filter Rate          : {(after.get('total_reads', 0) / max(before.get('total_reads', 1), 1) * 100):.2f}%")
print(f"• Q20 Bases Rate (Before ➔ After) : {before.get('q20_rate', 0)*100:.2f}% ➔ {after.get('q20_rate', 0)*100:.2f}%")
print(f"• Q30 Bases Rate (Before ➔ After) : {before.get('q30_rate', 0)*100:.2f}% ➔ {after.get('q30_rate', 0)*100:.2f}%")
print(f"• GC Content (Before ➔ After)     : {before.get('gc_content', 0)*100:.2f}% ➔ {after.get('gc_content', 0)*100:.2f}%")
print(f"• Generated Trimmed R1 Size   : {r1.stat().st_size:,} bytes")
print(f"• Generated Trimmed R2 Size   : {r2.stat().st_size:,} bytes")
print(f"• fastp HTML Report           : {fastp_html} ({fastp_html.stat().st_size:,} bytes)")
print("-" * 80)
print("✅ AssemblyFastpTrimNode Custom Node Test: PASSED (100% Success)")
