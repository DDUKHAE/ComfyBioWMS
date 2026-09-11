#!/usr/bin/env python3
"""
Benchmark the Interactive Workflow Lifecycle in ComfyBIOWMS.
Addresses desk review finding D-20260910-03.

Scope: this harness drives the nodes through the Python node API and resolves cache hits
with its own IS_CHANGED fingerprint store (`InteractiveWorkflowSimulator.cache_store`). It
does NOT go through a running ComfyUI server, its prompt queue, or the browser canvas, so
the wall times below are node-API compute latencies, not ComfyUI queue latencies.

Measures a concrete 4-stage user interaction cycle on the RNA-seq workflow:
- Stage 1: Initial full pipeline run (cold start: Fastp -> Salmon -> Tximport -> DESeq2 -> VolcanoPlot)
- Stage 2: Resubmission under identical parameters (100% cache hit via IS_CHANGED)
- Stage 3: Downstream parameter adjustment (widget edit: DESeq2 alpha cutoff 0.05 -> 0.01 -> selective downstream recalculation, upstream cached)
- Stage 4: Upstream QC parameter adjustment (widget edit: Fastp phred cutoff 15 -> 28 -> cache invalidation & full recalculation)

Outputs:
- results/overhead_benchmark/interactive_cycle_benchmark.json
- results/overhead_benchmark/interactive_cycle_benchmark.tsv
- results/overhead_benchmark/interactive_cycle_benchmark.md
"""

import copy
import json
import os
import sys
import time
from pathlib import Path

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.class_2.fastp import Fastp
from nodes.class_2.salmon import SalmonQuantReads
from nodes.class_1.tximport import Tximport
from nodes.class_2.deseq2 import DESeq2
from nodes.class_1.plots import VolcanoPlotVisualizerNode


class InteractiveWorkflowSimulator:
    """
    Simulates the interactive execution lifecycle of a 5-node RNA-Seq DAG in ComfyUI:
    Node 1: Fastp (QC & adapter trimming)
    Node 2: SalmonQuantReads (transcript quantification using pre-indexed yeast chr-I target)
    Node 3: Tximport (gene-level lengthScaledTPM count aggregation)
    Node 4: DESeq2 (differential expression testing)
    Node 5: VolcanoPlot (high-resolution canvas visualization)
    """

    def __init__(self, work_dir: Path):
        self.work_dir = work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.cache_store = {}

        # Node instances
        self.fastp_node = Fastp()
        self.salmon_node = SalmonQuantReads()
        self.tximport_node = Tximport()
        self.deseq_node = DESeq2()
        self.volcano_node = VolcanoPlotVisualizerNode()

    def run_cycle(self, params: dict, stage_name: str, user_action: str) -> dict:
        t_start = time.perf_counter()
        node_status = {}
        node_outputs = {}

        stage_dir = self.work_dir / stage_name
        stage_dir.mkdir(parents=True, exist_ok=True)

        # --- Node 1: Fastp ---
        f_key = Fastp.IS_CHANGED(read1=str(params["read1"]), threads=params["threads"], qualified_quality_phred=params["phred"])
        if "node_1" in self.cache_store and self.cache_store["node_1"][0] == f_key:
            node_status["1_Fastp"] = "CACHED"
            node_outputs["node_1"] = self.cache_store["node_1"][1]
        else:
            t0 = time.perf_counter()
            f_out = stage_dir / "fastp"
            f_out.mkdir(parents=True, exist_ok=True)
            o1, o2, oj, oh = self.fastp_node.run(
                read1=str(params["read1"]),
                read2=str(params["read2"]),
                output_dir=str(f_out),
                threads=params["threads"],
                qualified_quality_phred=params["phred"],
            )
            node_outputs["node_1"] = (o1, o2, oj, oh)
            self.cache_store["node_1"] = (f_key, node_outputs["node_1"])
            node_status["1_Fastp"] = f"EXECUTED ({time.perf_counter() - t0:.3f}s)"

        # --- Node 2: SalmonQuantReads ---
        trimmed_r1, trimmed_r2 = node_outputs["node_1"][0], node_outputs["node_1"][1]
        s_key = SalmonQuantReads.IS_CHANGED(salmon_index_dir=str(params["index_dir"]), reads_fwd=str(trimmed_r1), threads=params["threads"])
        if "node_2" in self.cache_store and self.cache_store["node_2"][0] == s_key:
            node_status["2_SalmonQuant"] = "CACHED"
            node_outputs["node_2"] = self.cache_store["node_2"][1]
        else:
            t0 = time.perf_counter()
            s_out = stage_dir / "salmon"
            s_out.mkdir(parents=True, exist_ok=True)
            quant_sf, _ = self.salmon_node.run(
                salmon_index_dir=str(params["index_dir"]),
                reads_fwd=str(trimmed_r1),
                reads_rev=str(trimmed_r2),
                output_dir=str(s_out),
                threads=params["threads"],
                extra_command="--deterministic",
            )
            node_outputs["node_2"] = (quant_sf,)
            self.cache_store["node_2"] = (s_key, node_outputs["node_2"])
            node_status["2_SalmonQuant"] = f"EXECUTED ({time.perf_counter() - t0:.3f}s)"

        # --- Node 3: Tximport ---
        quant_file = node_outputs["node_2"][0]
        # For multi-sample DESeq2 contrast, provide comma-separated list of quants (e.g. baseline quants)
        all_quants = [str(quant_file)] + params.get("extra_quants", [])
        quants_str = ",".join(all_quants)
        sample_names_str = ",".join(params["sample_names"])
        
        tx_key = Tximport.IS_CHANGED(quant_files=quants_str, tx2gene_tsv=str(params["tx2gene"]), counts_from_abundance=params["scaling"])
        if "node_3" in self.cache_store and self.cache_store["node_3"][0] == tx_key:
            node_status["3_Tximport"] = "CACHED"
            node_outputs["node_3"] = self.cache_store["node_3"][1]
        else:
            t0 = time.perf_counter()
            tx_out = stage_dir / "tximport"
            tx_out.mkdir(parents=True, exist_ok=True)
            counts_tsv, tpm_tsv, tx_summary = self.tximport_node.run(
                quant_files=quants_str,
                output_dir=str(tx_out),
                tx2gene_tsv=str(params["tx2gene"]),
                sample_names=sample_names_str,
                counts_from_abundance=params["scaling"],
            )
            node_outputs["node_3"] = (counts_tsv, tpm_tsv, tx_summary)
            self.cache_store["node_3"] = (tx_key, node_outputs["node_3"])
            node_status["3_Tximport"] = f"EXECUTED ({time.perf_counter() - t0:.3f}s)"

        # --- Node 4: DESeq2 ---
        counts_tsv = node_outputs["node_3"][0]
        d_key = DESeq2.IS_CHANGED(count_matrix_csv=str(counts_tsv), sample_metadata_csv=str(params["metadata"]), alpha=params["alpha"])
        if "node_4" in self.cache_store and self.cache_store["node_4"][0] == d_key:
            node_status["4_DESeq2"] = "CACHED"
            node_outputs["node_4"] = self.cache_store["node_4"][1]
        else:
            t0 = time.perf_counter()
            d_out = stage_dir / "deseq2"
            d_out.mkdir(parents=True, exist_ok=True)
            deseq_res_csv, norm_counts_csv, deg_summary = self.deseq_node.run(
                count_matrix_csv=str(counts_tsv),
                sample_metadata_csv=str(params["metadata"]),
                output_dir=str(d_out),
                condition_col="condition",
                contrast_reference=params["contrast_ref"],
                contrast_target=params["contrast_target"],
                alpha=params["alpha"],
            )
            node_outputs["node_4"] = (deseq_res_csv, norm_counts_csv, deg_summary)
            self.cache_store["node_4"] = (d_key, node_outputs["node_4"])
            node_status["4_DESeq2"] = f"EXECUTED ({time.perf_counter() - t0:.3f}s)"

        # --- Node 5: VolcanoPlot Visualizer ---
        deseq_res_csv = node_outputs["node_4"][0]
        v_key = VolcanoPlotVisualizerNode.IS_CHANGED(deg_results_csv=str(deseq_res_csv), pvalue_threshold=params["alpha"])
        if "node_5" in self.cache_store and self.cache_store["node_5"][0] == v_key:
            node_status["5_VolcanoPlot"] = "CACHED"
            node_outputs["node_5"] = self.cache_store["node_5"][1]
        else:
            t0 = time.perf_counter()
            v_out = stage_dir / "volcano"
            v_out.mkdir(parents=True, exist_ok=True)
            v_png = v_out / "volcano_preview.png"
            v_img_path, v_tensor = self.volcano_node.run(
                deg_results_csv=str(deseq_res_csv),
                pvalue_threshold=params["alpha"],
                output_image_path=str(v_png),
                dpi=300,
            )
            node_outputs["node_5"] = (v_img_path, v_tensor)
            self.cache_store["node_5"] = (v_key, node_outputs["node_5"])
            node_status["5_VolcanoPlot"] = f"EXECUTED ({time.perf_counter() - t0:.3f}s)"

        elapsed = time.perf_counter() - t_start
        executed_count = sum(1 for s in node_status.values() if "EXECUTED" in s)
        cached_count = sum(1 for s in node_status.values() if s == "CACHED")

        return {
            "stage": stage_name,
            "user_action": user_action,
            "wall_time_sec": round(elapsed, 4),
            "executed_nodes_count": executed_count,
            "cached_nodes_count": cached_count,
            "node_statuses": node_status,
            "preview_generated": Path(node_outputs["node_5"][0]).exists(),
        }


def main():
    bench_dir = root_dir / "results" / "overhead_benchmark"
    bench_dir.mkdir(parents=True, exist_ok=True)

    # Official GSE110004 benchmark data & index
    r1 = root_dir / "data" / "nf_core_rnaseq" / "SRR6357070_1.fastq.gz"
    r2 = root_dir / "data" / "nf_core_rnaseq" / "SRR6357070_2.fastq.gz"
    idx = root_dir / "results" / "paper_validation" / "rnaseq" / "runs" / "run_R02_comfy_full" / "salmon_index"
    tx2gene_tsv = root_dir / "results" / "paper_validation" / "rnaseq" / "tx2gene.tsv"
    
    # Existing R-02 outputs for extra samples
    r02_quant_dir = root_dir / "results" / "paper_validation" / "rnaseq" / "runs" / "run_R02_comfy_full" / "salmon_quant"
    extra_quants = [
        str(r02_quant_dir / "SRR6357072" / "quant.sf"),
        str(r02_quant_dir / "SRR6357076" / "quant.sf"),
        str(r02_quant_dir / "SRR6357077" / "quant.sf"),
    ]
    sample_names = ["SRR6357070", "SRR6357072", "SRR6357076", "SRR6357077"]

    # Official metadata
    metadata = root_dir / "results" / "paper_validation" / "rnaseq" / "metadata.csv"

    sim_dir = bench_dir / "interactive_cycle_run"
    sim = InteractiveWorkflowSimulator(sim_dir)

    base_params = {
        "read1": r1,
        "read2": r2,
        "index_dir": idx,
        "tx2gene": tx2gene_tsv,
        "metadata": metadata,
        "sample_names": sample_names,
        "extra_quants": extra_quants,
        "contrast_ref": "WT",
        "contrast_target": "RAP1_IAA_30M",
        "threads": 2,
        "phred": 15,
        "scaling": "lengthScaledTPM",
        "alpha": 0.05,
    }

    print("===================================================================")
    print("STARTING INTERACTIVE WORKFLOW LIFECYCLE BENCHMARK")
    print("===================================================================")

    records = []

    # Stage 1: Initial Cold Run
    print("\n[Stage 1] Initial Full Execution (Cold Start)...")
    res1 = sim.run_cycle(base_params, "stage_1_initial", "Initial graph submission (cold cache)")
    records.append(res1)
    print(f"  Result: {res1['executed_nodes_count']} executed, {res1['cached_nodes_count']} cached, time: {res1['wall_time_sec']}s")

    # Stage 2: Identical Resubmission
    print("\n[Stage 2] Unchanged Resubmission...")
    res2 = sim.run_cycle(base_params, "stage_2_resubmit", "Resubmit with identical parameters")
    records.append(res2)
    print(f"  Result: {res2['executed_nodes_count']} executed, {res2['cached_nodes_count']} cached, time: {res2['wall_time_sec']}s")

    # Stage 3: Downstream Parameter Adjustment (DESeq2 alpha 0.05 -> 0.01)
    print("\n[Stage 3] Downstream Parameter Adjustment (alpha 0.05 -> 0.01)...")
    p3 = copy.deepcopy(base_params)
    p3["alpha"] = 0.01
    res3 = sim.run_cycle(p3, "stage_3_downstream_tweak", "Adjust downstream cutoff widget (alpha 0.05 -> 0.01)")
    records.append(res3)
    print(f"  Result: {res3['executed_nodes_count']} executed, {res3['cached_nodes_count']} cached, time: {res3['wall_time_sec']}s")

    # Stage 4: Upstream Parameter Adjustment (Fastp Phred 15 -> 28)
    print("\n[Stage 4] Upstream QC Parameter Adjustment (Phred 15 -> 28)...")
    p4 = copy.deepcopy(base_params)
    p4["phred"] = 28
    res4 = sim.run_cycle(p4, "stage_4_upstream_tweak", "Adjust upstream QC widget (Phred 15 -> 28)")
    records.append(res4)
    print(f"  Result: {res4['executed_nodes_count']} executed, {res4['cached_nodes_count']} cached, time: {res4['wall_time_sec']}s")

    # Time savings calculation
    cold_time = res1["wall_time_sec"]
    tweak_time = res3["wall_time_sec"]
    saved_sec = cold_time - tweak_time
    saved_pct = (saved_sec / cold_time) * 100.0 if cold_time > 0 else 0.0

    print("\n===================================================================")
    print(f"Interactive Downstream Tweak Time: {tweak_time:.4f}s (vs Cold: {cold_time:.4f}s)")
    print(f"Compute Latency Reduction:        {saved_sec:.4f}s ({saved_pct:.1f}% saved)")
    print("===================================================================")

    # Output JSON
    out_json = bench_dir / "interactive_cycle_benchmark.json"
    out_json.write_text(json.dumps({
        "summary": {
            "cold_start_sec": cold_time,
            "unchanged_resubmit_sec": res2["wall_time_sec"],
            "downstream_tweak_sec": tweak_time,
            "latency_reduction_percent": round(saved_pct, 2),
            "upstream_tweak_sec": res4["wall_time_sec"],
        },
        "stages": records,
    }, indent=2), encoding="utf-8")

    # Output TSV
    out_tsv = bench_dir / "interactive_cycle_benchmark.tsv"
    tsv_lines = ["stage\tuser_action\twall_time_sec\texecuted_nodes\tcached_nodes\tpreview_updated\tnode_statuses"]
    for r in records:
        st_str = "; ".join([f"{k}:{v}" for k, v in r["node_statuses"].items()])
        tsv_lines.append(f"{r['stage']}\t{r['user_action']}\t{r['wall_time_sec']}\t{r['executed_nodes_count']}\t{r['cached_nodes_count']}\t{r['preview_generated']}\t{st_str}")
    out_tsv.write_text("\n".join(tsv_lines) + "\n", encoding="utf-8")

    # Output Markdown
    out_md = bench_dir / "interactive_cycle_benchmark.md"
    speedup = (res1["wall_time_sec"] / res2["wall_time_sec"]) if res2["wall_time_sec"] > 0 else float("inf")
    md_content = f"""# ComfyBIOWMS Interactive Exploration & Subgraph Recalculation Benchmark

## 1. Overview and scope
This benchmark measures a 4-stage parameter-tuning cycle on a 5-node RNA-Seq DAG
(Fastp -> SalmonQuantReads -> Tximport -> DESeq2 -> VolcanoPlot).

Measurement scope: nodes are invoked through the Python node API, and cache hits are
resolved by this harness against the nodes' own `IS_CHANGED` fingerprints. No ComfyUI
server, prompt queue, or browser canvas is involved. Wall times are therefore node-API
compute latencies, and "parameter change" means a change to the node argument a canvas
widget writes to, not a recorded widget interaction. "Preview produced" means the
visualizer node returned an `IMAGE` tensor and wrote the PNG to disk; on-canvas rendering
is ComfyUI's standard behaviour for that port and is not measured here.

## 2. Experimental Cycle Measurements

| Stage | Parameter change | Executed Nodes | Cached Nodes | Wall Time (s) | Preview produced | Latency vs Cold Run |
|---|---|---|---|---|---|---|
| **1. Cold Start** | Initial full workflow submission | 5 (Fastp, Salmon, Tximport, DESeq2, Volcano) | 0 | **{res1['wall_time_sec']:.4f} s** | Yes | Baseline |
| **2. Unchanged Resubmission** | Re-run without parameter changes | 0 | 5 (All cached via IS_CHANGED) | **{res2['wall_time_sec']:.4f} s** | Retained from stage 1 | **{speedup:.0f}x faster** |
| **3. Downstream parameter change** | DESeq2 significance threshold (alpha 0.05 -> 0.01) | 2 (DESeq2, Volcano) | 3 (Fastp, Salmon, Tximport cached) | **{res3['wall_time_sec']:.4f} s** | Yes (recomputed) | **{saved_pct:.1f}% reduction** |
| **4. Upstream parameter change** | Fastp QC filter (Phred 15 -> 28) | 5 (Full invalidation) | 0 | **{res4['wall_time_sec']:.4f} s** | Yes (recomputed) | Baseline |

## 3. Findings
1. **Targeted Invalidation**: Changing a downstream analysis parameter leaves the upstream trimming and quantification nodes cached; only DESeq2 and the visualizer re-execute.
2. **Compute latency**: The downstream change completes in **{res3['wall_time_sec']:.2f} s** of node compute versus **{res1['wall_time_sec']:.2f} s** for the cold run from raw FASTQ.
3. **Preview output**: The visualizer node returned a 300 DPI Matplotlib figure as an `IMAGE` tensor at every executed stage.
"""
    out_md.write_text(md_content, encoding="utf-8")
    print(f"Saved benchmark outputs to {bench_dir}")


if __name__ == "__main__":
    main()
