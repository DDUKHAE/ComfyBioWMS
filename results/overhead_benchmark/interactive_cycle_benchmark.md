# ComfyBIOWMS Interactive Exploration & Subgraph Recalculation Benchmark

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
| **1. Cold Start** | Initial full workflow submission | 5 (Fastp, Salmon, Tximport, DESeq2, Volcano) | 0 | **4.4519 s** | Yes | Baseline |
| **2. Unchanged Resubmission** | Re-run without parameter changes | 0 | 5 (All cached via IS_CHANGED) | **0.0015 s** | Retained from stage 1 | **2968x faster** |
| **3. Downstream parameter change** | DESeq2 significance threshold (alpha 0.05 -> 0.01) | 2 (DESeq2, Volcano) | 3 (Fastp, Salmon, Tximport cached) | **3.9201 s** | Yes (recomputed) | **11.9% reduction** |
| **4. Upstream parameter change** | Fastp QC filter (Phred 15 -> 28) | 5 (Full invalidation) | 0 | **4.3508 s** | Yes (recomputed) | Baseline |

## 3. Findings
1. **Targeted Invalidation**: Changing a downstream analysis parameter leaves the upstream trimming and quantification nodes cached; only DESeq2 and the visualizer re-execute.
2. **Compute latency**: The downstream change completes in **3.92 s** of node compute versus **4.45 s** for the cold run from raw FASTQ.
3. **Preview output**: The visualizer node returned a 300 DPI Matplotlib figure as an `IMAGE` tensor at every executed stage.
