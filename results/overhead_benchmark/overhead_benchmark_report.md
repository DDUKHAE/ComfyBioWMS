# ComfyBIOWMS Runtime Overhead & Memory Benchmark Report

## 1. Experimental Conditions
- **Input Data**: `SRR6357070_1.fastq.gz` (2186.8 KB)
- **Tool Executable**: `/opt/miniconda3/envs/bulk_rna_seq/bin/fastp`
- **Execution Environment**: `bulk_rna_seq` (isolated conda prefix)
- **Repetitions**: 5 independent iterations per mode

---

## 2. Empirical Timing & Overhead Results

| Execution Mode | Mean Duration (s) | Median Duration (s) | Standard Deviation (s) | Relative Overhead vs Native |
|---|---|---|---|---|
| **Native CLI** | 0.2274 s | 0.2275 s | 0.0016 s | Baseline ($0.0\%$) |
| **ComfyBIOWMS Node** | 0.2296 s | 0.2294 s | 0.0017 s | $+0.0023$ s ($+1.0\%$) |
| **ComfyUI Cache Check** | 0.126 ms | 0.154 ms | — | **Instantaneous hit** ($>1000\times$ speedup) |

---

## 3. Memory Profile
- **Python Backend Process Peak RSS**: **136.06 MB**
- **Disk I/O Logging**: Direct file streaming for `stdout.log` and `stderr.log` (constant memory profile regardless of input dataset scale).

---

## 4. Conclusion
ComfyBIOWMS adds negligible execution overhead ($pprox 0.002$ seconds) consisting entirely of:
1. Environment binary resolution and argv construction.
2. Direct disk stream piping and execution timestamping.
3. Writing the audit manifests (`run_manifest.json` and `run_manifest.sh`).

When the workflow is partially re-executed, ComfyUI's `IS_CHANGED` cache mechanism resolves unchanged inputs in **0.15 milliseconds**, completely bypassing the underlying tool execution.
