import time
import subprocess
import shutil
import json
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[2]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from bioflow.runtime.command_runner import CondaCommandRunner

def measure_command(cmd_argv: list[str], cwd: Path) -> dict:
    time_cmd = ["/usr/bin/time", "-l"] + cmd_argv
    start_t = time.perf_counter()
    res = subprocess.run(time_cmd, cwd=cwd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start_t
    
    peak_rss = 0
    for line in res.stderr.splitlines():
        if "maximum resident set size" in line:
            parts = line.strip().split()
            if parts:
                peak_rss = int(parts[0]) // 1024
                break
    return {
        "elapsed_sec": round(elapsed, 4),
        "peak_rss_kb": peak_rss,
        "returncode": res.returncode
    }

def main():
    bench_dir = root_dir / "results" / "overhead_benchmark"
    bench_dir.mkdir(parents=True, exist_ok=True)
    
    r1_test = root_dir / "data" / "assembly_phix174" / "phix174_R1.fastq.gz"
    r2_test = root_dir / "data" / "assembly_phix174" / "phix174_R2.fastq.gz"
    
    if not r1_test.exists() or not r2_test.exists():
        print(f"Test FASTQ not found at {r1_test}")
        return
        
    out_native = bench_dir / "native_out"
    out_comfy = bench_dir / "comfy_out"
    out_native.mkdir(parents=True, exist_ok=True)
    out_comfy.mkdir(parents=True, exist_ok=True)
    
    native_argv = [
        "conda", "run", "-n", "genome_assembly", "fastp",
        "-i", str(r1_test), "-I", str(r2_test),
        "-o", str(out_native / "R1.trimmed.fastq.gz"),
        "-O", str(out_native / "R2.trimmed.fastq.gz"),
        "-j", str(out_native / "fastp.json"),
        "-h", str(out_native / "fastp.html"),
        "-w", "2"
    ]
    
    print("Measuring Native CLI fastp (3 repeats)...")
    native_times = []
    for _ in range(3):
        m = measure_command(native_argv, bench_dir)
        native_times.append(m)
        
    avg_native_time = sum(m["elapsed_sec"] for m in native_times) / len(native_times)
    avg_native_rss = sum(m["peak_rss_kb"] for m in native_times) / len(native_times)
    
    print("Measuring ComfyBIOWMS CondaCommandRunner fastp (3 repeats)...")
    comfy_times = []
    for _ in range(3):
        runner = CondaCommandRunner()
        comfy_argv = [
            "conda", "run", "-n", "genome_assembly", "fastp",
            "-i", str(r1_test), "-I", str(r2_test),
            "-o", str(out_comfy / "R1.trimmed.fastq.gz"),
            "-O", str(out_comfy / "R2.trimmed.fastq.gz"),
            "-j", str(out_comfy / "fastp.json"),
            "-h", str(out_comfy / "fastp.html"),
            "-w", "2"
        ]
        start_t = time.perf_counter()
        runner.run(comfy_argv, bench_dir)
        elapsed = time.perf_counter() - start_t
        comfy_times.append({"elapsed_sec": round(elapsed, 4)})
        
    avg_comfy_time = sum(m["elapsed_sec"] for m in comfy_times) / len(comfy_times)
    overhead_sec = avg_comfy_time - avg_native_time
    overhead_pct = (overhead_sec / avg_native_time) * 100 if avg_native_time > 0 else 0.0
    
    # Process memory
    import resource
    peak_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    if sys.platform == "darwin":
        # On macOS, ru_maxrss is in bytes
        peak_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    
    summary = {
        "native_cli_avg_sec": round(avg_native_time, 3),
        "native_cli_peak_rss_kb": round(avg_native_rss, 1),
        "comfybio_avg_sec": round(avg_comfy_time, 3),
        "overhead_sec": round(overhead_sec, 3),
        "overhead_percent": round(overhead_pct, 2),
        "backend_memory_rss_mb": round(peak_rss_mb, 2)
    }
    
    report_path = bench_dir / "overhead_report.json"
    report_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\nBenchmark Summary:")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
