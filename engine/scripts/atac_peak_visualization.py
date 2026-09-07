#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def parse_narrowpeak(peak_file: Path):
    peaks = []
    with peak_file.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                signal = float(parts[6])
                p_val = float(parts[7])
                q_val = float(parts[8])
                summit_offset = int(parts[9]) if len(parts) > 9 else (end - start) // 2
                peaks.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "name": name,
                    "signal": signal,
                    "pval": p_val,
                    "qval": q_val,
                    "summit": start + summit_offset,
                })
    return peaks

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot publication-quality Genome Browser-style ATAC-seq peak tracks matching Buenrostro et al. 2013 Fig 2.")
    parser.add_argument("--peaks-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    files = list(args.peaks_dir.glob("*/*_peaks.narrowPeak")) + list(args.peaks_dir.glob("*_peaks.narrowPeak"))
    if not files:
        raise SystemExit(f"No narrowPeak files found under {args.peaks_dir}")

    peak_file = sorted(files)[0]
    peaks = parse_narrowpeak(peak_file)
    sample_name = peak_file.name.removesuffix("_peaks.narrowPeak")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Styling for Nature Methods publication aesthetic
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
    plt.rcParams["axes.edgecolor"] = "#334155"
    plt.rcParams["axes.linewidth"] = 1.0

    fig = plt.figure(figsize=(12, 6), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.8, 1.2], hspace=0.35)

    # Panel 1: Whole Chromosome 22 ATAC-seq Peak Signal Track (Genome Browser style)
    ax1 = fig.add_subplot(gs[0])
    chr22_length_mb = 51.3  # Human chr22 length ~51.3 Mb

    # Plot genomic coverage baseline
    x_grid = np.linspace(15, 51.3, 1000)
    y_baseline = np.zeros_like(x_grid)
    ax1.plot(x_grid, y_baseline, color="#94a3b8", linewidth=1.5)

    if peaks:
        for p in peaks:
            center_mb = p["summit"] / 1e6
            width_mb = max((p["end"] - p["start"]) / 1e6, 0.05)
            sig = p["signal"]
            
            # Local window for peak bell shape
            x_local = np.linspace(center_mb - 0.25, center_mb + 0.25, 100)
            y_local = sig * np.exp(-0.5 * ((x_local - center_mb) / (width_mb * 0.8)) ** 2)
            ax1.fill_between(x_local, y_local, color="#2563eb", alpha=0.75)
            ax1.plot(x_local, y_local, color="#1d4ed8", linewidth=1.2)
            
            # Top summit marker
            if sig > 15.0:
                ax1.scatter([center_mb], [sig + 1.5], color="#dc2626", s=25, zorder=5)
                ax1.annotate(f"{p['name']}\n(Signal: {sig:.1f})", xy=(center_mb, sig + 2.0),
                             xytext=(center_mb, sig + 7.0),
                             fontsize=7.5, fontweight="bold", color="#1e293b", ha="center",
                             arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.0))

        max_sig = max([p["signal"] for p in peaks])
        ax1.set_ylim(0, max_sig * 1.35 if max_sig > 0 else 10)
    else:
        ax1.text(33.0, 5.0, "No significant open chromatin peaks detected", ha="center", va="center", color="#64748b", fontweight="bold")
        ax1.set_ylim(0, 10)

    ax1.set_title(f"Epigenomics Open Chromatin Profile: {sample_name} (Buenrostro et al. 2013 Genome Browser Format)",
                  fontsize=12, fontweight="bold", pad=12, color="#0f172a")
    ax1.set_xlabel("Genomic Position on Chromosome 22 (Mb)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("ATAC-Seq Signal Value\n(MACS3 Pileup Intensity)", fontsize=10, fontweight="bold")
    ax1.set_xlim(15.0, 51.5)
    ax1.grid(True, linestyle="--", alpha=0.3, color="#cbd5e1")

    # Panel 2: Zoom-in on the High-Confidence Peak Region
    ax2 = fig.add_subplot(gs[1])
    if peaks:
        top_peak = max(peaks, key=lambda p: p["signal"])
        zoom_center_mb = top_peak["summit"] / 1e6
        zoom_x = np.linspace(zoom_center_mb - 0.05, zoom_center_mb + 0.05, 300)
        zoom_y = top_peak["signal"] * np.exp(-0.5 * ((zoom_x - zoom_center_mb) / 0.008) ** 2)

        ax2.fill_between(zoom_x, zoom_y, color="#3b82f6", alpha=0.6)
        ax2.plot(zoom_x, zoom_y, color="#1d4ed8", linewidth=2.0)
        ax2.axvline(zoom_center_mb, color="#dc2626", linestyle=":", linewidth=1.5, label="Tn5 Transposase Cleavage Summit")

        start_mb = top_peak["start"] / 1e6
        end_mb = top_peak["end"] / 1e6
        rect = patches.Rectangle((start_mb, -3), end_mb - start_mb, 2, linewidth=1, edgecolor="#1e40af", facecolor="#10b981", label="Called Peak Interval (MAPQ >= 30)")
        ax2.add_patch(rect)

        ax2.set_title(f"Local High-Confidence Open Chromatin Locus: Chr22: {top_peak['start']:,} – {top_peak['end']:,} bp (Signal = {top_peak['signal']:.1f}, q-val = 10^{{-{top_peak['qval']:.1f}}})",
                      fontsize=10, fontweight="bold", pad=8, color="#0f172a")
        ax2.set_xlim(zoom_center_mb - 0.04, zoom_center_mb + 0.04)
        ax2.set_ylim(-4, top_peak["signal"] * 1.25 if top_peak["signal"] > 0 else 10)
        ax2.legend(loc="upper right", fontsize=8, frameon=True)
    else:
        ax2.text(0.5, 0.5, "No peak locus available for zoom-in view", ha="center", va="center", color="#64748b", fontweight="bold", transform=ax2.transAxes)
        ax2.set_title("Local High-Confidence Open Chromatin Locus (No Peaks)", fontsize=10, fontweight="bold", pad=8, color="#0f172a")

    ax2.set_xlabel("Genomic Coordinates (Mb)", fontsize=9, fontweight="bold")
    ax2.set_ylabel("Signal", fontsize=9, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.3, color="#cbd5e1")

    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print(f"✅ Successfully generated Nature Methods-style ATAC-seq track plot: {args.output}")

if __name__ == "__main__":
    main()
