"""
Script to generate Figure 3: Interactive Canvas Usability, Provenance Tracking, and One-Drop PNG Replication.
Saves to paper_platform/figures/figure_3_platform_features.png at 300 DPI.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec

# Publication aesthetics
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig = plt.figure(figsize=(18, 10), dpi=300)
gs = GridSpec(1, 3, figure=fig, wspace=0.18)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

# -------------------------------------------------------------------------
# Subpanel (a): Interactive No-Code Parameter Tuning on Canvas
# -------------------------------------------------------------------------
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.axis("off")
ax1.set_title("(a) Interactive No-Code Parameter Tuning", fontsize=12, fontweight="bold", pad=12, loc="left", color="#111111")

# Node Mockup
node_box = patches.FancyBboxPatch((0.5, 2.5), 9.0, 7.0, boxstyle="round,pad=0.2", fc="#2d3748", ec="#4a5568", lw=1.5)
ax1.add_patch(node_box)
# Header
header_box = patches.FancyBboxPatch((0.5, 8.2), 9.0, 1.3, boxstyle="round,pad=0.1", fc="#3182ce", ec="none")
ax1.add_patch(header_box)
ax1.text(5.0, 8.85, "VolcanoPlotVisualizerNode", ha="center", va="center", color="white", fontsize=11, fontweight="bold")

# Sliders & Widgets
widgets = [
    ("Fold Change Cutoff (|log2FC|)", "1.5", "[ 0.5 ────●──── 5.0 ]", "#63b3ed"),
    ("Significance Cutoff (padj)", "0.01", "[ 1e-5 ──●────── 0.05 ]", "#63b3ed"),
    ("Top Genes to Label", "15", "[ 0 ──────●─── 50 ]", "#63b3ed"),
    ("Publication Style Preset", "Nature / Cell", "[ ▼ Nature (300 DPI) ]", "#9ae6b4"),
    ("Native Output Socket", "IMAGE Tensor", "● preview_image (Live)", "#fbd38d")
]

y_w = 7.4
for name, val, slider, col in widgets:
    ax1.text(0.9, y_w, f"{name}:", fontsize=8.5, fontweight="bold", color="#e2e8f0")
    ax1.text(9.0, y_w, val, ha="right", fontsize=8.5, fontweight="bold", color=col)
    ax1.text(5.0, y_w - 0.45, slider, ha="center", fontsize=8.0, color="#a0aec0")
    y_w -= 1.1

# Bottom reaction banner
react_box = patches.FancyBboxPatch((0.5, 0.4), 9.0, 1.6, boxstyle="round,pad=0.15", fc="#ebf8ff", ec="#bee3f8", lw=1.0)
ax1.add_patch(react_box)
ax1.text(5.0, 1.2, "⚡ Live Canvas Interaction Loop", ha="center", fontsize=9.5, fontweight="bold", color="#2b6cb0")
ax1.text(5.0, 0.7, "Slider change triggers instant sub-second\nplot refresh without full pipeline re-run", 
         ha="center", fontsize=7.8, color="#2d3748")

# -------------------------------------------------------------------------
# Subpanel (b): Automated Provenance Auditing (artifacts.sidecar.json)
# -------------------------------------------------------------------------
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis("off")
ax2.set_title("(b) Automated Provenance Auditing (Sidecar)", fontsize=12, fontweight="bold", pad=12, loc="left", color="#111111")

sidecar_box = patches.FancyBboxPatch((0.5, 0.4), 9.0, 9.1, boxstyle="round,pad=0.2", fc="#1a202c", ec="#4a5568", lw=1.5)
ax2.add_patch(sidecar_box)

header_sidecar = patches.FancyBboxPatch((0.5, 8.2), 9.0, 1.3, boxstyle="round,pad=0.1", fc="#805ad5", ec="none")
ax2.add_patch(header_sidecar)
ax2.text(5.0, 8.85, "artifacts.sidecar.json (Audit Log)", ha="center", va="center", color="white", fontsize=11, fontweight="bold")

json_text = """{
  "execution_id": "exec_2026_seqc_004",
  "pipeline": "bulk_rna_seq_salmon_ref",
  "tool_provenance": {
    "tool_id": "salmon_quant_01",
    "version": "Salmon v1.10.2",
    "conda_env": "comfybio_rnaseq_env",
    "cli_command": "salmon quant -i index -l A -1 r1.fq -2 r2.fq --validateMappings -o quant_out",
    "execution_time_sec": 42.15,
    "max_memory_mb": 2450.8
  },
  "input_hashes": {
    "sample_1_R1.fastq.gz": "sha256:7f83b1657ff1...",
    "sample_1_R2.fastq.gz": "sha256:8a1b2c3d4e5f..."
  },
  "output_artifacts": {
    "quant.sf": "sha256:4b9a7c8e0d1f...",
    "volcano_figure.png": "300 DPI, 7.0x5.5 in"
  },
  "status": "SUCCESS_AUDITED"
}"""

ax2.text(0.8, 4.3, json_text, fontsize=7.2, family="monospace", color="#68d391", va="center")

# -------------------------------------------------------------------------
# Subpanel (c): One-Drop PNG Workflow Serialization & Replication
# -------------------------------------------------------------------------
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis("off")
ax3.set_title("(c) One-Drop PNG Workflow Replication", fontsize=12, fontweight="bold", pad=12, loc="left", color="#111111")

# Top: PNG figure file representation
fig_box = patches.FancyBboxPatch((1.5, 5.8), 7.0, 3.7, boxstyle="round,pad=0.15", fc="#fffaf0", ec="#dd6b20", lw=1.5)
ax3.add_patch(fig_box)
ax3.text(5.0, 8.8, "volcano_plot_nature.png", ha="center", fontsize=9.5, fontweight="bold", color="#c05621")
ax3.text(5.0, 7.7, "[Rendered 300 DPI Figure Content]\n+ Embedded PNG 'tEXt' Chunk:\n{ full_dag_nodes, edges, parameters }", 
         ha="center", va="center", fontsize=7.8, color="#7b341e")

# Drag & Drop Action Arrow
ax3.annotate("Drag & Drop PNG onto Canvas", xy=(5.0, 3.8), xytext=(5.0, 5.4),
             arrowprops=dict(arrowstyle="->", color="#dd6b20", lw=2.5, mutation_scale=18),
             ha="center", fontsize=9.0, fontweight="bold", color="#c05621")

# Bottom: Reconstructed Canvas Graph representation
canvas_box = patches.FancyBboxPatch((1.0, 0.4), 8.0, 3.2, boxstyle="round,pad=0.15", fc="#f7fafc", ec="#4a5568", lw=1.2)
ax3.add_patch(canvas_box)
ax3.text(5.0, 3.0, "ComfyBioWMS Interactive Canvas", ha="center", fontsize=9.5, fontweight="bold", color="#2d3748")

# Mini nodes inside canvas
mini_n1 = patches.FancyBboxPatch((1.5, 0.8), 1.8, 1.4, boxstyle="round,pad=0.08", fc="#3182ce", ec="none")
mini_n2 = patches.FancyBboxPatch((4.1, 0.8), 1.8, 1.4, boxstyle="round,pad=0.08", fc="#38a169", ec="none")
mini_n3 = patches.FancyBboxPatch((6.7, 0.8), 1.8, 1.4, boxstyle="round,pad=0.08", fc="#dd6b20", ec="none")
ax3.add_patch(mini_n1)
ax3.add_patch(mini_n2)
ax3.add_patch(mini_n3)
ax3.text(2.4, 1.5, "Salmon\nQuant", ha="center", va="center", color="white", fontsize=7.0, fontweight="bold")
ax3.text(5.0, 1.5, "DESeq2\nAnalysis", ha="center", va="center", color="white", fontsize=7.0, fontweight="bold")
ax3.text(7.6, 1.5, "Volcano\nPlot", ha="center", va="center", color="white", fontsize=7.0, fontweight="bold")

# Connecting lines
ax3.plot([3.3, 4.1], [1.5, 1.5], color="#4a5568", lw=1.5)
ax3.plot([5.9, 6.7], [1.5, 1.5], color="#4a5568", lw=1.5)

ax3.text(5.0, 0.55, "✓ 100% Graph & Parameter Reconstruction", ha="center", fontsize=7.8, fontweight="bold", color="#276749")

# Save
output_path = "paper_platform/figures/figure_3_platform_features.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Successfully generated Figure 3 at {output_path}")
