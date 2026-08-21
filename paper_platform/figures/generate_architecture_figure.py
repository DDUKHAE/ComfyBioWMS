"""
Script to generate Figure 1: 3-Tier System Architecture & End-to-End Execution Dataflow of ComfyBioflow.
Saves to paper_platform/figures/figure_1_system_architecture.png at 300 DPI.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Publication aesthetics
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig = plt.figure(figsize=(16, 11), dpi=300)
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(0, 16)
ax.set_ylim(0, 11)
ax.axis("off")

# Title
ax.text(8.0, 10.5, "ComfyBioflow: Three-Tier System Architecture & Dataflow Wiring", 
        ha="center", va="center", fontsize=15, fontweight="bold", color="#111111")

# =========================================================================
# LAYER 1: Visual Interface & Interaction Layer (Canvas)
# =========================================================================
l1_box = patches.FancyBboxPatch((0.8, 7.1), 14.4, 2.8, boxstyle="round,pad=0.25", 
                                fc="#f0f4f8", ec="#1d3557", lw=2.0)
ax.add_patch(l1_box)

# Header
h1_box = patches.FancyBboxPatch((1.0, 9.1), 14.0, 0.65, boxstyle="round,pad=0.1", fc="#1d3557", ec="none")
ax.add_patch(h1_box)
ax.text(8.0, 9.42, "Layer 1: Visual Interface & Interactive Presentation Layer (Web Canvas UI)", 
        ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")

# Sub-components of Layer 1
c1_1 = patches.FancyBboxPatch((1.2, 7.3), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#a8dadc", lw=1.2)
c1_2 = patches.FancyBboxPatch((4.7, 7.3), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#a8dadc", lw=1.2)
c1_3 = patches.FancyBboxPatch((8.2, 7.3), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#a8dadc", lw=1.2)
c1_4 = patches.FancyBboxPatch((11.7, 7.3), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#a8dadc", lw=1.2)
ax.add_patch(c1_1)
ax.add_patch(c1_2)
ax.add_patch(c1_3)
ax.add_patch(c1_4)

ax.text(2.8, 8.5, "Interactive DAG Canvas", ha="center", fontsize=9.5, fontweight="bold", color="#1d3557")
ax.text(2.8, 7.8, "• Drag-and-drop nodes\n• Mouse wire port linking\n• Zoom & Pan controls", ha="center", fontsize=8.0, color="#444444")

ax.text(6.3, 8.5, "Live Parameter Controls", ha="center", fontsize=9.5, fontweight="bold", color="#1d3557")
ax.text(6.3, 7.8, "• Numeric slider widgets\n• Dynamic cutoffs (padj, FC)\n• Color palette selectors", ha="center", fontsize=8.0, color="#444444")

ax.text(9.8, 8.5, "Native IMAGE Preview", ha="center", fontsize=9.5, fontweight="bold", color="#1d3557")
ax.text(9.8, 7.8, "• PyTorch tensor rendering\n• 300+ DPI live plot display\n• Nature/Cell style presets", ha="center", fontsize=8.0, color="#444444")

ax.text(13.3, 8.5, "One-Drop Replication", ha="center", fontsize=9.5, fontweight="bold", color="#1d3557")
ax.text(13.3, 7.8, "• PNG 'tEXt' metadata embedding\n• Drag image to restore graph\n• 100% parameter fidelity", ha="center", fontsize=8.0, color="#444444")

# Flow Arrow Layer 1 -> Layer 2
ax.annotate("", xy=(8.0, 6.7), xytext=(8.0, 7.1),
            arrowprops=dict(arrowstyle="->", color="#1d3557", lw=2.5, mutation_scale=18))
ax.text(8.3, 6.9, "Graph Topology & Slider Parameters (JSON)", va="center", fontsize=8.5, fontweight="bold", color="#1d3557")

# =========================================================================
# LAYER 2: Graph Execution & Contract Engine Layer (Backend)
# =========================================================================
l2_box = patches.FancyBboxPatch((0.8, 3.7), 14.4, 2.8, boxstyle="round,pad=0.25", 
                                fc="#f4fbf7", ec="#2a9d8f", lw=2.0)
ax.add_patch(l2_box)

# Header
h2_box = patches.FancyBboxPatch((1.0, 5.7), 14.0, 0.65, boxstyle="round,pad=0.1", fc="#2a9d8f", ec="none")
ax.add_patch(h2_box)
ax.text(8.0, 6.02, "Layer 2: Graph Execution & Contract Engine Layer (Platform Backend)", 
        ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")

# Sub-components of Layer 2
c2_1 = patches.FancyBboxPatch((1.2, 3.9), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#b7e4c7", lw=1.2)
c2_2 = patches.FancyBboxPatch((4.7, 3.9), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#b7e4c7", lw=1.2)
c2_3 = patches.FancyBboxPatch((8.2, 3.9), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#b7e4c7", lw=1.2)
c2_4 = patches.FancyBboxPatch((11.7, 3.9), 3.2, 1.6, boxstyle="round,pad=0.15", fc="white", ec="#b7e4c7", lw=1.2)
ax.add_patch(c2_1)
ax.add_patch(c2_2)
ax.add_patch(c2_3)
ax.add_patch(c2_4)

ax.text(2.8, 5.1, "157 Standard Nodes", ha="center", fontsize=9.5, fontweight="bold", color="#2a9d8f")
ax.text(2.8, 4.4, "• 8 Biological domains\n• 24 Gold-standard pipelines\n• Community standard curated", ha="center", fontsize=8.0, color="#444444")

ax.text(6.3, 5.1, "Socket Type Contracts", ha="center", fontsize=9.5, fontweight="bold", color="#2a9d8f")
ax.text(6.3, 4.4, "• FASTQ / BAM / VCF typed\n• Count Matrix / BED sockets\n• Pre-execution crash block", ha="center", fontsize=8.0, color="#444444")

ax.text(9.8, 5.1, "Topological Smart Cache", ha="center", fontsize=9.5, fontweight="bold", color="#2a9d8f")
ax.text(9.8, 4.4, "• Upstream file hash retention\n• Partial downstream execution\n• Sub-second plot refresh", ha="center", fontsize=8.0, color="#444444")

ax.text(13.3, 5.1, "20 Tier-1 Visualizers", ha="center", fontsize=9.5, fontweight="bold", color="#2a9d8f")
ax.text(13.3, 4.4, "• Volcano, Clustered Heatmap\n• UMAP, Coverage Track\n• Manhattan, Taxonomic Bar", ha="center", fontsize=8.0, color="#444444")

# Flow Arrow Layer 2 -> Layer 3
ax.annotate("", xy=(8.0, 3.3), xytext=(8.0, 3.7),
            arrowprops=dict(arrowstyle="->", color="#2a9d8f", lw=2.5, mutation_scale=18))
ax.text(8.3, 3.5, "Validated Execution Jobs & Domain Environment IDs", va="center", fontsize=8.5, fontweight="bold", color="#2a9d8f")

# =========================================================================
# LAYER 3: Isolated Runtime Execution & Provenance Layer (System)
# =========================================================================
l3_box = patches.FancyBboxPatch((0.8, 0.4), 14.4, 2.7, boxstyle="round,pad=0.25", 
                                fc="#fbf5f4", ec="#e76f51", lw=2.0)
ax.add_patch(l3_box)

# Header
h3_box = patches.FancyBboxPatch((1.0, 2.3), 14.0, 0.65, boxstyle="round,pad=0.1", fc="#e76f51", ec="none")
ax.add_patch(h3_box)
ax.text(8.0, 2.62, "Layer 3: Isolated Runtime Execution & Provenance Layer (System Infrastructure)", 
        ha="center", va="center", color="white", fontsize=11.5, fontweight="bold")

# Sub-components of Layer 3
c3_1 = patches.FancyBboxPatch((1.2, 0.6), 4.3, 1.5, boxstyle="round,pad=0.15", fc="white", ec="#f4a261", lw=1.2)
c3_2 = patches.FancyBboxPatch((5.8, 0.6), 4.3, 1.5, boxstyle="round,pad=0.15", fc="white", ec="#f4a261", lw=1.2)
c3_3 = patches.FancyBboxPatch((10.4, 0.6), 4.5, 1.5, boxstyle="round,pad=0.15", fc="white", ec="#f4a261", lw=1.2)
ax.add_patch(c3_1)
ax.add_patch(c3_2)
ax.add_patch(c3_3)

ax.text(3.35, 1.75, "Domain-Isolated Conda Environments", ha="center", fontsize=9.5, fontweight="bold", color="#e76f51")
ax.text(3.35, 1.15, "• env_rnaseq (fastp, Salmon, DESeq2)\n• env_variant (BWA-MEM2, BCFtools)\n• env_metag, env_assembly, env_cadd", ha="center", fontsize=8.0, color="#444444")

ax.text(7.95, 1.75, "Subprocess Execution Dispatcher", ha="center", fontsize=9.5, fontweight="bold", color="#e76f51")
ax.text(7.95, 1.15, "• Deterministic binary invocation\n• Background non-blocking execution\n• Process monitoring & exit code checking", ha="center", fontsize=8.0, color="#444444")

ax.text(12.65, 1.75, "Sidecar Provenance Logger", ha="center", fontsize=9.5, fontweight="bold", color="#e76f51")
ax.text(12.65, 1.15, "• artifacts.sidecar.json generation\n• Exact CLI command + tool versions\n• SHA-256 input/output hashes & time/memory", ha="center", fontsize=8.0, color="#444444")

# Return Flow Arrow from Layer 3 back to Layer 1 (IMAGE tensors and Sidecar logs)
ax.annotate("", xy=(0.5, 8.5), xytext=(0.5, 1.8),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-0.15", color="#8338ec", lw=2.2, mutation_scale=16))
ax.text(0.1, 5.2, "Native IMAGE Tensors &\nAudited Sidecar Artifacts", 
        ha="center", va="center", rotation=90, fontsize=8.5, fontweight="bold", color="#8338ec")

# Save figure
output_path = "paper_platform/figures/figure_1_system_architecture.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Successfully generated Figure 1 at {output_path}")
