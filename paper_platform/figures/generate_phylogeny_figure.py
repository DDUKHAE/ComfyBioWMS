"""
Script to generate publication-grade Figure: Node Phylogeny Tree of ComfyBioWMS Ecosystem.
Saves to paper_platform/figures/figure_node_phylogeny_tree.png at 300 DPI.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import numpy as np

# Set high-quality styling
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

fig = plt.figure(figsize=(18, 12), dpi=300)
gs = GridSpec(2, 2, width_ratios=[1.25, 0.75], height_ratios=[1, 1], figure=fig, wspace=0.15, hspace=0.2)

ax_tree = fig.add_subplot(gs[:, 0])
ax_donut = fig.add_subplot(gs[0, 1])
ax_socket = fig.add_subplot(gs[1, 1])

# ---------------------------------------------------------
# Panel A: Phylogenetic Tree (Cladogram) of 157 Nodes
# ---------------------------------------------------------
ax_tree.set_xlim(-1, 10.5)
ax_tree.set_ylim(-0.5, 9.5)
ax_tree.axis("off")
ax_tree.set_title("(a) Taxonomic Phylogeny of ComfyBioWMS Node Ecosystem (157 Nodes)", 
                   fontsize=14, fontweight="bold", pad=15, loc="left", color="#111111")

kingdoms = [
    {
        "name": "1. Transcriptomics & Spatial",
        "count": 25,
        "color": "#1d3557",
        "phyla": ["Bulk Core (13)", "Single-Cell (12)", "Spatial (4)", "Isoform (3)"],
        "nodes": "Fastp, Salmon, DESeq2, Scanpy, scVI, Squidpy, StringTie2..."
    },
    {
        "name": "2. Genomics & Variants",
        "count": 16,
        "color": "#2a9d8f",
        "phyla": ["Short-Read (5)", "Variant Calling (4)", "Long-Read/SV (4)", "Intervals (3)"],
        "nodes": "BWA-MEM2, MarkDup, BCFtools, DeepVariant, Minimap2, Sniffles2..."
    },
    {
        "name": "3. Epigenomics & 3D Chromatin",
        "count": 17,
        "color": "#e76f51",
        "phyla": ["ATAC-Seq (7)", "ChIP/Motif (4)", "3D Hi-C (3)", "Methyl/CRISPR (3)"],
        "nodes": "MACS3, PyGenrich, HOMER, MEME, Cooler, Cooltools, MAGeCK..."
    },
    {
        "name": "4. Metagenomics & Virology",
        "count": 14,
        "color": "#e9c46a",
        "phyla": ["Taxonomic (4)", "Amplicon (4)", "Functional (1)", "Virome/Phage (3)"],
        "nodes": "Kraken2, Bracken, MetaPhlAn4, DADA2, FastUniFrac, HUMAnN3, VirSorter2..."
    },
    {
        "name": "5. De Novo Assembly & Annotation",
        "count": 8,
        "color": "#457b9d",
        "phyla": ["Isolate Assembly (2)", "Long-Read Assembly (2)", "QC (2)", "Annotation (2)"],
        "nodes": "SPAdes, Flye, Hifiasm, Racon, QUAST, Bakta, Prokka..."
    },
    {
        "name": "6. Proteomics & Metabolomics",
        "count": 16,
        "color": "#8338ec",
        "phyla": ["LC-MS/MS DDA/DIA (6)", "MS Features (4)", "Metabolomics (5)", "Neoantigen (1)"],
        "nodes": "MaxQuant, DIA-NN, MSFragger, PyOpenMS, SIRIUS, MetaboAnalystR..."
    },
    {
        "name": "7. Structural Biology & CADD",
        "count": 18,
        "color": "#d90429",
        "phyla": ["AI 3D Structure (3)", "Docking/Pockets (6)", "MD Simulation (6)", "Cheminformatics (3)"],
        "nodes": "AlphaFold2, ESMFold, AutoDock Vina, GNINA, OpenMM, MDAnalysis, RDKit..."
    },
    {
        "name": "8. BioPython & Sequence Utilities",
        "count": 23,
        "color": "#3a86ff",
        "phyla": ["Sequence I/O (5)", "Homology/BLAST (5)", "Cloning/Primers (4)", "Circos/Phylo (9)"],
        "nodes": "SeqIO, BLAST, PyHMMER, Primer3, Edlib, PyCircos, DnaFeaturesViewer..."
    },
    {
        "name": "9. Tier-1 Publication Visualizers",
        "count": 20,
        "color": "#fb5607",
        "phyla": ["Expression/DEG (5)", "Single-Cell (3)", "Genomic Tracks (5)", "Microbiome (4)", "CADD (3)"],
        "nodes": "Volcano, Heatmap, UMAP, Manhattan, CoverageTrack, OncoPrint, 3D Render..."
    }
]

# Draw Root
root_x, root_y = 0.0, 4.5
ax_tree.scatter([root_x], [root_y], s=220, color="#222222", zorder=5)
ax_tree.text(root_x - 0.2, root_y, "ROOT\nComfyBioWMS", ha="right", va="center", 
             fontsize=10, fontweight="bold", color="#111111")

y_positions = np.linspace(8.8, 0.2, len(kingdoms))

for idx, (k, y_pos) in enumerate(zip(kingdoms, y_positions)):
    color = k["color"]
    # Draw cladogram branch from root to Kingdom
    ax_tree.plot([root_x, root_x + 0.8, root_x + 0.8, root_x + 1.8], 
                 [root_y, root_y, y_pos, y_pos], 
                 color="#555555", lw=1.8, zorder=2)
    
    # Kingdom Node
    ax_tree.scatter([root_x + 1.8], [y_pos], s=140, color=color, edgecolors="#222222", lw=1.2, zorder=5)
    
    # Kingdom Box
    bbox_props = dict(boxstyle="round,pad=0.35", fc=color, ec="none", alpha=0.92)
    ax_tree.text(root_x + 2.1, y_pos, f"{k['name']} (n={k['count']})", 
                 va="center", fontsize=9.5, fontweight="bold", color="white", bbox=bbox_props)
    
    # Draw sub-branches (Phyla)
    ax_tree.plot([root_x + 5.2, root_x + 5.6], [y_pos, y_pos], color=color, lw=1.5, zorder=2)
    
    phyla_text = " • ".join(k["phyla"])
    ax_tree.text(root_x + 5.8, y_pos + 0.16, f"Phyla: {phyla_text}", 
                 va="center", fontsize=8.2, fontweight="semibold", color="#2b2b2b")
    
    ax_tree.text(root_x + 5.8, y_pos - 0.22, f"Nodes: {k['nodes']}", 
                 va="center", fontsize=7.6, fontstyle="italic", color="#555555")

# ---------------------------------------------------------
# Panel B: Node Census Distribution (Donut Chart)
# ---------------------------------------------------------
ax_donut.set_title("(b) Multi-Domain Node Census (Total 157 Nodes)", 
                   fontsize=12, fontweight="bold", pad=12, loc="left", color="#111111")

counts = [k["count"] for k in kingdoms]
labels = [f"{k['name'].split('.')[1].split('(')[0].strip()} ({k['count']})" for k in kingdoms]
colors = [k["color"] for k in kingdoms]

wedges, texts, autotexts = ax_donut.pie(
    counts, labels=None, autopct="%1.1f%%", startangle=140,
    colors=colors, pctdistance=0.78,
    wedgeprops=dict(width=0.4, edgecolor="white", linewidth=1.5)
)

for autotext in autotexts:
    autotext.set_color("white")
    autotext.set_fontsize(7.5)
    autotext.set_fontweight("bold")

ax_donut.legend(wedges, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), 
                fontsize=8, frameon=False)

ax_donut.text(0, 0, "157\nNodes\n100% Conda", ha="center", va="center", 
              fontsize=11, fontweight="bold", color="#222222")

# ---------------------------------------------------------
# Panel C: Strict Socket Type Contract & Dataflow Wiring
# ---------------------------------------------------------
ax_socket.set_xlim(0, 10)
ax_socket.set_ylim(0, 10)
ax_socket.axis("off")
ax_socket.set_title("(c) Strict Socket Type Contract & Unified Native Dataflow", 
                     fontsize=12, fontweight="bold", pad=12, loc="left", color="#111111")

socket_types = [
    ("FASTQ_PAIR", "#2a9d8f", "Raw Sequenced Reads (Paired/Single)"),
    ("SAM_BAM_INDEXED", "#1d3557", "Coordinate-Sorted & Deduplicated Alignments"),
    ("VCF_FILTERED", "#e76f51", "High-Confidence Genomic Variants (SNVs/InDels)"),
    ("COUNT_MATRIX_CSV", "#8338ec", "Normalized Transcript/Gene Expression Counts"),
    ("PEAK_BED", "#e9c46a", "Enriched Chromatin & Transcription Factor Peaks"),
    ("CONTIG_FASTA", "#457b9d", "De Novo Assembled Genomic Scaffolds"),
    ("IMAGE_TENSOR", "#fb5607", "Native 300+ DPI PyTorch Tensors for Live Preview")
]

y_starts = np.linspace(8.8, 1.2, len(socket_types))

for (sock, col, desc), y_s in zip(socket_types, y_starts):
    # Socket Pill
    rect = patches.FancyBboxPatch((0.2, y_s - 0.4), 3.4, 0.75, boxstyle="round,pad=0.15", 
                                  fc=col, ec="#333333", lw=0.8, alpha=0.9)
    ax_socket.add_patch(rect)
    ax_socket.text(1.9, y_s, sock, ha="center", va="center", color="white", 
                   fontsize=8.5, fontweight="bold")
    
    # Arrow
    ax_socket.annotate("", xy=(4.2, y_s), xytext=(3.7, y_s),
                        arrowprops=dict(arrowstyle="->", color="#444444", lw=1.5))
    
    # Description
    ax_socket.text(4.4, y_s, desc, va="center", fontsize=8.2, color="#222222")

# Summary banner box at bottom of Panel C
banner = patches.FancyBboxPatch((0.2, 0.05), 9.6, 0.75, boxstyle="round,pad=0.15",
                                fc="#f1faee", ec="#a8dadc", lw=1.2)
ax_socket.add_patch(banner)
ax_socket.text(5.0, 0.42, "✓ Zero Type-Mismatch Crashing   ✓ Upstream Smart Caching   ✓ One-Drop PNG Replication",
               ha="center", va="center", fontsize=8.5, fontweight="bold", color="#1d3557")

# Save figure
output_path = "paper_platform/figures/figure_node_phylogeny_tree.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Successfully generated figure at {output_path}")
