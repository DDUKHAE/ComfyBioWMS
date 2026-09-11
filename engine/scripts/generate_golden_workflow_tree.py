"""Generate Figure 2: Hierarchical Tree Diagram of ComfyBIOWMS 10 Functional Modules (98 Nodes).

This script constructs a publication-grade hierarchical tree visualization (Figure 2)
faithfully representing all 98 active custom nodes registered in ComfyBIOWMS,
categorized into 10 bioinformatic functional modules with exact node counts.
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

DOMAINS_DATA = [
    {
        "name": "Quality Control & Read Preprocessing",
        "total_count": 17,
        "color": "#0d9488",      # Teal
        "bg_color": "#ccfbf1",
        "nodes": [
            "FastQC & MultiQC",
            "Fastp, Trimmomatic & TrimGalore",
            "SortMeRNA & BBSplit",
            "CatFastq, InferStrandedness & UmiToolsExtract",
            "Preseq, QualimapRNASeq & DupRadar",
            "DESeq2SampleQC & RSeQC (GeneBodyCoverage, InferExperiment, JunctionSaturation)",
        ]
    },
    {
        "name": "Sequence Alignment & Coordinate Processing",
        "total_count": 13,
        "color": "#2563eb",      # Blue
        "bg_color": "#dbeafe",
        "nodes": [
            "STARAlignReads & STARGenomeGenerate",
            "HISAT2Align & HISAT2Build",
            "Bowtie2Align & Bowtie2Build",
            "BwaMem2Align & BwaMem2Index",
            "SamtoolsSort, SamtoolsIndex & SamtoolsMarkdup",
            "PicardMarkDuplicates & UmiToolsDedup",
        ]
    },
    {
        "name": "De Novo Assembly & Structure Assessment",
        "total_count": 5,
        "color": "#d97706",      # Amber
        "bg_color": "#fef3c7",
        "nodes": [
            "Spades",
            "FlyeAssemble & HifiasmAssemble",
            "StringTie",
            "Quast",
        ]
    },
    {
        "name": "Expression Quantification & Statistical Modeling",
        "total_count": 8,
        "color": "#4f46e5",      # Indigo
        "bg_color": "#e0e7ff",
        "nodes": [
            "SalmonQuantReads, SalmonQuantAlignment & SalmonIndex",
            "KallistoQuant & RSEMCalculateExpression",
            "Tximport & DESeq2",
            "GSEAPathway",
        ]
    },
    {
        "name": "Variant Calling & Genomic Interval Arithmetic",
        "total_count": 8,
        "color": "#059669",      # Emerald
        "bg_color": "#d1fae5",
        "nodes": [
            "BcftoolsMpileup, BcftoolsCall & BcftoolsFilter",
            "Cyvcf2Stats, Mosdepth & PysamStats",
            "PybedtoolsIntersect & SeqKitStats",
        ]
    },
    {
        "name": "Taxonomic Profiling & Functional Annotation",
        "total_count": 7,
        "color": "#65a30d",      # Lime
        "bg_color": "#ecfccb",
        "nodes": [
            "Kraken2Classify, Bracken & SylphProfile",
            "MetaPhlAn & HUMAnN",
            "Prokka & Bakta",
        ]
    },
    {
        "name": "Single-Cell Processing & Manifold Learning",
        "total_count": 4,
        "color": "#9333ea",      # Purple
        "bg_color": "#f3e8ff",
        "nodes": [
            "AnnDataInspect",
            "ScanpyQC & ScanpyNormalize",
            "ScanpyCluster",
        ]
    },
    {
        "name": "3D Structure Prediction & Molecular Docking",
        "total_count": 4,
        "color": "#e11d48",      # Rose
        "bg_color": "#ffe4e6",
        "nodes": [
            "ColabFold & ESMFold",
            "AutoDockVina & RDKitDescriptor",
        ]
    },
    {
        "name": "Sequence Analytics & Workflow Utilities",
        "total_count": 9,
        "color": "#475569",      # Slate
        "bg_color": "#f1f5f9",
        "nodes": [
            "BiopythonSeqIOStats & BiopythonGCContent",
            "BiopythonPairwiseAlign & BiopythonAlignmentStats",
            "BiopythonProtParam & BiopythonRestrictionDigest",
            "BiopythonSeqFilter & BiopythonSeqTransform",
            "JoinPaths",
        ]
    },
    {
        "name": "Multi-Omics Data Visualization & Genomic Tracks",
        "total_count": 23,
        "color": "#0284c7",      # Sky
        "bg_color": "#e0f2fe",
        "nodes": [
            "VolcanoPlot, MaPlot & ClustermapHeatmap",
            "ManhattanPlot, QqPlot & LinkageDisequilibrium",
            "SyntenyGenome & OncoPrint",
            "MicrobiomeStackedBar, PcoaScatter & PhylogeneticTree",
            "UmapScatter, SankeyCellFate & SpatialTissueOverlay",
            "ProteinLigandInteraction, RamachandranPlot & MdTrajectoryPlotter",
            "GseaEnrichmentPlot & KaplanMeierSurvival",
            "DeeptoolsMatrix, ChipAtacCoverageProfile, BedGraphToBigWig & BedtoolsGenomeCoverage",
        ]
    }
]


def render_aligned_tree(out_path='figures/fig2_node_domain_distribution.png'):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    total_ecosystem_nodes = sum(d["total_count"] for d in DOMAINS_DATA)
    total_domains = len(DOMAINS_DATA)

    fig_w, fig_h = 16.0, 13.0
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    domain_gap = 1.35
    slot_allocations = []
    current_y = 0.0

    # Layout domains from bottom to top
    for d in reversed(DOMAINS_DATA):
        n_leaves = len(d["nodes"])
        leaf_slots = [current_y + (n_leaves - 1 - i) * 1.0 for i in range(n_leaves)]
        domain_center_y = np.mean(leaf_slots)
        slot_allocations.append({
            "name": d["name"],
            "total_count": d["total_count"],
            "color": d["color"],
            "bg_color": d["bg_color"],
            "nodes": d["nodes"],
            "leaf_ys_raw": leaf_slots,
            "domain_y_raw": domain_center_y
        })
        current_y += n_leaves * 1.0 + domain_gap

    slot_allocations = list(reversed(slot_allocations))
    max_raw_y = current_y - domain_gap

    y_min_target, y_max_target = 0.038, 0.962
    def norm_y(val):
        return y_min_target + (val / max_raw_y) * (y_max_target - y_min_target)

    for d in slot_allocations:
        d["leaf_ys"] = [norm_y(y) for y in d["leaf_ys_raw"]]
        d["domain_y"] = norm_y(d["domain_y_raw"])

    x_root = 0.115
    x_domain_in = 0.285
    x_domain_text = 0.297
    y_root = norm_y(max_raw_y / 2.0)

    curve_color_root = '#94a3b8'   # Slate 400
    curve_color_leaf = '#cbd5e1'   # Slate 300
    root_border = '#0f172a'        # Slate 900
    root_fill = '#bfdbfe'          # Soft Blue

    def draw_smooth_spline(x1, y1, x2, y2, color, lw=1.2, alpha=0.9):
        dx = x2 - x1
        c1 = (x1 + dx * 0.48, y1)
        c2 = (x2 - dx * 0.48, y2)
        verts = [(x1, y1), c1, c2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        patch = patches.PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # 1. Connect Root to Domain Input Circles
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain_in, d["domain_y"], color=curve_color_root, lw=1.5, alpha=0.9)

    # 2. Draw Root Node & Text
    ax.scatter([x_root], [y_root], s=250, facecolor=root_fill, edgecolor=root_border, linewidth=2.4, zorder=5)
    ax.text(x_root - 0.012, y_root + 0.017, "ComfyBIOWMS", 
            ha='right', va='center', fontsize=13.5, fontweight='bold', color='#0f172a', zorder=6)
    ax.text(x_root - 0.012, y_root - 0.008, f"{total_domains} Functional Modules", 
            ha='right', va='center', fontsize=10.0, fontweight='bold', color='#1e293b', zorder=6)
    ax.text(x_root - 0.012, y_root - 0.024, f"({total_ecosystem_nodes} Active Registered Nodes)", 
            ha='right', va='center', fontsize=8.8, color='#475569', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

    # 3. Draw Domain Nodes & Text and Measure Bounding Boxes
    domain_text_objs = []
    for d in slot_allocations:
        ax.scatter([x_domain_in], [d["domain_y"]], s=115, facecolor=d["bg_color"], edgecolor=d["color"], linewidth=1.8, zorder=5)
        cat_text = f"{d['name']} ({d['total_count']} nodes)"
        t_obj = ax.text(x_domain_text, d["domain_y"], cat_text, 
                        ha='left', va='center', fontsize=9.6, fontweight='bold', color='#0f172a', zorder=6)
        domain_text_objs.append((d, t_obj))

    # Canvas draw to get true text bounding boxes
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv_trans = ax.transAxes.inverted()

    # Find maximum right edge of domain text
    max_text_right = max(t_obj.get_window_extent(renderer).transformed(inv_trans).x1 for _, t_obj in domain_text_objs)
    x_aligned_branch = max_text_right + 0.016
    x_aligned_leaf = x_aligned_branch + 0.050

    # 4. Draw horizontal connector and fanned out leaf nodes
    for d, t_obj in domain_text_objs:
        bbox = t_obj.get_window_extent(renderer).transformed(inv_trans)
        text_right = bbox.x1
        
        # Line from domain text to branch point
        ax.plot([text_right + 0.004, x_aligned_branch], [d["domain_y"], d["domain_y"]], color='#cbd5e1', lw=1.2, zorder=2)
        
        # Splines from branch point to leaves
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_aligned_branch, d["domain_y"], x_aligned_leaf, y_l, color=curve_color_leaf, lw=1.0, alpha=0.85)

        # Draw Leaf circles & Clean Node Labels
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_aligned_leaf], [y_l], s=42, facecolor='#f8fafc', edgecolor=d["color"], linewidth=1.4, zorder=5)
            ax.text(x_aligned_leaf + 0.009, y_l, node_label, 
                    ha='left', va='center', fontsize=8.8, color='#1e293b', fontweight='normal', zorder=6)

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Successfully generated functional tree at {out_path}")


if __name__ == '__main__':
    render_aligned_tree('figures/fig2_node_domain_distribution.png')
    shutil.copyfile('figures/fig2_node_domain_distribution.png', 'figures/fig2_golden_workflow_tree.png')
    print("Copied to figures/fig2_golden_workflow_tree.png")
