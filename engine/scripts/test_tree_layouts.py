import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

GOLDEN_WORKFLOW_DATA = [
    {
        "name": "Biopython & Sequence Ops",
        "total_count": 26,
        "nodes": [
            "BiopythonSeqIONode (Multi-Format Sequence Parser)",
            "BiopythonBlastNode (NCBI Local/Remote BLAST Search)",
            "Primer3DesignNode (Automated PCR Primer Designer)",
            "BiopythonBioPDBNode (3D Macromolecular Structure Parser)",
            "LogomakerVisualizerNode (Sequence Logo & Motif Visualizer)"
        ]
    },
    {
        "name": "Single-Cell & Spatial Omics",
        "total_count": 19,
        "nodes": [
            "ScanpyQCNode (Mitochondrial & Gene Count Quality Filter)",
            "ScanpyNormalizeNode (Total Counts & Log1p Normalization)",
            "ScanpyClusterNode (Leiden & Louvain Community Detection)",
            "ScanpyMarkerGenesNode (Differential Marker Gene Ranking)",
            "ScRNAVisualizationNode & Report (UMAP Plot & QC Report)"
        ]
    },
    {
        "name": "Genome Assembly & Long-Read",
        "total_count": 21,
        "nodes": [
            "AssemblyInputValidatorNode (FASTQ/FASTA Integrity Check)",
            "AssemblyFastpTrimNode (Paired-End Quality & Adapter Trim)",
            "SpadesAssembleNode (Multi k-mer De Bruijn Graph Assembly)",
            "QuastQcNode (Contig Metrics, N50 & Misassembly QC)",
            "AssemblyVisualizationNode & Report (Bandage Graph & Summary)"
        ]
    },
    {
        "name": "Epigenomics & Functional Screening",
        "total_count": 22,
        "nodes": [
            "AtacInputValidatorNode (Paired FASTQ & Ref Genome Check)",
            "AtacFastpTrimNode (Poly-G & Adapter Read Trimming)",
            "AtacBwaMem2AlignNode (Paired-End Alignment Engine)",
            "AtacMarkDuplicatesNode & Filter (Optical Duplicates & MAPQ)",
            "Macs3PeakCallingNode (BAMPE Open Chromatin Peak Caller)",
            "AtacPeakVisualizationNode & Report (Signal Track & Summary)"
        ]
    },
    {
        "name": "Metagenome & Pathogen Virome",
        "total_count": 20,
        "nodes": [
            "MetagenomeInputValidatorNode (FASTQ & DB Integrity Check)",
            "MetagenomeFastpTrimNode (Quality & Length Filtering)",
            "Kraken2ClassifyNode (Exact k-mer Taxonomic Classifier)",
            "BrackenAbundanceNode (Bayesian Abundance Re-Estimation)",
            "MetagenomeVisualizationNode & Report (Krona / Stacked & Summary)"
        ]
    },
    {
        "name": "Bulk RNA-Seq & Core Pipeline",
        "total_count": 20,
        "nodes": [
            "SampleMetadataValidatorNode (Samplesheet & Path Validator)",
            "FastpTrimNode (Automated Quality & Poly-G Read Trimming)",
            "SalmonIndexNode & SalmonQuantNode (Quasi-Mapping Quantifier)",
            "TximportNode (Gene-Level Count Matrix Summarizer)",
            "DESeq2AnalysisNode (Negative Binomial GLM Differential Expression)",
            "DESeq2VisualizationNode & Report (Volcano/PCA & Audit Report)"
        ]
    },
    {
        "name": "Publication Quality Visualizers",
        "total_count": 20,
        "nodes": [
            "VolcanoPlotVisualizerNode (DEG Significance Volcano Plot)",
            "ManhattanPlotVisualizerNode (GWAS P-Value Manhattan Plot)",
            "ClustermapHeatmapVisualizerNode (Hierarchical Clustermap)",
            "UmapScatterVisualizerNode (Single-Cell UMAP Embeddings)",
            "GseaEnrichmentPlotVisualizerNode (Running Enrichment Score Curves)"
        ]
    },
    {
        "name": "CADD & Structural Biology",
        "total_count": 18,
        "nodes": [
            "ColabFoldAlphaFoldNode (MMseqs2 Accelerated 3D Structure)",
            "ESMFoldNode (Large Protein Language Model Folding)",
            "DiffDockPredictNode (Generative Blind Ligand Docking)",
            "AutoDockVinaDockingNode (Physics-Based Molecular Docking)",
            "OpenMMSimulationNode (GPU Molecular Dynamics Simulation)"
        ]
    },
    {
        "name": "Proteomics & Metabolomics",
        "total_count": 15,
        "nodes": [
            "PyOpenMSFeatureNode (LC-MS 2D/3D Feature Detection Engine)",
            "DiaNNQuantNode (Neural Network DIA Proteomics Quantifier)",
            "MSFraggerSearchNode (Ultra-Fast Peptide Database Search)",
            "MatchmsSpectrumNode (Mass Spectral Similarity Engine)",
            "SiriusStructureNode (Metabolite MS/MS Molecular Formula)"
        ]
    },
    {
        "name": "DNA Variant Calling",
        "total_count": 8,
        "nodes": [
            "VariantInputValidatorNode (FASTQ & Ref Genome Validator)",
            "BwaMem2AlignNode (Burrows-Wheeler DNA Sequence Aligner)",
            "MarkDuplicatesNode (Picard Optical Duplicate Tagger)",
            "BcftoolsCallNode & Filter (Multiallelic Variant Calling & Filter)",
            "VariantVisualizationNode & Report (SNV/Indel Spectrum & Ti/Tv)"
        ]
    }
]

def test_layouts():
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    fig_w, fig_h = 19.5, 13.5
    dpi = 300

    # Layout Option 1: Clean Sequential Column Layout (Zero Overlap Guaranteed)
    # [Root: 0.08] ====> [Cat Circle: 0.28] Cat Label [Cat Branch: 0.52] ====> [Leaf Circle: 0.62] Leaf Label
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    domain_gap = 1.35
    slot_allocations = []
    current_y = 0.0

    for d in reversed(GOLDEN_WORKFLOW_DATA):
        n_leaves = len(d["nodes"])
        leaf_slots = [current_y + (n_leaves - 1 - i) * 1.0 for i in range(n_leaves)]
        domain_center_y = np.mean(leaf_slots)
        slot_allocations.append({
            "name": d["name"],
            "total_count": d["total_count"],
            "nodes": d["nodes"],
            "leaf_ys_raw": leaf_slots,
            "domain_y_raw": domain_center_y
        })
        current_y += n_leaves * 1.0 + domain_gap

    slot_allocations = list(reversed(slot_allocations))
    max_raw_y = current_y - domain_gap

    y_min_target, y_max_target = 0.045, 0.955
    def norm_y(val):
        return y_min_target + (val / max_raw_y) * (y_max_target - y_min_target)

    for d in slot_allocations:
        d["leaf_ys"] = [norm_y(y) for y in d["leaf_ys_raw"]]
        d["domain_y"] = norm_y(d["domain_y_raw"])

    x_root = 0.06
    x_domain_in = 0.24     # Where Root curve connects to Domain
    x_domain_text = 0.255  # Domain text start
    x_domain_out = 0.52    # Where Domain curves branch to Leaves
    x_leaf = 0.62          # Leaf circles
    y_root = norm_y(max_raw_y / 2.0)

    curve_color_root = '#94a3b8'   # Slate 400
    curve_color_leaf = '#cbd5e1'   # Slate 300
    circle_border = '#2563eb'      # Royal Blue
    circle_fill = '#e0f2fe'        # Light Sky Blue
    root_border = '#1d4ed8'        # Deep Blue
    root_fill = '#bfdbfe'          # Soft Blue

    def draw_smooth_spline(x1, y1, x2, y2, color, lw=1.2, alpha=0.9):
        dx = x2 - x1
        c1 = (x1 + dx * 0.5, y1)
        c2 = (x2 - dx * 0.5, y2)
        verts = [(x1, y1), c1, c2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        patch = patches.PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # 1. Connect Root to Domain Input Nodes
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain_in, d["domain_y"], color=curve_color_root, lw=1.4, alpha=0.9)

    # 2. Connect Domain Output to Leaf Nodes
    for d in slot_allocations:
        # Subtle horizontal link connecting domain text to domain output circle
        ax.plot([x_domain_in, x_domain_out], [d["domain_y"], d["domain_y"]], color='#e2e8f0', lw=1.2, zorder=1)
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_domain_out, d["domain_y"], x_leaf, y_l, color=curve_color_leaf, lw=1.0, alpha=0.85)

    # 3. Draw Root Node & Label
    ax.scatter([x_root], [y_root], s=190, facecolor=root_fill, edgecolor=root_border, linewidth=2.2, zorder=5)
    ax.text(x_root - 0.012, y_root + 0.016, "ComfyBIOWMS", 
            ha='right', va='center', fontsize=12.5, fontweight='bold', color='#0f172a', zorder=6)
    ax.text(x_root - 0.012, y_root - 0.016, "(10 Domains · 189 Nodes)", 
            ha='right', va='center', fontsize=9.2, fontweight='bold', color='#475569', zorder=6)

    # 4. Draw Domain Nodes & Labels
    for d in slot_allocations:
        # Input Circle
        ax.scatter([x_domain_in], [d["domain_y"]], s=100, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.8, zorder=5)
        # Output Circle
        ax.scatter([x_domain_out], [d["domain_y"]], s=70, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.4, zorder=5)
        
        # Domain Label positioned cleanly between Input and Output circles
        cat_text = f"{d['name']} ({d['total_count']} nodes)"
        ax.text(x_domain_text, d["domain_y"], cat_text, 
                ha='left', va='center', fontsize=9.8, fontweight='bold', color='#1e293b', zorder=6,
                bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='none'))

    # 5. Draw Leaf Nodes & Labels
    for d in slot_allocations:
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf], [y_l], s=44, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.3, zorder=5)
            ax.text(x_leaf + 0.010, y_l, node_label, 
                    ha='left', va='center', fontsize=8.2, color='#334155', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    out_opt1 = 'figures/fig2_test_clean_columns.png'
    plt.savefig(out_opt1, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved Option 1 to {out_opt1}")

    # Layout Option 2: Pure Reference Style with Curve Masking & Left Category Circle
    # Domain Circle at x=0.48, Category Label right-aligned at x=0.465, with white halo/masked root curves terminating at text boundary or circle
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    x_root = 0.08
    x_domain_circle = 0.44
    x_leaf = 0.58
    
    # 1. Connect Root to Domain Circle
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain_circle, d["domain_y"], color=curve_color_root, lw=1.4, alpha=0.9)

    # 2. Connect Domain Circle to Leaves
    for d in slot_allocations:
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_domain_circle, d["domain_y"], x_leaf, y_l, color=curve_color_leaf, lw=1.0, alpha=0.85)

    # 3. Root Node
    ax.scatter([x_root], [y_root], s=190, facecolor=root_fill, edgecolor=root_border, linewidth=2.2, zorder=5)
    ax.text(x_root - 0.012, y_root + 0.016, "ComfyBIOWMS", 
            ha='right', va='center', fontsize=12.5, fontweight='bold', color='#0f172a', zorder=6)
    ax.text(x_root - 0.012, y_root - 0.016, "(10 Domains · 189 Nodes)", 
            ha='right', va='center', fontsize=9.2, fontweight='bold', color='#475569', zorder=6)

    # 4. Domain Nodes & Labels with Clean White Halo Bounding Box
    for d in slot_allocations:
        ax.scatter([x_domain_circle], [d["domain_y"]], s=115, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.8, zorder=7)
        cat_text = f"{d['name']} ({d['total_count']} nodes)"
        # Using a crisp white background box to cleanly mask any crossing root curves!
        ax.text(x_domain_circle - 0.012, d["domain_y"], cat_text, 
                ha='right', va='center', fontsize=10.0, fontweight='bold', color='#1e293b', zorder=6,
                bbox=dict(boxstyle='round,pad=0.25,rounding_size=0.2', facecolor='white', edgecolor='#e2e8f0', lw=0.8))

    # 5. Leaf Nodes
    for d in slot_allocations:
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf], [y_l], s=44, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.3, zorder=5)
            ax.text(x_leaf + 0.010, y_l, node_label, 
                    ha='left', va='center', fontsize=8.2, color='#334155', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    out_opt2 = 'figures/fig2_test_masked_reference.png'
    plt.savefig(out_opt2, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved Option 2 to {out_opt2}")

if __name__ == '__main__':
    test_layouts()
