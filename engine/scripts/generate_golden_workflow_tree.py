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
        "total_count": 19,
        "nodes": [
            "SampleMetadataValidatorNode (FASTQ/FASTA Integrity Check)",
            "FastpTrimNode (Paired-End Quality & Adapter Trim)",
            "SpadesAssembleNode (Multi k-mer De Bruijn Graph Assembly)",
            "QuastQcNode (Contig Metrics, N50 & Misassembly QC)",
            "AssemblyVisualizationNode & Report (Bandage Graph & Summary)"
        ]
    },
    {
        "name": "Epigenomics & Functional Screening",
        "total_count": 17,
        "nodes": [
            "SampleMetadataValidatorNode (Paired FASTQ & Ref Genome Check)",
            "FastpTrimNode (Poly-G & Adapter Read Trimming)",
            "BwaMem2AlignNode (Paired-End Alignment Engine)",
            "MarkDuplicatesNode & Filter (Optical Duplicates & MAPQ)",
            "Macs3PeakCallingNode (BAMPE Open Chromatin Peak Caller)",
            "AtacPeakVisualizationNode & Report (Signal Track & Summary)"
        ]
    },
    {
        "name": "Metagenome & Pathogen Virome",
        "total_count": 18,
        "nodes": [
            "SampleMetadataValidatorNode (FASTQ & DB Integrity Check)",
            "FastpTrimNode (Quality & Length Filtering)",
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
        "total_count": 7,
        "nodes": [
            "SampleMetadataValidatorNode (FASTQ & Ref Genome Validator)",
            "BwaMem2AlignNode (Burrows-Wheeler DNA Sequence Aligner)",
            "MarkDuplicatesNode (Picard Optical Duplicate Tagger)",
            "BcftoolsCallNode & Filter (Multiallelic Variant Calling & Filter)",
            "VariantVisualizationNode & Report (SNV/Indel Spectrum & Ti/Tv)"
        ]
    }
]

def render_aligned_tree(out_path='figures/fig2_node_domain_distribution.png'):
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    total_ecosystem_nodes = sum(d["total_count"] for d in GOLDEN_WORKFLOW_DATA)

    fig_w, fig_h = 18.5, 13.5
    dpi = 300

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

    x_root = 0.075
    x_domain_in = 0.23
    x_domain_text = 0.242
    y_root = norm_y(max_raw_y / 2.0)

    curve_color_root = '#94a3b8'   # Slate 400
    curve_color_leaf = '#cbd5e1'   # Slate 300
    circle_border = '#2563eb'      # Royal Blue
    circle_fill = '#e0f2fe'        # Light Sky Blue
    root_border = '#1d4ed8'        # Deep Blue
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
    ax.scatter([x_root], [y_root], s=200, facecolor=root_fill, edgecolor=root_border, linewidth=2.2, zorder=5)
    ax.text(x_root - 0.012, y_root + 0.016, "ComfyBIOWMS", 
            ha='right', va='center', fontsize=13.0, fontweight='bold', color='#0f172a', zorder=6)
    ax.text(x_root - 0.012, y_root - 0.016, f"(10 Domains · {total_ecosystem_nodes} Nodes)", 
            ha='right', va='center', fontsize=9.4, fontweight='bold', color='#475569', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

    # 3. Draw Domain Nodes & Text and Measure Bounding Boxes
    domain_text_objs = []
    for d in slot_allocations:
        ax.scatter([x_domain_in], [d["domain_y"]], s=115, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.8, zorder=5)
        cat_text = f"{d['name']} ({d['total_count']} nodes)"
        t_obj = ax.text(x_domain_text, d["domain_y"], cat_text, 
                        ha='left', va='center', fontsize=10.2, fontweight='bold', color='#1e293b', zorder=6)
        domain_text_objs.append((d, t_obj))

    # Canvas draw to get true text bounding boxes
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv_trans = ax.transAxes.inverted()

    # Find the maximum right edge of the longest domain text
    max_text_right = max(t_obj.get_window_extent(renderer).transformed(inv_trans).x1 for _, t_obj in domain_text_objs)
    print(f"Max domain text right edge: {max_text_right:.4f}")

    # Set aligned branch starting point based on the longest domain name (tight margin of +0.016)
    x_aligned_branch = max_text_right + 0.016
    x_aligned_leaf = x_aligned_branch + 0.088
    print(f"Aligned branch point: {x_aligned_branch:.4f}, Aligned leaf column: {x_aligned_leaf:.4f}")

    # 4. Draw horizontal connector from each domain text to the aligned branch point, and fan out to aligned leaf column
    for d, t_obj in domain_text_objs:
        bbox = t_obj.get_window_extent(renderer).transformed(inv_trans)
        text_right = bbox.x1
        
        # Horizontal line from the end of the text to the aligned branch point
        ax.plot([text_right + 0.003, x_aligned_branch], [d["domain_y"], d["domain_y"]], color='#cbd5e1', lw=1.2, zorder=2)
        
        # Draw Bézier curves fanning out from the unified branch point to the unified leaf column
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_aligned_branch, d["domain_y"], x_aligned_leaf, y_l, color=curve_color_leaf, lw=1.0, alpha=0.85)

        # Draw Leaf circles & Labels
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_aligned_leaf], [y_l], s=45, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.3, zorder=5)
            ax.text(x_aligned_leaf + 0.009, y_l, node_label, 
                    ha='left', va='center', fontsize=8.5, color='#334155', zorder=6)

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Successfully generated aligned tree at {out_path}")

if __name__ == '__main__':
    render_aligned_tree('figures/fig2_node_domain_distribution.png')
    import shutil
    shutil.copyfile('figures/fig2_node_domain_distribution.png', 'figures/fig2_golden_workflow_tree.png')
