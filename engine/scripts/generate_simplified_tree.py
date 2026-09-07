import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# Curated Representative Structure for ComfyBIOWMS Tree
SIMPLIFIED_DOMAIN_DATA = [
    {
        "name": "Biopython & Sequence Ops",
        "total_count": 26,
        "representative_nodes": [
            "SeqIO & AlignIO (Multi-Format Sequence & Alignment Parsers)",
            "BLAST & Entrez (NCBI Remote Search & Query Retrieval)",
            "Primer3 (Automated PCR Primer Design Engine)",
            "BioPDB & Biotite (3D Macromolecular Parsing & RMSD)",
            "Logomaker & PyCircos (Sequence Logo & Circular Visualizers)"
        ]
    },
    {
        "name": "Single-Cell & Spatial Omics",
        "total_count": 19,
        "representative_nodes": [
            "Scanpy (Single-Cell QC, Normalization & Leiden Clustering)",
            "scVI & SCANVI (Deep Generative Latent & Cell Type Annotations)",
            "scVelo & CellRank (RNA Velocity Kinetics & Lineage Trajectory)",
            "Squidpy & Tangram (Spatial Transcriptomics & Cell Mapping)",
            "CellPhoneDB & PySCENIC (Ligand-Receptor & Regulon Networks)"
        ]
    },
    {
        "name": "Genome Assembly & Long-Read",
        "total_count": 21,
        "representative_nodes": [
            "SPAdes & QUAST (De Novo Assembly & Quality Metrics)",
            "Flye & Hifiasm (Long-Read Repeat & PacBio HiFi Assemblers)",
            "Sniffles2 & CuteSV (Long-Read Structural Variant Callers)",
            "DeepVariant (CNN-Based High-Precision Germline Caller)",
            "Pysam & pybedtools (Genomic Interval & BAM Processing)"
        ]
    },
    {
        "name": "Epigenomics & Functional Screening",
        "total_count": 22,
        "representative_nodes": [
            "MACS3 (BAMPE Open Chromatin & Peak Calling)",
            "deepTools (Coverage Profiler, Metagene & Heatmaps)",
            "TOBIAS (Transcription Factor Footprinting & Occupancy)",
            "CRISPResso2 & MAGeCK (CRISPR Editing & Screen Analysis)",
            "Cooler & Chromosight (Hi-C Contact Matrices & Loop Detection)"
        ]
    },
    {
        "name": "Metagenome & Microbiome Virome",
        "total_count": 20,
        "representative_nodes": [
            "Kraken2 & Bracken (Exact k-mer Taxon Classifier & Abundance)",
            "DADA2 (Error-Corrected Amplicon ASV Denoiser)",
            "HUMAnN3 & MetaPhlAn4 (Functional Metabolic Pathway Profilers)",
            "geNomad & VirSorter2 (Viral & Plasmid Mining Engines)",
            "Prokka & Bakta (Rapid Prokaryotic & Bacterial Annotators)"
        ]
    },
    {
        "name": "Bulk RNA-Seq & Core Pipeline",
        "total_count": 20,
        "representative_nodes": [
            "Fastp & FastQC (Automated Quality Assessment & Trimming)",
            "Salmon & Tximport (Decoy-Aware Quant & Count Aggregator)",
            "DESeq2 & edgeR (Negative Binomial Differential Expression)",
            "10x Cell Ranger & STARsolo (High-Throughput Droplet Processors)",
            "StringTie2 & FLAIR (Transcript Assembly & Isoform Analysis)"
        ]
    },
    {
        "name": "Publication Quality Visualizers",
        "total_count": 20,
        "representative_nodes": [
            "Volcano & MA Plots (Differential Expression Significance)",
            "Manhattan & Q-Q Plots (Genome-Wide Association & Quantiles)",
            "UMAP & Clustermap (Cell Embeddings & Hierarchical Heatmaps)",
            "GSEA & OncoPrint (Gene Set Enrichment & Genomic Alterations)",
            "Kaplan-Meier & PCoA (Survival Analysis & Beta Diversity)"
        ]
    },
    {
        "name": "CADD & Structural Biology",
        "total_count": 18,
        "representative_nodes": [
            "AlphaFold2 & ColabFold (MMseqs2-Accelerated 3D Structure)",
            "ESMFold (Large Protein Language Model Folding)",
            "DiffDock & GNINA (Generative Blind & Deep Learning Docking)",
            "AutoDock Vina & Smina (Physics-Based Scoring & Docking)",
            "OpenMM & MDAnalysis (GPU Molecular Dynamics & Trajectories)"
        ]
    },
    {
        "name": "Proteomics & Metabolomics",
        "total_count": 15,
        "representative_nodes": [
            "PyOpenMS & Pyteomics (LC-MS 2D/3D Feature Detection & Parsers)",
            "DIA-NN & MaxQuant (Neural Network DIA & Label-Free Proteomics)",
            "MSFragger (Ultra-Fast High-Throughput Peptide Search)",
            "Matchms & Spec2Vec (Spectral Similarity & NLP Word Matching)",
            "SIRIUS & MS-DIAL (Metabolite MS/MS Formula & Lipidomics)"
        ]
    },
    {
        "name": "DNA Variant Calling",
        "total_count": 8,
        "representative_nodes": [
            "BWA-MEM2 (Burrows-Wheeler Short-Read DNA Aligner)",
            "MarkDuplicates (Picard Optical & PCR Duplicate Tagger)",
            "BCFtools Call & Filter (Multiallelic Variant Calling & Filtering)",
            "Variant Visualization (SNV/Indel Spectrum, Ti/Tv & Quality QC)"
        ]
    }
]

def render_simplified_tree(out_path='figures/fig2_node_domain_distribution.png'):
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    total_repr_leaves = sum(len(d["representative_nodes"]) for d in SIMPLIFIED_DOMAIN_DATA)
    total_all_nodes = sum(d["total_count"] for d in SIMPLIFIED_DOMAIN_DATA)
    print(f"Rendering simplified tree: 10 domains, {total_repr_leaves} representative nodes (representing all {total_all_nodes} nodes)...")

    # Dimensions: 17 inches wide x 13 inches tall for standard paper / slide landscape layout
    fig_w, fig_h = 17.5, 12.5
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    # Vertical spacing allocation
    domain_gap = 1.3
    slot_allocations = []
    current_y = 0.0

    for d in reversed(SIMPLIFIED_DOMAIN_DATA):
        n_leaves = len(d["representative_nodes"])
        leaf_slots = [current_y + i * 1.0 for i in range(n_leaves)]
        domain_center_y = np.mean(leaf_slots)
        slot_allocations.append({
            "name": d["name"],
            "total_count": d["total_count"],
            "nodes": d["representative_nodes"],
            "leaf_ys_raw": leaf_slots,
            "domain_y_raw": domain_center_y
        })
        current_y += n_leaves * 1.0 + domain_gap

    slot_allocations = list(reversed(slot_allocations))
    max_raw_y = current_y - domain_gap

    y_min_target, y_max_target = 0.05, 0.95
    def norm_y(val):
        return y_min_target + (val / max_raw_y) * (y_max_target - y_min_target)

    for d in slot_allocations:
        d["leaf_ys"] = [norm_y(y) for y in d["leaf_ys_raw"]]
        d["domain_y"] = norm_y(d["domain_y_raw"])

    # X coordinates (Root -> Domains -> Leaves)
    x_root = 0.10
    x_domain = 0.40
    x_leaf = 0.58
    y_root = norm_y(max_raw_y / 2.0)

    # Reference-matched aesthetic styling
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

    # 1. Draw smooth Bézier curves from Root to Domains
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain, d["domain_y"], color=curve_color_root, lw=1.5, alpha=0.9)

    # 2. Draw smooth Bézier curves from Domains to Representative Leaves
    for d in slot_allocations:
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_domain, d["domain_y"], x_leaf, y_l, color=curve_color_leaf, lw=1.0, alpha=0.85)

    # 3. Draw Root Node
    ax.scatter([x_root], [y_root], s=180, facecolor=root_fill, edgecolor=root_border, linewidth=2.2, zorder=5)
    ax.text(x_root - 0.015, y_root + 0.015, "ComfyBIOWMS", 
            ha='right', va='center', fontsize=12.5, fontweight='bold', color='#0f172a', zorder=6)
    ax.text(x_root - 0.015, y_root - 0.015, f"(10 Domains · {total_all_nodes} Nodes)", 
            ha='right', va='center', fontsize=9.5, fontweight='semibold', color='#475569', zorder=6)

    # 4. Draw Domain Nodes & Labels
    for d in slot_allocations:
        ax.scatter([x_domain], [d["domain_y"]], s=110, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.8, zorder=5)
        # Category label to left of circle
        cat_text = f"{d['name']} ({d['total_count']} nodes)"
        ax.text(x_domain - 0.012, d["domain_y"], cat_text, 
                ha='right', va='center', fontsize=10.0, fontweight='bold', color='#1e293b', zorder=6)

    # 5. Draw Representative Leaf Nodes & Labels
    for d in slot_allocations:
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf], [y_l], s=42, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.3, zorder=5)
            # Leaf label to right of circle
            ax.text(x_leaf + 0.010, y_l, node_label, 
                    ha='left', va='center', fontsize=8.2, color='#334155', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved simplified tree figure to {out_path}")

if __name__ == '__main__':
    render_simplified_tree('figures/fig2_node_domain_distribution.png')
    # Also save a dedicated copy
    import shutil
    shutil.copyfile('figures/fig2_node_domain_distribution.png', 'figures/fig2_node_domain_tree_simplified.png')
