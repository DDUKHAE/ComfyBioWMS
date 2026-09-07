import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# Complete 10 Domain data with 189 toolchain nodes
DOMAIN_DATA = [
    {
        "name": "Biopython & Sequence Ops",
        "nodes": [
            "SeqIO (Multi-Format Sequence Parser)",
            "SeqTranscribe (DNA/RNA/Protein Translation)",
            "AlignIO (Multiple Sequence Alignment)",
            "BLAST (NCBI Remote & Local Search)",
            "Entrez (NCBI E-Utilities Query Engine)",
            "BioPDB (3D Macromolecular Structure Parser)",
            "Restriction (Enzyme Digestion Simulator)",
            "BioMotifs (PSSM & IUPAC Motif Search)",
            "Logomaker Visualizer (Sequence Logo Plotter)",
            "ProtParam (Physicochemical Property Calculator)",
            "Biotite (Structural Superimposition & RMSD)",
            "Squiggle (Raw Nanopore Signal Visualizer)",
            "PyHMMER (Hidden Markov Model Profile Search)",
            "DNA Features Viewer (Genomic Feature Annotator)",
            "Primer3 (Automated PCR Primer Designer)",
            "BioSeq Analysis (GC-Skew & Sequence Complexity)",
            "Pyfaidx Index (High-Speed FASTA Substring)",
            "PyCircos Plot (Circular Genomic Visualizer)",
            "Biopython Phylo (Phylogenetic Tree Parser)",
            "Edlib Align (Ultra-Fast Levenshtein Match)",
            "CodonW Analysis (Codon Usage Bias Metrics)",
            "PyTFBS Scan (Transcription Factor Binding Sites)",
            "BioKEGG Pathway (Metabolic Pathway Mapping)",
            "Helical Wheel (Amphipathic Helices Plotter)",
            "SeqLogo Generator (Information Content Engine)",
            "Sequence Info (Format & Integrity Inspector)"
        ]
    },
    {
        "name": "Single-Cell & Spatial Omics",
        "nodes": [
            "Scanpy QC (Mitochondrial & Gene Count Filter)",
            "Scanpy Normalize (Total Counts & Log1p Transformation)",
            "Scanpy Cluster (Leiden & Louvain Community Detection)",
            "Scanpy Marker Genes (Differential Expression Ranking)",
            "ScRNA Visualization (UMAP & t-SNE Embedding Plotter)",
            "ScRNA Report (Automated Single-Cell QC Summary)",
            "AnnData I/O (In-Memory H5AD Matrix Bridge)",
            "scVI Tools (Deep Generative Latent Embedding)",
            "SCANVI (Semi-Supervised Cell-Type Annotation)",
            "scVelo Dynamics (RNA Velocity Kinetics Engine)",
            "CellRank Fate (Lineage Trajectory & Fate Mapping)",
            "Squidpy Spatial (Spatial Transcriptomics Graphs)",
            "Tangram Spatial (Single-Cell Spatial Alignment)",
            "Cell2location (Cell-Type Deconvolution Engine)",
            "PySCENIC (Gene Regulatory Network Regulon)",
            "CellPhoneDB (Ligand-Receptor Cell Crosstalk)",
            "GSEApy Enrichment (Single-Cell Pathway Enrichment)",
            "Muon Multimodal (CITE-seq & Multiome Integrator)",
            "CITE-seq-Count (Antibody-Derived Tag Quantifier)"
        ]
    },
    {
        "name": "Genome Assembly & Long-Read",
        "nodes": [
            "Assembly Input Validator (FASTQ/FASTA Integrity)",
            "Assembly Fastp Trim (Paired-End Read Trimmer)",
            "SPAdes Assemble (Multi k-mer De Bruijn Assembler)",
            "QUAST QC (N50 & Assembly Metric Evaluator)",
            "Assembly Visualization (Bandage Graph Renderer)",
            "Assembly Report (Automated Assembly QC Summary)",
            "Pysam Analysis (SAM/BAM/CRAM Alignment Engine)",
            "cyvcf2 Variant (High-Speed VCF Parser & Filter)",
            "pybedtools Interval (Genomic Interval Operations)",
            "Mappy Align (Minimap2 Python Binding Engine)",
            "pyfastx Index (High-Speed FASTA/FASTQ Accessor)",
            "SeqKit Tool (Multi-Threaded Sequence Toolkit)",
            "Bowtie2 Align (Gap-Aware Short-Read Aligner)",
            "Mosdepth Coverage (High-Speed BAM Depth Profiler)",
            "Sniffles2 SV (Long-Read Structural Variant Caller)",
            "CuteSV (High-Precision Breakpoint SV Caller)",
            "Flye Assemble (Repeat-Graph Long-Read Assembler)",
            "Hifiasm Assemble (PacBio HiFi De Novo Assembler)",
            "Racon Polish (Consensus Polishing Module)",
            "Medaka Consensus (Oxford Nanopore Polish Engine)",
            "DeepVariant Call (CNN Germline Variant Caller)"
        ]
    },
    {
        "name": "Epigenomics & Functional Screening",
        "nodes": [
            "ATAC Input Validator (Paired FASTQ Check)",
            "ATAC Fastp Trim (Poly-G & Adapter Trimmer)",
            "ATAC BWA-MEM2 Index (Reference Genome Indexer)",
            "ATAC BWA-MEM2 Align (Paired-End Alignment Engine)",
            "ATAC MarkDuplicates (Picard Duplicate Tagging)",
            "ATAC Quality Filter (Mitochondrial & Low-MAPQ Filter)",
            "MACS3 Peak Calling (BAMPE Open Chromatin Caller)",
            "ATAC Peak Visualization (Coverage Profile Tracks)",
            "ATAC Report (Automated Chromatin QC Summary)",
            "deepTools Profile (Metagene & Heatmap Profiler)",
            "TOBIAS Footprint (ATAC Footprinting & TF Occupancy)",
            "PyGenrich (ATAC-seq Peak Calling Engine)",
            "CRISPResso2 (HDR/NHEJ Genome Editing Analyzer)",
            "MAGeCK Screen (CRISPR Screen Guide Enrichment)",
            "Scikit-Fusion (Multi-Omics Relational Data Fusion)",
            "SEACR Peak (CUT&Tag / CUT&RUN Peak Caller)",
            "HOMER Motif (ChIP/ATAC Motif Discovery Engine)",
            "MEME Suite (De Novo Sequence Motif Discovery)",
            "MethylDackel (Bisulfite Methylation Extractor)",
            "Cooler Matrix (Hi-C Genomic Contact Matrix)",
            "Cooltools TAD (Insulation & Domain Boundary Caller)",
            "Chromosight Loop (Chromatin Loop Detection Engine)"
        ]
    },
    {
        "name": "Metagenome & Microbiome Virome",
        "nodes": [
            "Metagenome Validator (FASTQ Integrity Inspector)",
            "Metagenome Fastp Trim (Adapter & Quality Filter)",
            "Kraken2 Classify (Exact k-mer Taxonomic Classifier)",
            "Bracken Abundance (Bayesian Abundance Re-Estimator)",
            "Metagenome Visualization (Krona & Stacked Bar)",
            "Metagenome Report (Taxonomic Profiling Summary)",
            "scikit-bio Diversity (Alpha & Beta Diversity Metrics)",
            "BIOM Format (Sparse Contingency Table Bridge)",
            "ETE3 Tree Parser (Phylogenetic Tree Engine)",
            "FastUniFrac (Phylogenetic Distance Matrix)",
            "DADA2 Amplicon (Error-Corrected ASV Denoiser)",
            "HUMAnN3 Pathway (Functional Metabolic Profiler)",
            "MetaPhlAn4 Profile (Clade-Specific Marker Profiler)",
            "geNomad Virome (Virus & Plasmid Identification)",
            "VirSorter2 (Mining Viral Sequences Engine)",
            "CheckV Quality (Viral Genome Completeness Evaluator)",
            "Prokka Annotation (Rapid Prokaryotic Gene Annotator)",
            "Bakta Annotation (Modern Bacterial Genome Annotator)",
            "AMRFinderPlus (Antimicrobial Resistance Scanner)",
            "Nextstrain Augur (Phylodynamic Tracking Engine)"
        ]
    },
    {
        "name": "Bulk RNA-Seq & Core Pipeline",
        "nodes": [
            "Sample Metadata Validator (Samplesheet Inspector)",
            "FastQC Node (Base Quality Metrics Engine)",
            "Fastp QC Node (Automated Pre-Check Reporter)",
            "Fastp Trim Node (Poly-G & Adapter Read Trimmer)",
            "Trimmomatic Trim (Sliding Window Quality Trimmer)",
            "Salmon Index (Decoy-Aware Transcriptome Indexer)",
            "Salmon Quant (Quasi-Mapping Abundance Quantifier)",
            "Tximport Node (Gene-Level Count Matrix Summarizer)",
            "DESeq2 Analysis (Negative Binomial Generalized Linear)",
            "DESeq2 Visualization (Volcano, PCA & Heatmap Plots)",
            "ComfyBIO Report (Markdown Audit & Provenance Logger)",
            "10x Cell Ranger Count (Droplet Matrix Processor)",
            "Kallisto Quant (Pseudoalignment Abundance Quantifier)",
            "Alevin-Fry Quant (High-Speed Single-Cell Engine)",
            "STARsolo Quant (Ultra-Fast Single-Cell RNA-Seq)",
            "StringTie2 Assemble (Transcriptome Assembly Engine)",
            "FLAIR Isoform (Full-Length Nanopore Isoform Caller)",
            "IsoTools Splicing (Alternative Splicing Event Analyzer)",
            "edgeR Analysis (Empirical Bayes Dispersion DEG)",
            "limma-voom Analysis (Precision-Weight Linear Modeling)"
        ]
    },
    {
        "name": "Publication Quality Visualizers",
        "nodes": [
            "Volcano Plot Visualizer (DEG Statistical Significance)",
            "Manhattan Plot Visualizer (GWAS Genome-Wide Association)",
            "UMAP Scatter Visualizer (High-D Cell Embeddings)",
            "Clustermap Heatmap Visualizer (Hierarchical Clustering)",
            "GSEA Enrichment Plot (Running Enrichment Score Curves)",
            "OncoPrint Visualizer (Genomic Alteration Heatmaps)",
            "Q-Q Plot Visualizer (Quantile-Quantile Calibration)",
            "Sankey Cell Fate (Cell Lineage Transition Flux)",
            "Microbiome Stacked Bar (Relative Taxonomic Abundance)",
            "PCoA Scatter Visualizer (UniFrac Distance Projection)",
            "Synteny Genome Visualizer (Synteny & Orthology Blocks)",
            "ChIP/ATAC Coverage Profile (Signal Track Heatmaps)",
            "MA Plot Visualizer (Log Fold-Change vs Mean Expression)",
            "Spatial Tissue Overlay (H&E Histology Image Alignment)",
            "MD Trajectory Plotter (RMSD & RMSF Dynamic Curves)",
            "Ramachandran Plot Visualizer (Peptide Dihedral Angles)",
            "Protein-Ligand Interaction (2D/3D Contact Profiler)",
            "Phylogenetic Tree Visualizer (Circular & Linear Trees)",
            "Linkage Disequilibrium Visualizer (LD Heatmap Matrix)",
            "Kaplan-Meier Survival Plot (Log-Rank Survival Curves)"
        ]
    },
    {
        "name": "CADD & Structural Biology",
        "nodes": [
            "ColabFold AlphaFold (MMseqs2-Accelerated Structure)",
            "ESMFold Node (Large Protein Language Model)",
            "DiffDock Predict (Generative Diffusion Blind Docking)",
            "GNINA Docking (Deep Learning Molecular Docking)",
            "RDKit Cheminformatics (SMILES & Fingerprint Descriptors)",
            "OpenMM Simulation (GPU Molecular Dynamics Simulation)",
            "MDAnalysis Node (Trajectory Topology & Coordinates)",
            "MDTraj Node (High-Throughput Trajectory Analyzer)",
            "TorchDrug Node (Drug Discovery Graph Deep Learning)",
            "OpenFold Node (Trainable PyTorch-Native AlphaFold2)",
            "ProDy Dynamics (Elastic Network & Normal Mode Analysis)",
            "AutoDock Vina (Physics-Based Molecular Docking)",
            "P2Rank Pocket (Machine Learning Binding Site Predictor)",
            "FPocket Node (Voronoi Tessellation Pocket Detector)",
            "PLUMED Node (Metadynamics Enhanced Sampling Engine)",
            "pmx FEP (Alchemical Free Energy Perturbation)",
            "OpenBabel Convert (Multi-Format Chemical Converter)",
            "Smina Docking (Custom Scoring Function Optimizer)"
        ]
    },
    {
        "name": "Proteomics & Metabolomics",
        "nodes": [
            "Pyteomics MS (Peptide & Mass Spectrometry Parser)",
            "PyOpenMS Feature (LC-MS 2D/3D Feature Detection Engine)",
            "Matchms Spectrum (Mass Spectral Similarity & Filtering)",
            "Spec2Vec Embedding (NLP-Driven Word Spectral Matcher)",
            "MassQL Query (Mass Spectrometry Pattern Query Engine)",
            "msdeisotope (High-Resolution Charge Deconvolution)",
            "DIA-NN Quant (Neural Network Data-Independent Proteomics)",
            "MSFragger Search (Ultra-Fast Peptide Database Search)",
            "msconvert Convert (Vendor RAW to Open mzML Converter)",
            "MS-DIAL Lipid (Lipidomics Identification & Annotation)",
            "SIRIUS Structure (Metabolite Formula & De Novo MS/MS)",
            "MaxQuant Quant (Label-Free Proteomics Quantification)",
            "Perseus CLI (Statistical Omics Analysis Workbench)",
            "MetaboAnalystR (Functional Pathway Enrichment Analysis)",
            "MHCquant Neoantigen (Immunopeptidomics Target Pipeline)"
        ]
    },
    {
        "name": "DNA Variant Calling",
        "nodes": [
            "Variant Input Validator (FASTQ & Reference Validator)",
            "BWA-MEM2 Index Node (Burrows-Wheeler Genome Indexer)",
            "BWA-MEM2 Align Node (Paired-End DNA Sequence Aligner)",
            "MarkDuplicates Node (Picard Optical Duplicate Tagger)",
            "BCFtools Call Node (Multiallelic Genomic Variant Caller)",
            "BCFtools Filter Node (Depth, QUAL & Strand Bias Filter)",
            "Variant Visualization Node (SNV/Indel Spectrum & Ti/Tv)",
            "Variant Report Node (Automated Variant QC Summary)"
        ]
    }
]

def render_vertical_tree(out_path='figures/fig2_node_domain_tree.png'):
    fig_w, fig_h = 19.0, 28.0
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    domain_gap = 1.8
    slot_allocations = []
    current_y = 0.0

    for d in reversed(DOMAIN_DATA):
        n_leaves = len(d["nodes"])
        leaf_slots = [current_y + i * 1.0 for i in range(n_leaves)]
        domain_center_y = np.mean(leaf_slots)
        slot_allocations.append({
            "name": d["name"],
            "count": n_leaves,
            "nodes": d["nodes"],
            "leaf_ys_raw": leaf_slots,
            "domain_y_raw": domain_center_y
        })
        current_y += n_leaves * 1.0 + domain_gap

    slot_allocations = list(reversed(slot_allocations))
    max_raw_y = current_y - domain_gap

    y_min_target, y_max_target = 0.025, 0.975
    def norm_y(val):
        return y_min_target + (val / max_raw_y) * (y_max_target - y_min_target)

    for d in slot_allocations:
        d["leaf_ys"] = [norm_y(y) for y in d["leaf_ys_raw"]]
        d["domain_y"] = norm_y(d["domain_y_raw"])

    x_root = 0.14
    x_domain = 0.44
    x_leaf = 0.63
    y_root = norm_y(max_raw_y / 2.0)

    curve_color_root = '#94a3b8'   # Slate 400
    curve_color_leaf = '#cbd5e1'   # Slate 300
    circle_border = '#2563eb'      # Blue 600
    circle_fill = '#e0f2fe'        # Sky 100
    root_border = '#1e40af'        # Blue 800
    root_fill = '#bfdbfe'          # Blue 200

    def draw_smooth_spline(x1, y1, x2, y2, color, lw, alpha=0.9):
        dx = x2 - x1
        c1 = (x1 + dx * 0.48, y1)
        c2 = (x2 - dx * 0.48, y2)
        verts = [(x1, y1), c1, c2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        patch = patches.PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # 1. Connect Root to Domains
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain, d["domain_y"], color=curve_color_root, lw=1.35, alpha=0.9)

    # 2. Connect Domains to Leaves
    for d in slot_allocations:
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_domain, d["domain_y"], x_leaf, y_l, color=curve_color_leaf, lw=0.9, alpha=0.85)

    # 3. Draw Root Node
    ax.scatter([x_root], [y_root], s=170, facecolor=root_fill, edgecolor=root_border, linewidth=2.0, zorder=5)
    ax.text(x_root - 0.015, y_root, "ComfyBIOWMS Nodes", 
            ha='right', va='center', fontsize=13, fontweight='bold', color='#0f172a', zorder=6)

    # 4. Draw Domain Nodes & Labels
    for d in slot_allocations:
        ax.scatter([x_domain], [d["domain_y"]], s=105, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.75, zorder=5)
        cat_text = f"{d['name']} ({d['count']})"
        ax.text(x_domain - 0.012, d["domain_y"], cat_text, 
                ha='right', va='center', fontsize=10.0, fontweight='bold', color='#1e293b', zorder=6)

    # 5. Draw Leaf Nodes & Labels
    for d in slot_allocations:
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf], [y_l], s=40, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.2, zorder=5)
            ax.text(x_leaf + 0.008, y_l, node_label, 
                    ha='left', va='center', fontsize=7.2, color='#334155', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved vertical tree to {out_path}")

def render_bilateral_landscape_tree(out_path='figures/fig2_node_domain_tree_landscape.png'):
    # 5 domains on left, 5 domains on right, root in the exact center!
    # Creates a perfect landscape ratio (24 x 14 inches)
    fig_w, fig_h = 26.0, 15.0
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    # Left domains (5) and Right domains (5)
    left_domains_raw = DOMAIN_DATA[:5]
    right_domains_raw = DOMAIN_DATA[5:]

    def process_wing(domains, direction="right"):
        domain_gap = 1.4
        allocations = []
        current_y = 0.0
        for d in reversed(domains):
            n_leaves = len(d["nodes"])
            leaf_slots = [current_y + i * 1.0 for i in range(n_leaves)]
            domain_center_y = np.mean(leaf_slots)
            allocations.append({
                "name": d["name"],
                "count": n_leaves,
                "nodes": d["nodes"],
                "leaf_ys_raw": leaf_slots,
                "domain_y_raw": domain_center_y
            })
            current_y += n_leaves * 1.0 + domain_gap
        allocations = list(reversed(allocations))
        max_raw_y = current_y - domain_gap
        return allocations, max_raw_y

    left_alloc, max_l = process_wing(left_domains_raw, "left")
    right_alloc, max_r = process_wing(right_domains_raw, "right")

    y_min_target, y_max_target = 0.04, 0.96
    for d in left_alloc:
        d["leaf_ys"] = [y_min_target + (y / max_l) * (y_max_target - y_min_target) for y in d["leaf_ys_raw"]]
        d["domain_y"] = y_min_target + (d["domain_y_raw"] / max_l) * (y_max_target - y_min_target)

    for d in right_alloc:
        d["leaf_ys"] = [y_min_target + (y / max_r) * (y_max_target - y_min_target) for y in d["leaf_ys_raw"]]
        d["domain_y"] = y_min_target + (d["domain_y_raw"] / max_r) * (y_max_target - y_min_target)

    x_root = 0.50
    y_root = 0.50

    # Left coordinates
    x_domain_l = 0.33
    x_leaf_l = 0.20

    # Right coordinates
    x_domain_r = 0.67
    x_leaf_r = 0.80

    curve_root_color = '#94a3b8'
    curve_leaf_color = '#cbd5e1'
    circle_border = '#2563eb'
    circle_fill = '#e0f2fe'
    root_border = '#1e40af'
    root_fill = '#bfdbfe'

    def draw_spline(x1, y1, x2, y2, color, lw, alpha=0.9):
        dx = x2 - x1
        c1 = (x1 + dx * 0.48, y1)
        c2 = (x2 - dx * 0.48, y2)
        verts = [(x1, y1), c1, c2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        patch = patches.PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # --- Draw Left Wing ---
    for d in left_alloc:
        draw_spline(x_root, y_root, x_domain_l, d["domain_y"], color=curve_root_color, lw=1.35)
        for y_l in d["leaf_ys"]:
            draw_spline(x_domain_l, d["domain_y"], x_leaf_l, y_l, color=curve_leaf_color, lw=0.9)

    for d in left_alloc:
        ax.scatter([x_domain_l], [d["domain_y"]], s=100, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.75, zorder=5)
        cat_text = f"({d['count']}) {d['name']}"
        ax.text(x_domain_l + 0.012, d["domain_y"], cat_text, 
                ha='left', va='center', fontsize=9.8, fontweight='bold', color='#1e293b', zorder=6)

        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf_l], [y_l], s=36, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.2, zorder=5)
            ax.text(x_leaf_l - 0.006, y_l, node_label, 
                    ha='right', va='center', fontsize=6.8, color='#334155', zorder=6)

    # --- Draw Right Wing ---
    for d in right_alloc:
        draw_spline(x_root, y_root, x_domain_r, d["domain_y"], color=curve_root_color, lw=1.35)
        for y_l in d["leaf_ys"]:
            draw_spline(x_domain_r, d["domain_y"], x_leaf_r, y_l, color=curve_leaf_color, lw=0.9)

    for d in right_alloc:
        ax.scatter([x_domain_r], [d["domain_y"]], s=100, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.75, zorder=5)
        cat_text = f"{d['name']} ({d['count']})"
        ax.text(x_domain_r - 0.012, d["domain_y"], cat_text, 
                ha='right', va='center', fontsize=9.8, fontweight='bold', color='#1e293b', zorder=6)

        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf_r], [y_l], s=36, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.2, zorder=5)
            ax.text(x_leaf_r + 0.006, y_l, node_label, 
                    ha='left', va='center', fontsize=6.8, color='#334155', zorder=6)

    # --- Draw Center Root Node ---
    ax.scatter([x_root], [y_root], s=180, facecolor=root_fill, edgecolor=root_border, linewidth=2.2, zorder=7)
    ax.text(x_root, y_root + 0.035, "ComfyBIOWMS Nodes", 
            ha='center', va='bottom', fontsize=12.5, fontweight='bold', color='#0f172a', zorder=8)
    ax.text(x_root, y_root - 0.035, "(189 Multi-Omics Tools)", 
            ha='center', va='top', fontsize=9.5, color='#475569', zorder=8)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved landscape bilateral tree to {out_path}")

if __name__ == '__main__':
    render_vertical_tree('figures/fig2_node_domain_tree.png')
    render_bilateral_landscape_tree('figures/fig2_node_domain_tree_landscape.png')
    # Also update fig2_node_domain_distribution.png to the vertical tree or as needed
    import shutil
    shutil.copyfile('figures/fig2_node_domain_tree.png', 'figures/fig2_node_domain_distribution.png')
    print("Updated figures/fig2_node_domain_distribution.png")
