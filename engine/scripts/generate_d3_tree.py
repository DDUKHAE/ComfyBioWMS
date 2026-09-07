import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np

def generate_d3_style_tree():
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    # 10 Domains and their comprehensive curated node toolchain
    domain_data = [
        {
            "name": "Biopython & Sequence Ops",
            "nodes": [
                "SeqIO multi-format sequence parser",
                "SeqTranscribe DNA/RNA/protein translation",
                "AlignIO multiple sequence alignment",
                "BLAST NCBI remote & local search",
                "Entrez NCBI E-utilities query",
                "BioPDB 3D macromolecular parser",
                "Restriction enzyme cleavage simulator",
                "BioMotifs PSSM & IUPAC motif search",
                "Logomaker visualizer & sequence logo",
                "ProtParam physicochemical calculator",
                "Biotite structural superimposition",
                "Squiggle raw nanopore signal plotter",
                "PyHMMER hidden Markov model search",
                "DNA Features Viewer genomic annotator",
                "Primer3 PCR primer automated designer",
                "BioSeq GC skew & complexity metrics",
                "Pyfaidx ultra-fast FASTA indexer",
                "PyCircos circular genomic visualizer",
                "Biopython Phylo dendrogram engine",
                "Edlib fast Levenshtein aligner",
                "CodonW codon usage bias analyzer",
                "PyTFBS transcription factor scanner",
                "BioKEGG metabolic pathway visualizer",
                "Helical Wheel amphipathic plotter",
                "SeqLogo information content engine",
                "Biopython Sequence Info format validator"
            ]
        },
        {
            "name": "Single-Cell & Spatial Omics",
            "nodes": [
                "Scanpy QC mitochondrial & gene filter",
                "Scanpy Normalize total counts & log1p",
                "Scanpy Cluster Leiden & Louvain graphs",
                "Scanpy Marker Genes differential ranking",
                "ScRNA Visualization UMAP / t-SNE plotter",
                "ScRNA Report automated workflow summary",
                "AnnData I/O in-memory H5AD bridge",
                "scVI deep generative latent embedding",
                "SCANVI semi-supervised cell annotator",
                "scVelo RNA velocity kinetics engine",
                "CellRank lineage trajectory fate mapper",
                "Squidpy spatial transcriptomics graphs",
                "Tangram spatial single-cell mapper",
                "Cell2location cell-type deconvolution",
                "PySCENIC gene regulatory network regulon",
                "CellPhoneDB ligand-receptor crosstalk",
                "GSEApy single-cell pathway enrichment",
                "Muon multimodal CITE-seq integrator",
                "CITE-seq-Count antibody tag quantifier"
            ]
        },
        {
            "name": "Genome Assembly & Long-Read",
            "nodes": [
                "Assembly Input Validator FASTQ/FASTA",
                "Assembly Fastp Trim paired-end trimmer",
                "SPAdes de novo multi k-mer assembler",
                "QUAST assembly quality evaluator",
                "Assembly Visualization Bandage graph",
                "Assembly Report automated QC summary",
                "Pysam SAM/BAM/CRAM alignment parser",
                "cyvcf2 high-speed VCF record parser",
                "pybedtools genomic interval operations",
                "Mappy Minimap2 Python alignment engine",
                "pyfastx indexed FASTQ/FASTA accessor",
                "SeqKit multi-threaded sequence toolkit",
                "Bowtie2 gap-aware short-read aligner",
                "Mosdepth high-speed BAM depth profiler",
                "Sniffles2 long-read structural caller",
                "CuteSV breakpoint structural caller",
                "Flye repeat-graph long-read assembler",
                "Hifiasm PacBio HiFi de novo assembler",
                "Racon consensus polishing module",
                "Medaka Oxford Nanopore polish engine",
                "DeepVariant CNN germline caller"
            ]
        },
        {
            "name": "Epigenomics & Functional Genomics",
            "nodes": [
                "ATAC Input Validator paired FASTQ check",
                "ATAC Fastp Trim poly-G & adapter trimmer",
                "ATAC BWA-MEM2 Reference genome indexer",
                "ATAC BWA-MEM2 Paired alignment engine",
                "ATAC MarkDuplicates Picard duplicate tag",
                "ATAC Quality Filter mitochondrial filter",
                "MACS3 BAMPE open chromatin peak caller",
                "ATAC Peak Visualization coverage track",
                "ATAC Report automated chromatin summary",
                "deepTools metagene & profile heatmaps",
                "TOBIAS transcription factor footprinting",
                "PyGenrich ATAC-seq peak caller",
                "CRISPResso2 HDR/NHEJ editing analyzer",
                "MAGeCK CRISPR screen enrichment",
                "Scikit-Fusion multi-omics data fusion",
                "SEACR CUT&Tag / CUT&RUN peak caller",
                "HOMER motif enrichment analyzer",
                "MEME Suite de novo motif discovery",
                "MethylDackel bisulfite methylation caller",
                "Cooler Hi-C contact matrix generator",
                "Cooltools TAD insulation & boundary caller",
                "Chromosight chromatin loop detection"
            ]
        },
        {
            "name": "Metagenome & Pathogen Virome",
            "nodes": [
                "Metagenome Validator raw read inspector",
                "Metagenome Fastp Trim quality filter",
                "Kraken2 exact k-mer taxonomic classifier",
                "Bracken Bayesian abundance re-estimator",
                "Metagenome Visualization Krona/stacked",
                "Metagenome Report automated taxon summary",
                "scikit-bio alpha & beta diversity metrics",
                "BIOM Format sparse contingency matrix",
                "ETE3 phylogenetic tree parser",
                "FastUniFrac phylogenetic distance matrix",
                "DADA2 error-corrected ASV denoiser",
                "HUMAnN3 functional pathway profiler",
                "MetaPhlAn4 clade-specific marker profiler",
                "geNomad virus & plasmid identification",
                "VirSorter2 viral sequence miner",
                "CheckV viral genome completeness evaluator",
                "Prokka rapid prokaryotic gene annotator",
                "Bakta modern bacterial genome annotator",
                "AMRFinderPlus antimicrobial resistance",
                "Nextstrain Augur phylodynamic tracker"
            ]
        },
        {
            "name": "Bulk RNA-Seq & Core Pipeline",
            "nodes": [
                "Sample Metadata Validator samplesheet check",
                "FastQC per-base quality metric engine",
                "Fastp QC automated read quality reporter",
                "Fastp Trim adapter & quality trimmer",
                "Trimmomatic sliding-window quality trimmer",
                "Salmon Index transcript decoy-aware indexer",
                "Salmon Quant quasi-mapping quantifier",
                "Tximport gene-level count summarizer",
                "DESeq2 Analysis negative binomial GLM",
                "DESeq2 Visualization volcano/PCA plots",
                "ComfyBIO Report execution audit logger",
                "10x Cell Ranger droplet matrix processor",
                "Kallisto pseudoalignment RNA quantifier",
                "Alevin-Fry fast single-cell quantification",
                "STARsolo ultra-fast single-cell aligner",
                "StringTie2 transcript assembly engine",
                "FLAIR full-length Nanopore isoform caller",
                "IsoTools alternative splicing event analyzer",
                "edgeR empirical Bayes dispersion DEG",
                "limma-voom precision-weight linear models"
            ]
        },
        {
            "name": "Publication Visualizers",
            "nodes": [
                "Volcano Plot visualizer (DEG significance)",
                "Manhattan Plot visualizer (GWAS p-values)",
                "UMAP Scatter visualizer (cell embeddings)",
                "Clustermap Heatmap visualizer (hierarchical)",
                "GSEA Enrichment Plot (running score curves)",
                "OncoPrint visualizer (genomic alteration)",
                "Q-Q Plot visualizer (uniform calibration)",
                "Sankey Cell Fate (lineage transition flux)",
                "Microbiome Stacked Bar (taxonomic composition)",
                "PCoA Scatter visualizer (UniFrac distance)",
                "Synteny Genome visualizer (orthology blocks)",
                "ChIP/ATAC Coverage Profile (signal tracks)",
                "MA Plot visualizer (fold-change vs intensity)",
                "Spatial Tissue Overlay (H&E alignment)",
                "MD Trajectory Plotter (RMSD / RMSF curves)",
                "Ramachandran Plot visualizer (dihedral angles)",
                "Protein-Ligand Interaction (2D/3D contacts)",
                "Phylogenetic Tree visualizer (cladograms)",
                "Linkage Disequilibrium visualizer (LD matrix)",
                "Kaplan-Meier Survival Plot (log-rank test)"
            ]
        },
        {
            "name": "CADD & Structural Biology",
            "nodes": [
                "ColabFold AlphaFold2 monomer & multimer",
                "ESMFold protein language structure model",
                "DiffDock generative blind docking engine",
                "GNINA deep learning molecular docking",
                "RDKit cheminformatics descriptor engine",
                "OpenMM GPU molecular dynamics simulation",
                "MDAnalysis trajectory topology parser",
                "MDTraj high-throughput coordinate analyzer",
                "TorchDrug drug discovery deep learning",
                "OpenFold PyTorch-native AlphaFold2",
                "ProDy elastic network normal mode model",
                "AutoDock Vina physics-based docking",
                "P2Rank machine learning pocket predictor",
                "FPocket Voronoi geometry pocket detector",
                "PLUMED metadynamics enhanced sampling",
                "pmx alchemical free energy perturbation",
                "OpenBabel chemical format converter",
                "Smina custom scoring docking engine"
            ]
        },
        {
            "name": "Proteomics & Metabolomics",
            "nodes": [
                "Pyteomics peptide & MS data parser",
                "PyOpenMS LC-MS 2D/3D feature detection",
                "Matchms mass spectral similarity engine",
                "Spec2Vec word-embedding spectra matcher",
                "MassQL mass spectrometry query language",
                "msdeisotope high-res charge deconvolution",
                "DIA-NN neural network DIA proteomics",
                "MSFragger ultra-fast peptide database search",
                "msconvert RAW to mzML format converter",
                "MS-DIAL lipidomics & untargeted profiling",
                "SIRIUS metabolite structure & formula MS/MS",
                "MaxQuant label-free quantification engine",
                "Perseus CLI statistical omics workbench",
                "MetaboAnalystR functional pathway enrichment",
                "MHCquant neoantigen immunopeptidomics"
            ]
        },
        {
            "name": "DNA Variant Calling",
            "nodes": [
                "Variant Input Validator FASTQ/FASTA check",
                "BWA-MEM2 Index Burrows-Wheeler indexer",
                "BWA-MEM2 Align paired-end DNA aligner",
                "MarkDuplicates Picard optical duplicate tag",
                "BCFtools Call multiallelic variant caller",
                "BCFtools Filter depth & quality filter",
                "Variant Visualization SNV/Indel spectrum",
                "Variant Report Ti/Tv ratio & QC summary"
            ]
        }
    ]

    total_leaves = sum(len(d["nodes"]) for d in domain_data)
    print(f"Generating D3-style tree with {len(domain_data)} domains and {total_leaves} leaf nodes...")

    # Calculate optimal canvas dimensions
    # Using 28 inches height and 18 inches width for generous spacing
    fig_w, fig_h = 18.0, 26.0
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    # Domain spacing: assign vertical slots with padding between domains
    # Leaf heights calculated with intra-domain spacing and inter-domain gaps
    domain_gap = 1.6  # gap multiplier between domain blocks
    slot_allocations = []
    
    current_y = 0.0
    for d in reversed(domain_data):  # bottom to top
        n_leaves = len(d["nodes"])
        # each leaf gets 1.0 unit
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

    # Reverse back so top domain is at the top
    slot_allocations = list(reversed(slot_allocations))
    max_raw_y = current_y - domain_gap

    # Normalize to Y range [0.03, 0.97]
    y_min_target, y_max_target = 0.03, 0.97
    def norm_y(val):
        return y_min_target + (val / max_raw_y) * (y_max_target - y_min_target)

    for d in slot_allocations:
        d["leaf_ys"] = [norm_y(y) for y in d["leaf_ys_raw"]]
        d["domain_y"] = norm_y(d["domain_y_raw"])

    # X coordinates
    x_root = 0.13
    x_domain = 0.44
    x_leaf = 0.65

    y_root = norm_y(max_raw_y / 2.0)

    # Styling constants matching the reference image
    # Soft slate gray curves
    curve_color_root = '#94a3b8'   # Slate 400
    curve_color_leaf = '#cbd5e1'   # Slate 300
    
    # Blue circle nodes
    circle_border = '#3b82f6'      # Blue 500
    circle_fill = '#e0f2fe'        # Sky 100 / light blue
    root_border = '#1d4ed8'        # Blue 700
    root_fill = '#bfdbfe'          # Blue 200

    def draw_smooth_spline(x1, y1, x2, y2, color, lw, alpha=0.9):
        # Perfect horizontal tangent cubic bezier (D3 tree style)
        dx = x2 - x1
        c1 = (x1 + dx * 0.48, y1)
        c2 = (x2 - dx * 0.48, y2)
        verts = [(x1, y1), c1, c2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        patch = patches.PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # 1. Connect Root to Domains
    for d in slot_allocations:
        draw_smooth_spline(x_root, y_root, x_domain, d["domain_y"], color=curve_color_root, lw=1.3, alpha=0.9)

    # 2. Connect Domains to Leaves
    for d in slot_allocations:
        for y_l in d["leaf_ys"]:
            draw_smooth_spline(x_domain, d["domain_y"], x_leaf, y_l, color=curve_color_leaf, lw=0.9, alpha=0.85)

    # 3. Draw Root Node
    ax.scatter([x_root], [y_root], s=160, facecolor=root_fill, edgecolor=root_border, linewidth=2.0, zorder=5)
    ax.text(x_root - 0.015, y_root, "ComfyBIOWMS Nodes", 
            ha='right', va='center', fontsize=12.5, fontweight='bold', color='#0f172a', zorder=6)

    # 4. Draw Domain Nodes & Labels
    for d in slot_allocations:
        ax.scatter([x_domain], [d["domain_y"]], s=100, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.7, zorder=5)
        # Category label to left of circle
        cat_text = f"{d['name']} ({d['count']})"
        ax.text(x_domain - 0.012, d["domain_y"], cat_text, 
                ha='right', va='center', fontsize=9.8, fontweight='bold', color='#1e293b', zorder=6)

    # 5. Draw Leaf Nodes & Labels
    for d in slot_allocations:
        for y_l, node_label in zip(d["leaf_ys"], d["nodes"]):
            ax.scatter([x_leaf], [y_l], s=38, facecolor=circle_fill, edgecolor=circle_border, linewidth=1.2, zorder=5)
            # Leaf label to right of circle
            ax.text(x_leaf + 0.008, y_l, node_label, 
                    ha='left', va='center', fontsize=6.9, color='#334155', zorder=6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

    out_png = 'figures/fig2_node_domain_tree_d3.png'
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight', pad_inches=0.15, facecolor='#ffffff')
    plt.close()
    print(f"Saved D3-style tree figure to {out_png}")

if __name__ == '__main__':
    generate_d3_style_tree()
