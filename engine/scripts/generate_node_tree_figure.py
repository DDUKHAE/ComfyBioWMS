import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np

def create_tree_diagram():
    # Set font family
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
    plt.rcParams['svg.fonttype'] = 'none'

    # Domain structure and nodes
    domain_data = [
        {
            "name": "Biopython & Sequence Ops",
            "count": 26,
            "nodes": [
                "Biopython SeqIO (Multi-Format Parser)",
                "Biopython SeqTranscribe (DNA/RNA/Protein)",
                "Biopython AlignIO (Multiple Alignment)",
                "Biopython BLAST (NCBI Remote/Local)",
                "Biopython Entrez (E-utilities NCBI)",
                "Biopython BioPDB (3D Structure Parser)",
                "Biopython Restriction (Enzyme Digestion)",
                "Biopython BioMotifs (PSSM & IUPAC)",
                "Logomaker Visualizer (Sequence Logo)",
                "Biopython ProtParam (Physicochemical)",
                "Biotite Structure Align (Superimposition)",
                "Squiggle Waveform (Nanopore Raw Signal)",
                "PyHMMER Search (HMM Profile Search)",
                "DNA Features Viewer (Genomic Annotation)",
                "Primer3 Design (PCR Primer Designer)",
                "BioSeq Analysis (GC/Skew/Complexity)",
                "Pyfaidx Index (Fast FASTA Substring)",
                "PyCircos Plot (Circos Genomic View)",
                "Biopython Phylo (Phylogenetic Trees)",
                "Edlib Align (Fast Levenshtein Match)",
                "CodonW Analysis (Codon Usage Bias)",
                "PyTFBS Scan (Transcription Factors)",
                "BioKEGG Pathway (Metabolic Mapping)",
                "Helical Wheel (Amphipathic Helices)",
                "SeqLogo Generator (Information Content)",
                "Biopython Sequence Info (Format Inspector)"
            ]
        },
        {
            "name": "Single-Cell & Spatial Omics",
            "count": 19,
            "nodes": [
                "Scanpy QC (Mitochondrial & Gene Filter)",
                "Scanpy Normalize (Total Count & Log1p)",
                "Scanpy Cluster (Leiden & Louvain)",
                "Scanpy Marker Genes (Wilcoxon / t-test)",
                "ScRNA Visualization (UMAP / t-SNE Plot)",
                "ScRNA Report (Automated scRNA Summary)",
                "AnnData I/O (In-Memory H5AD Matrix)",
                "scVI Tools (Deep Latent Representation)",
                "SCANVI (Semi-Supervised Cell Annotation)",
                "scVelo Dynamics (RNA Velocity Kinetics)",
                "CellRank Fate (Lineage Trajectory)",
                "Squidpy Spatial (Spatial Graph Analysis)",
                "Tangram Spatial (Single-Cell Mapping)",
                "Cell2location (Cell-Type Deconvolution)",
                "PySCENIC (Gene Regulatory Networks)",
                "CellPhoneDB (Ligand-Receptor Crosstalk)",
                "GSEApy Enrichment (Single-Cell Pathways)",
                "Muon Multimodal (CITE-seq / Multiome)",
                "CITE-seq-Count (Antibody-Derived Tags)"
            ]
        },
        {
            "name": "Genome Assembly & Long-Read",
            "count": 21,
            "nodes": [
                "Assembly Input Validator (FASTA/FASTQ)",
                "Assembly Fastp Trim (Paired-End Quality)",
                "SPAdes Assemble (Multi k-mer De Bruijn)",
                "QUAST QC (N50 / Contig Metrics)",
                "Assembly Visualization (Bandage Graph)",
                "Assembly Report (Automated Summary)",
                "Pysam Analysis (BAM/CRAM Operations)",
                "cyvcf2 Variant (High-Speed VCF Parser)",
                "pybedtools Interval (BED Overlaps)",
                "Mappy Align (Minimap2 Python Binding)",
                "pyfastx Index (Ultra-Fast Fasta/Fastq)",
                "SeqKit Tool (Multi-Threaded Sequence CLI)",
                "Bowtie2 Align (Fast Gap-Aware Aligner)",
                "Mosdepth Coverage (Fast BAM Depth)",
                "Sniffles2 SV (Long-Read Structural Variants)",
                "CuteSV (Accurate Long-Read SVs)",
                "Flye Assemble (Repeat-Graph Long-Read)",
                "Hifiasm Assemble (PacBio HiFi Assembler)",
                "Racon Polish (Consensus Module Polish)",
                "Medaka Consensus (Oxford Nanopore Polish)",
                "DeepVariant Call (CNN Germline Caller)"
            ]
        },
        {
            "name": "Epigenomics & Functional Genomics",
            "count": 22,
            "nodes": [
                "ATAC Input Validator (Paired FASTQ Check)",
                "ATAC Fastp Trim (Adapter Trimming)",
                "ATAC BWA-MEM2 Index (Reference Genome)",
                "ATAC BWA-MEM2 Align (Paired-End Alignment)",
                "ATAC MarkDuplicates (Picard Duplicate Tag)",
                "ATAC Quality Filter (Mitochondrial Filter)",
                "MACS3 Peak Calling (BAMPE Open Chromatin)",
                "ATAC Peak Visualization (Coverage Profile)",
                "ATAC Report (QC & Peak Summary)",
                "deepTools Profile (Heatmap & Metagene)",
                "TOBIAS Footprint (ATAC Footprinting)",
                "PyGenrich (ATAC-seq Peak Calling)",
                "CRISPResso2 (CRISPR Editing Analysis)",
                "MAGeCK Screen (CRISPR Screen Analysis)",
                "Scikit-Fusion (Multi-Omics Data Fusion)",
                "SEACR Peak (Cleave & CUT&Tag Peaks)",
                "HOMER Motif (ChIP-seq / ATAC Motif)",
                "MEME Suite (De Novo Motif Discovery)",
                "MethylDackel (Bisulfite-Seq Methylation)",
                "Cooler Matrix (Hi-C Genomic Contacts)",
                "Cooltools TAD (Insulation & Boundaries)",
                "Chromosight Loop (Chromatin Loop Caller)"
            ]
        },
        {
            "name": "Metagenome & Pathogen Virome",
            "count": 20,
            "nodes": [
                "Metagenome Validator (FASTQ Format Check)",
                "Metagenome Fastp Trim (Quality Filter)",
                "Kraken2 Classify (k-mer Taxon Classifier)",
                "Bracken Abundance (Bayesian Re-Estimation)",
                "Metagenome Visualization (Krona / Stacked)",
                "Metagenome Report (Taxonomic Summary)",
                "scikit-bio Diversity (Alpha & Beta Metrics)",
                "BIOM Format (Sparse Contingency Table)",
                "ETE3 Tree Parser (Phylogenetic Trees)",
                "FastUniFrac (Phylogenetic Distance)",
                "DADA2 Amplicon (Error-Corrected ASVs)",
                "HUMAnN3 Pathway (Functional Profiling)",
                "MetaPhlAn4 Profile (Clade-Specific Markers)",
                "geNomad Virome (Virus & Plasmid Identification)",
                "VirSorter2 (Mining Viral Sequences)",
                "CheckV Quality (Viral Genome Completeness)",
                "Prokka Annotation (Prokaryotic Annotation)",
                "Bakta Annotation (Fast Bacterial Annotation)",
                "AMRFinderPlus (Antimicrobial Resistance)",
                "Nextstrain Augur (Phylodynamic Tracking)"
            ]
        },
        {
            "name": "Bulk RNA-Seq & Core Pipeline",
            "count": 20,
            "nodes": [
                "Sample Metadata Validator (Samplesheet Check)",
                "FastQC Node (Base Quality Metrics)",
                "Fastp QC Node (Automated Pre-Check)",
                "Fastp Trim Node (Poly-G & Adapter Trim)",
                "Trimmomatic Trim (Sliding Window Quality)",
                "Salmon Index (Transcriptome Reference)",
                "Salmon Quant (Decoy-Aware Quasi-Mapping)",
                "Tximport Node (Gene-Level Summarization)",
                "DESeq2 Analysis (Negative Binomial GLM)",
                "DESeq2 Visualization (Volcano/PCA/Heatmap)",
                "ComfyBIO Report (Markdown Execution Audit)",
                "10x Cell Ranger Count (Single-Cell Droplet)",
                "Kallisto Quant (Pseudoalignment Quantification)",
                "Alevin-Fry Quant (Fast Single-Cell Engine)",
                "STARsolo Quant (High-Speed Single-Cell STAR)",
                "StringTie2 Assemble (Transcript Assembly)",
                "FLAIR Isoform (Nanopore Full-Length Isoforms)",
                "IsoTools Splicing (Alternative Splicing Events)",
                "edgeR Analysis (Exact Negative Binomial Test)",
                "limma-voom Analysis (Precision Weight Linear)"
            ]
        },
        {
            "name": "Publication Visualizers",
            "count": 20,
            "nodes": [
                "Volcano Plot Visualizer (DEG Significance)",
                "Manhattan Plot Visualizer (GWAS P-Values)",
                "UMAP Scatter Visualizer (Cell Embeddings)",
                "Clustermap Heatmap Visualizer (Expression)",
                "GSEA Enrichment Plot (Running ES Curves)",
                "OncoPrint Visualizer (Genomic Alterations)",
                "Q-Q Plot Visualizer (P-Value Calibration)",
                "Sankey Cell Fate (Lineage Transition)",
                "Microbiome Stacked Bar (Taxonomic Abundance)",
                "PCoA Scatter Visualizer (UniFrac Distances)",
                "Synteny Genome Visualizer (Synteny Blocks)",
                "ChIP/ATAC Coverage Profile (Signal Tracks)",
                "MA Plot Visualizer (Fold-Change vs Mean)",
                "Spatial Tissue Overlay (H&E Histology)",
                "MD Trajectory Plotter (RMSD / RMSF Curves)",
                "Ramachandran Plot Visualizer (Phi/Psi Angles)",
                "Protein-Ligand Interaction (2D/3D Contacts)",
                "Phylogenetic Tree Visualizer (Dendrograms)",
                "Linkage Disequilibrium Visualizer (LD Matrix)",
                "Kaplan-Meier Survival Plot (Survival Curves)"
            ]
        },
        {
            "name": "CADD & Structural Biology",
            "count": 18,
            "nodes": [
                "ColabFold AlphaFold (MMseqs2 Accelerated)",
                "ESMFold Node (Large Protein Language Model)",
                "DiffDock Predict (Generative Docking Model)",
                "GNINA Docking (Deep Learning Scoring)",
                "RDKit Cheminformatics (SMILES & Fingerprints)",
                "OpenMM Simulation (GPU Molecular Dynamics)",
                "MDAnalysis Node (Trajectory Topology Parser)",
                "MDTraj Node (High-Throughput Trajectory)",
                "TorchDrug Node (Drug Discovery Deep Learning)",
                "OpenFold Node (Trainable AlphaFold2)",
                "ProDy Dynamics (Elastic Network Models)",
                "AutoDock Vina (Classical Physics Docking)",
                "P2Rank Pocket (Ligand Binding Site ML)",
                "FPocket Node (Geometry Pocket Detection)",
                "PLUMED Node (Metadynamics Enhanced Sampling)",
                "pmx FEP (Free Energy Perturbation)",
                "OpenBabel Convert (Chemical Format Converter)",
                "Smina Docking (Scoring Function Optimizer)"
            ]
        },
        {
            "name": "Proteomics & Metabolomics",
            "count": 15,
            "nodes": [
                "Pyteomics MS (Peptide & Proteomics Parser)",
                "PyOpenMS Feature (LC-MS 2D/3D Peak Picking)",
                "Matchms Spectrum (Mass Spectral Comparison)",
                "Spec2Vec Embedding (NLP Spectral Similarity)",
                "MassQL Query (Mass Spectrometry Query)",
                "msdeisotope (High-Res Charge Deconvolution)",
                "DIA-NN Quant (Neural Network DIA Proteomics)",
                "MSFragger Search (Ultrafast Peptide Search)",
                "msconvert Convert (RAW to mzML Converter)",
                "MS-DIAL Lipid (Lipidomics Identification)",
                "SIRIUS Structure (Metabolite Formula MS/MS)",
                "MaxQuant Quant (Label-Free Quantification)",
                "Perseus CLI (Statistical Omics Workbench)",
                "MetaboAnalystR (Metabolomic Enrichment)",
                "MHCquant Neoantigen (Immunopeptidomics)"
            ]
        },
        {
            "name": "DNA Variant Calling",
            "count": 8,
            "nodes": [
                "Variant Input Validator (Fastq & Ref Validation)",
                "BWA-MEM2 Index Node (Burrows-Wheeler Index)",
                "BWA-MEM2 Align Node (Paired-End DNA Alignment)",
                "MarkDuplicates Node (Picard Optical Duplicate)",
                "BCFtools Call Node (Multiallelic Variant Calling)",
                "BCFtools Filter Node (Depth & Quality Threshold)",
                "Variant Visualization Node (SNV/Indel Distribution)",
                "Variant Report Node (Ti/Tv & Variant Summary)"
            ]
        }
    ]

    total_leaf_nodes = sum(len(d["nodes"]) for d in domain_data)
    print(f"Total domains: {len(domain_data)}, Total leaf nodes: {total_leaf_nodes}")

    # Set up coordinates
    # Height per leaf node: ~18-20 pixels
    fig_height = 24.0  # inches
    fig_width = 16.0   # inches
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    # Coordinate ranges
    # X: 0.0 (Root label) to 1.0 (Right edge)
    # Root node X: 0.12
    # Domain nodes X: 0.42
    # Leaf nodes X: 0.65
    x_root = 0.12
    x_domain = 0.42
    x_leaf = 0.62

    # Y coordinates: evenly space leaf nodes from bottom (0.02) to top (0.98)
    y_min, y_max = 0.02, 0.98
    y_leaves = np.linspace(y_max, y_min, total_leaf_nodes)

    # Compute domain Y coordinates as the center of their leaf nodes
    leaf_idx = 0
    for domain in domain_data:
        num_leaves = len(domain["nodes"])
        domain_leaf_ys = y_leaves[leaf_idx : leaf_idx + num_leaves]
        domain["y"] = np.mean(domain_leaf_ys)
        domain["leaf_ys"] = domain_leaf_ys
        leaf_idx += num_leaves

    # Root Y coordinate is center of all domains
    y_root = np.mean([d["y"] for d in domain_data])

    # Styling colors matching reference
    # Reference image uses soft slate/gray curves and blue-bordered circles
    curve_color = '#b0bec5'
    curve_alpha = 0.85
    curve_width = 1.0

    circle_edge = '#2563eb'       # Royal Blue border
    circle_fill = '#dbeafe'       # Soft light blue fill
    root_circle_edge = '#1d4ed8'  # Deep Blue
    root_circle_fill = '#bfdbfe'

    def draw_bezier(x1, y1, x2, y2, color=curve_color, width=curve_width, alpha=curve_alpha):
        # Cubic Bézier curve between (x1, y1) and (x2, y2)
        dx = x2 - x1
        ctrl1 = (x1 + dx * 0.5, y1)
        ctrl2 = (x2 - dx * 0.5, y2)
        
        verts = [(x1, y1), ctrl1, ctrl2, (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor='none', edgecolor=color, lw=width, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # 1. Draw curves from Root to Domains
    for domain in domain_data:
        draw_bezier(x_root, y_root, x_domain, domain["y"], color='#94a3b8', width=1.3, alpha=0.9)

    # 2. Draw curves from Domains to Leaves
    leaf_counter = 0
    for domain in domain_data:
        for y_l in domain["leaf_ys"]:
            draw_bezier(x_domain, domain["y"], x_leaf, y_l, color='#cbd5e1', width=0.85, alpha=0.85)
            leaf_counter += 1

    # 3. Draw Root Node
    ax.scatter([x_root], [y_root], s=140, facecolor=root_circle_fill, edgecolor=root_circle_edge, linewidth=2.0, zorder=5)
    ax.text(x_root - 0.015, y_root, "ComfyBIOWMS Nodes", 
            ha='right', va='center', fontsize=12, fontweight='bold', color='#0f172a', zorder=6)

    # 4. Draw Domain Nodes & Labels
    for domain in domain_data:
        ax.scatter([x_domain], [domain["y"]], s=90, facecolor=circle_fill, edgecolor=circle_edge, linewidth=1.6, zorder=5)
        # Domain label to the left of the circle
        label_text = f"{domain['name']} ({domain['count']})"
        ax.text(x_domain - 0.012, domain["y"], label_text, 
                ha='right', va='center', fontsize=9.5, fontweight='bold', color='#1e293b', zorder=6)

    # 5. Draw Leaf Nodes & Labels
    for domain in domain_data:
        for idx, (y_l, leaf_name) in enumerate(zip(domain["leaf_ys"], domain["nodes"])):
            ax.scatter([x_leaf], [y_l], s=35, facecolor=circle_fill, edgecolor=circle_edge, linewidth=1.1, zorder=5)
            # Leaf label to the right of the circle
            ax.text(x_leaf + 0.008, y_l, leaf_name, 
                    ha='left', va='center', fontsize=6.8, color='#334155', zorder=6)

    # Set limits and clean axes
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    
    out_path = 'figures/fig2_node_domain_tree.png'
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', pad_inches=0.1, facecolor='#ffffff')
    plt.close()
    print(f"Tree diagram successfully saved to {out_path}")

if __name__ == '__main__':
    create_tree_diagram()
