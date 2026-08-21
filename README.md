# ComfyBIOWMS: High-Throughput Bioinformatics Workflow Management System in ComfyUI

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![ComfyUI Custom Node](https://img.shields.io/badge/ComfyUI-Custom%20Node-green.svg)](https://github.com/comfyanonymous/ComfyUI)

**ComfyBIOWMS** is a full-featured, high-throughput Bioinformatics Workflow Management System natively integrated as a ComfyUI custom node pack. It transforms ComfyUI's interactive node graph into an enterprise-grade scientific workflow execution engine capable of running end-to-end multi-omics, structural biology, and publication-ready visualization pipelines.

---

## 🌟 Key Features

- **170+ Native Bioinformatics & CADD Nodes**: Covering genomics, transcriptomics, epigenomics, single-cell/spatial, metagenomics, proteomics, metabolomics, and structural biology (AlphaFold / ESMFold / DiffDock).
- **20+ Publication-Ready Visualizers**: Built-in 300+ DPI scientific plot renderers (Volcano plots, Manhattan plots, UMAPs, Clustermaps, OncoPrints, Kaplan-Meier curves, Phylogenetic trees, Ramachandran plots, etc.) emitting standard ComfyUI `IMAGE` tensors.
- **Isolated Conda/Mamba Runtime Bridge**: Automatic process isolation across 11+ specialized Conda environments without dependency conflicts.
- **ComfyUI Web Extension**: Domain-aware node color themes, badges, and real-time visualization previewers on the canvas.
- **Reproducible Artifact Tracking & Sidecars**: Manifest logging, QC tracking, and execution metrics for every pipeline stage.

---

## 📂 Project Structure

```text
ComfyBIOWMS/
├── __init__.py                     # ComfyUI root custom node entry point
├── pyproject.toml                  # Packaging & Comfy Registry v2 specifications
├── requirements.txt                # Core Python dependencies
├── README.md                       # Comprehensive documentation
├── LICENSE                         # MIT License
├── .gitignore                      # Git ignore rules for bio data and caches
│
├── nodes/                          # 170+ Bioinformatics custom node classes
│   ├── __init__.py                 # Subpackage mappings
│   ├── registry.py                 # Node registry & dynamic generation loader
│   ├── execution.py                # Runtime bridge & Conda execution
│   ├── sample_loading.py           # Sample metadata validator & parser
│   ├── visualizer_rendering.py     # Publication-grade plot renderer
│   ├── ref_nodes.py                # Core Bulk RNA-Seq & scRNA nodes
│   ├── variant_nodes.py            # DNA Variant Calling & BCFtools nodes
│   ├── atac_nodes.py               # ATAC-seq & MACS3 peak calling nodes
│   ├── metagenome_nodes.py         # Metagenomics & Kraken2/Bracken nodes
│   ├── assembly_nodes.py           # De Novo Genome Assembly & SPAdes/QUAST
│   ├── publication_visualizer_nodes.py # 20 Publication Plot Visualizers
│   ├── biopython_nodes.py          # 25 Biopython & Sequence Operation nodes
│   ├── genomics_longread_nodes.py  # 18 Long-Read Genomics & SV nodes
│   ├── transcriptomics_spatial_nodes.py # 24 Single-Cell & Spatial nodes
│   ├── epigenomics_nodes.py        # 14 Epigenomics & Functional Screen nodes
│   ├── proteomics_metabolomics_nodes.py # 15 Proteomics & Mass Spec nodes
│   ├── cadd_structural_nodes.py    # 18 AlphaFold, ESMFold, Docking & MD nodes
│   └── microbiome_nodes.py         # 16 Microbiome, Virome & Pathogen nodes
│
├── web/                            # ComfyUI Frontend Extensions
│   ├── js/comfybio_main.js         # LiteGraph node styling & badges
│   ├── css/comfybio.css            # Custom UI stylesheet
│   └── README.md
│
├── workflows/                      # Ready-to-use ComfyUI workflow JSONs
│   ├── 01_bulk_rnaseq_deseq2.json
│   ├── 02_variant_calling_bcftools.json
│   ├── 03_publication_visualizers.json
│   └── README.md
│
├── engine/                         # BioFlow Execution Engine Framework
│   ├── envs/                       # 12 Conda environment specifications
│   ├── scripts/                    # R and Python backend scripts
│   └── src/bioflow/                # Conda execution runtime & benchmark engine
│
├── tests/                          # Unit and integration test suite
│   ├── test_graph_checks.py
│   ├── test_node_mappings.py
│   └── test_web_directory.py
│
├── data/                           # Sample datasets and references
└── paper_platform/                 # Benchmark scripts & platform figures
```

---

## 🧬 Node Domains & Classification (170+ Nodes)

| Category | Node Count | Key Tools & Nodes |
| :--- | :---: | :--- |
| **Core & Reference RNA-Seq** | 20 | FastQC, Fastp, Salmon, Tximport, DESeq2, Scanpy, 10x |
| **DNA Variant Calling** | 8 | BWA-MEM2, MarkDuplicates, BCFtools (Call/Filter), VariantQC |
| **Epigenomics & ATAC-Seq** | 23 | MACS3, DeepTools, TOBIAS, Homer, MemeSuite, SEACR, Hi-C Cooler |
| **Metagenome & Microbiome** | 22 | Kraken2, Bracken, Scikit-Bio, DADA2, HUMAnN3, MetaPhlAn4, GeNomad |
| **Genome Assembly & Long-Read**| 24 | SPAdes, QUAST, Flye, Hifiasm, Sniffles2, CuteSV, DeepVariant, Pysam |
| **Single-Cell & Spatial** | 24 | AnnData, scVI, SCANVI, scVelo, CellRank, Squidpy, Tangram, Cell2location |
| **Publication Visualizers** | 20 | Volcano Plot, Manhattan Plot, Clustermap, GSEA, OncoPrint, Q-Q Plot, Sankey |
| **Biopython & Sequence Ops** | 25 | SeqIO, AlignIO, BLAST, Entrez, BioPDB, Primer3, Logomaker, Biotite |
| **Proteomics & Metabolomics** | 15 | PyOpenMS, Pyteomics, Matchms, Spec2Vec, MassQL, DIA-NN, MSFragger |
| **CADD & Structural Biology** | 18 | ColabFold (AlphaFold2), ESMFold, DiffDock, GNINA, RDKit, OpenMM, AutoDock Vina |

---

## 🚀 Quick Start

### 1. Installation in ComfyUI
Clone this repository directly into your ComfyUI `custom_nodes/` folder:

```bash
cd ComfyUI/custom_nodes
git clone https://github.com/DDUKHAE/ComfyBioWMS.git
cd ComfyBioWMS
pip install -r requirements.txt
```

### 2. Conda/Mamba Managed Environments
ComfyBIOWMS executes heavy command-line bioinformatics binaries inside isolated Conda environments. You can create the needed environments from `engine/envs/`:

```bash
# Example: Create Bulk RNA-Seq environment
conda env create -f engine/envs/bulk_rna_seq.yaml

# Example: Create Variant Analysis environment
conda env create -f engine/envs/variant_analysis.yaml
```

### 3. Launch ComfyUI & Load Workflows
1. Start your ComfyUI server:
   ```bash
   python main.py
   ```
2. Open `http://127.0.0.1:8188` in your browser.
3. Drag & drop any workflow from `workflows/` (e.g., `01_bulk_rnaseq_deseq2.json`) to run your first pipeline!

---

## 🧪 Testing

Run pytest from the repository root:

```bash
pytest
```

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
