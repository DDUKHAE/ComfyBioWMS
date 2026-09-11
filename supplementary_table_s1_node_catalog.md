# Supplementary Table S1: Catalog of ComfyBIOWMS Custom Nodes (98 Active Registered Nodes)

This table provides the functional classification and the evidence-backed verification tier for all 98 registered custom nodes in ComfyBIOWMS.
All registered nodes are functionally active with genuine computational or visualization execution (0 synthetic fallbacks):
- **CLI Execution Wrappers (58 nodes)**: Dispatches external bioinformatics CLI binaries with execution logging.
- **In-Process Python/R Scientific Nodes (40 nodes)**: Executes authentic scientific algorithms (Scanpy, AnnData, Biopython, RDKit, Bioconductor DESeq2) or renders publication figure tensors.

**Verification tiers (re-mapped 2026-09-11).** Each tier states the evidence actually obtained, not the node's execution class:

- **L2 · 독립 기준선 수치 동등성 (Independent numerical concordance)** — 7 nodes. Node outputs were compared numerically against an independent baseline (native `fastp`/`salmon` CLI and Bioconductor `tximport`/`DESeq2`, or an independent Scanpy reference script). These are the nodes taking part in the two quantitatively validated transcriptome case studies.
- **L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke)** — 55 nodes. The node was executed on fixture data and its outputs were checked for existence, non-emptiness, and format structure (SAM header/records, FASTA contigs, VCF/BED fields, report columns, rendered PNG), on top of the registration contract. No independent numerical baseline was run, so these do not support an equivalence claim.
- **L1 · 등록 계약만 검증 (Registration contract only)** — 36 nodes. Covered by the registry-wide contract test (`tests/test_comfyui_integration.py`: `INPUT_TYPES`/`RETURN_TYPES`/`FUNCTION`/`CATEGORY` declarations for all 98 nodes), with no dedicated execution test. Execution outputs are unevaluated.

A previous revision of this table labelled all 40 in-process nodes "Native Concordance Verified" and all 58 CLI wrappers "Execution & Interface Verified". That column restated the execution class rather than the evidence and overstated numerical verification; the tiers above replace it.

| Node Class Name | Display Name | Category | Functional Type | Verification Tier | Description |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `AnnDataInspect` | AnnData: Inspect H5AD Dataset | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `AutoDockVina` | AutoDock Vina: Molecular Docking | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `BBSplit` | BBTools: BBSplit Host/Contaminant Filter | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Bakta` | Bakta: Bacterial Genome Annotation | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `BcftoolsCall` | bcftools: Call Variants | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BcftoolsFilter` | bcftools: Filter Variants | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BcftoolsMpileup` | bcftools: Generate Pileup | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BedGraphToBigWig` | UCSC: bedGraph to BigWig | `ComfyBIO/Visualization` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BedtoolsGenomeCoverage` | BEDTools: BAM to bedGraph Coverage | `ComfyBIO/Visualization` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonAlignmentStats` | Biopython: Alignment Statistics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonGCContent` | Biopython: GC Content & Sequence Metrics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonPairwiseAlign` | Biopython: Pairwise Sequence Alignment | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonProtParam` | Biopython: Protein Physicochemical Properties | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonRestrictionDigest` | Biopython: Restriction Enzyme Digestion | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonSeqFilter` | Biopython: Sequence Length & GC Filter | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonSeqIOStats` | Biopython: Sequence File Statistics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BiopythonSeqTransform` | Biopython: Sequence Transformation (RevComp/Transcribe/Translate) | `ComfyBIO/Biopython` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Bowtie2Align` | Bowtie2: Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Bowtie2Build` | Bowtie2: Build Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `Bracken` | Bracken: Re-estimate Abundance | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BwaMem2Align` | BWA-MEM2: Align Reads | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `BwaMem2Index` | BWA-MEM2: Index Reference | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `CatFastq` | FastQ: Concatenate Split Files | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `ChipAtacCoverageProfile` | Visualization: ChIP/ATAC Coverage Profile | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `ClustermapHeatmap` | Visualization: Clustermap Heatmap | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `ColabFold` | ColabFold: AlphaFold2 Structure Prediction | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `Cyvcf2Stats` | cyvcf2: Variant File Statistics | `ComfyBIO/Genomics` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `DESeq2` | DESeq2: Differential Gene Expression | `ComfyBIO/RNA-Seq` | Subprocess CLI Execution Wrapper | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `DESeq2SampleQC` | RNA-Seq: Sample PCA & Correlation QC | `ComfyBIO/Quality Control` | In-Process Python/Bio Library | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `DeeptoolsMatrix` | deepTools: Compute Matrix & Profile | `ComfyBIO/Epigenomics` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `DupRadar` | dupRadar: Expression vs Duplication QC | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `ESMFold` | ESMFold: Fast Protein Structure Prediction | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `FastQC` | FastQC: Read Quality Report | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Fastp` | fastp: Trim and QC Reads | `ComfyBIO/Read Preprocessing` | Subprocess CLI Execution Wrapper | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `FlyeAssemble` | Flye: De Novo Long-read Assembler | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `GSEAPathway` | GSEApy: Gene Set Enrichment Analysis | `ComfyBIO/Pathways` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `GseaEnrichmentPlot` | Visualization: GSEA Enrichment Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `HISAT2Align` | HISAT2: Spliced Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `HISAT2Build` | HISAT2: Build Reference Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `HUMAnN` | HUMAnN: Functional Metabolic Profiler | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `HifiasmAssemble` | Hifiasm: HiFi Haplotype-resolved Assembler | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `InferStrandedness` | RNA-Seq: Infer Library Strandedness | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `JoinPaths` | Join Paths (multi-sample collector) | `ComfyBIO/Utility` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `KallistoQuant` | Kallisto: Pseudoalignment & Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `KaplanMeierSurvival` | Visualization: Kaplan-Meier Survival Curve | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Kraken2Classify` | Kraken2: Taxonomic Classification | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `LinkageDisequilibrium` | Visualization: Linkage Disequilibrium Heatmap | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `MaPlot` | Visualization: MA Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `ManhattanPlot` | Visualization: Manhattan Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `MdTrajectoryPlotter` | Visualization: MD Trajectory RMSD/RMSF Plotter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `MetaPhlAn` | MetaPhlAn: Taxonomic Profiler | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `MicrobiomeStackedBar` | Visualization: Microbiome Stacked Bar | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `Mosdepth` | Mosdepth: Fast BAM Coverage & Depth | `ComfyBIO/Genomics` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `MultiQC` | MultiQC: Comprehensive Analysis QC Report | `ComfyBIO/Reporting` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `OncoPrint` | Visualization: OncoPrint Mutation Landscape | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `PcoaScatter` | Visualization: Microbiome PCoA Scatter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `PhylogeneticTree` | Visualization: Phylogenetic Tree | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `PicardMarkDuplicates` | Picard: Mark Duplicates | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Preseq` | Preseq: Library Complexity Estimation | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Prokka` | Prokka: Rapid Prokaryotic Genome Annotation | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `ProteinLigandInteraction` | Visualization: Protein-Ligand Interaction | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `PybedtoolsIntersect` | pybedtools: Interval Intersect | `ComfyBIO/Genomics` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `PysamStats` | Pysam: Alignment Statistics & Summary | `ComfyBIO/Genomics` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `QqPlot` | Visualization: QQ Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `QualimapRNASeq` | Qualimap: RNA-seq Quality Control | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Quast` | QUAST: Evaluate Assembly | `ComfyBIO/Genome Assembly` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `RDKitDescriptor` | RDKit: Chemical Descriptors & Lipinski Rules | `ComfyBIO/CADD` | In-Process Python/Bio Library | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `RSEMCalculateExpression` | RSEM: Calculate Expression from BAM | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `RSeQCGeneBodyCoverage` | RSeQC: Gene Body Coverage | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `RSeQCInferExperiment` | RSeQC: Infer Experiment Strandedness | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `RSeQCJunctionSaturation` | RSeQC: Junction Saturation | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `RamachandranPlot` | Visualization: Ramachandran Dihedral Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `STARAlignReads` | STAR: Spliced Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `STARGenomeGenerate` | STAR: Generate Genome Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SalmonIndex` | Salmon: Transcriptome Index | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SalmonQuantAlignment` | Salmon: Transcriptome BAM Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SalmonQuantReads` | Salmon: FastQ Quasi-mapping Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `SamtoolsIndex` | samtools: Index Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SamtoolsMarkdup` | samtools: Mark Duplicates | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SamtoolsSort` | samtools: Sort Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SankeyCellFate` | Visualization: Sankey Cell Fate Alluvial | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `ScanpyCluster` | Scanpy: PCA, UMAP & Leiden Clustering | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `ScanpyNormalize` | Scanpy: Normalize & Highly Variable Genes | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `ScanpyQC` | Scanpy: Single-cell Quality Control | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `SeqKitStats` | SeqKit: Sequence File Statistics | `ComfyBIO/Genomics` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `SortMeRNA` | SortMeRNA: Ribosomal RNA Filter | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Spades` | SPAdes: Assemble Genome | `ComfyBIO/Genome Assembly` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SpatialTissueOverlay` | Visualization: Spatial Tissue Overlay | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `StringTie` | StringTie: Transcript Assembly & Quant | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SylphProfile` | Sylph: Fast Metagenomic Profiler | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `SyntenyGenome` | Visualization: Synteny Genome Alignment | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `TrimGalore` | Trim Galore!: Quality & Adapter Trimming | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `Trimmomatic` | Trimmomatic: Read Trimming | `ComfyBIO/Read Preprocessing` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `Tximport` | tximport: Transcript to Gene Count Matrix | `ComfyBIO/Quantification` | In-Process Python/Bio Library | L2 · 독립 기준선 수치 동등성 (Independent numerical concordance) |  |
| `UmapScatter` | Visualization: UMAP Scatter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `UmiToolsDedup` | UMI-tools: Deduplicate BAM | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | L1 · 등록 계약만 검증 (Registration contract only) |  |
| `UmiToolsExtract` | UMI-tools: Extract UMI to Header | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
| `VolcanoPlot` | Visualization: Volcano Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | L1+S · 실행 계약 + 출력 스모크 (Contract + output smoke) |  |
