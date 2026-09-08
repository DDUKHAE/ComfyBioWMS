# Supplementary Table S1: Catalog of ComfyBIOWMS Custom Nodes (98 Active Registered Nodes)

This table provides the definitive functional classification and verification status for all 98 registered custom nodes in ComfyBIOWMS.
All registered nodes are functionally active with genuine computational or visualization execution (0 synthetic fallbacks):
- **CLI Execution Wrappers (58 nodes)**: Dispatches external bioinformatics CLI binaries with execution logging.
- **In-Process Python/R Scientific Nodes (40 nodes)**: Executes authentic scientific algorithms (Scanpy, AnnData, Biopython, RDKit, Bioconductor DESeq2) or renders publication figure tensors.

| Node Class Name | Display Name | Category | Functional Type | Verification Tier | Description |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `AnnDataInspect` | AnnData: Inspect H5AD Dataset | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `AutoDockVina` | AutoDock Vina: Molecular Docking | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `BBSplit` | BBTools: BBSplit Host/Contaminant Filter | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Bakta` | Bakta: Bacterial Genome Annotation | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `BcftoolsCall` | bcftools: Call Variants | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BcftoolsFilter` | bcftools: Filter Variants | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BcftoolsMpileup` | bcftools: Generate Pileup | `ComfyBIO/Variants` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BedGraphToBigWig` | UCSC: bedGraph to BigWig | `ComfyBIO/Visualization` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `BedtoolsGenomeCoverage` | BEDTools: BAM to bedGraph Coverage | `ComfyBIO/Visualization` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `BiopythonAlignmentStats` | Biopython: Alignment Statistics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonGCContent` | Biopython: GC Content & Sequence Metrics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonPairwiseAlign` | Biopython: Pairwise Sequence Alignment | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonProtParam` | Biopython: Protein Physicochemical Properties | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonRestrictionDigest` | Biopython: Restriction Enzyme Digestion | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonSeqFilter` | Biopython: Sequence Length & GC Filter | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonSeqIOStats` | Biopython: Sequence File Statistics | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BiopythonSeqTransform` | Biopython: Sequence Transformation (RevComp/Transcribe/Translate) | `ComfyBIO/Biopython` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `Bowtie2Align` | Bowtie2: Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Bowtie2Build` | Bowtie2: Build Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Bracken` | Bracken: Re-estimate Abundance | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `BwaMem2Align` | BWA-MEM2: Align Reads | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `BwaMem2Index` | BWA-MEM2: Index Reference | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `CatFastq` | FastQ: Concatenate Split Files | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `ChipAtacCoverageProfile` | Visualization: ChIP/ATAC Coverage Profile | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `ClustermapHeatmap` | Visualization: Clustermap Heatmap | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `ColabFold` | ColabFold: AlphaFold2 Structure Prediction | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Cyvcf2Stats` | cyvcf2: Variant File Statistics | `ComfyBIO/Genomics` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `DESeq2` | DESeq2: Differential Gene Expression | `ComfyBIO/RNA-Seq` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `DESeq2SampleQC` | RNA-Seq: Sample PCA & Correlation QC | `ComfyBIO/Quality Control` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `DeeptoolsMatrix` | deepTools: Compute Matrix & Profile | `ComfyBIO/Epigenomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `DupRadar` | dupRadar: Expression vs Duplication QC | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `ESMFold` | ESMFold: Fast Protein Structure Prediction | `ComfyBIO/CADD` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `FastQC` | FastQC: Read Quality Report | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `Fastp` | fastp: Trim and QC Reads | `ComfyBIO/Read Preprocessing` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `FlyeAssemble` | Flye: De Novo Long-read Assembler | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `GSEAPathway` | GSEApy: Gene Set Enrichment Analysis | `ComfyBIO/Pathways` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `GseaEnrichmentPlot` | Visualization: GSEA Enrichment Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `HISAT2Align` | HISAT2: Spliced Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `HISAT2Build` | HISAT2: Build Reference Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `HUMAnN` | HUMAnN: Functional Metabolic Profiler | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `HifiasmAssemble` | Hifiasm: HiFi Haplotype-resolved Assembler | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `InferStrandedness` | RNA-Seq: Infer Library Strandedness | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `JoinPaths` | Join Paths (multi-sample collector) | `ComfyBIO/Utility` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `KallistoQuant` | Kallisto: Pseudoalignment & Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `KaplanMeierSurvival` | Visualization: Kaplan-Meier Survival Curve | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `Kraken2Classify` | Kraken2: Taxonomic Classification | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `LinkageDisequilibrium` | Visualization: Linkage Disequilibrium Heatmap | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `MaPlot` | Visualization: MA Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `ManhattanPlot` | Visualization: Manhattan Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `MdTrajectoryPlotter` | Visualization: MD Trajectory RMSD/RMSF Plotter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `MetaPhlAn` | MetaPhlAn: Taxonomic Profiler | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `MicrobiomeStackedBar` | Visualization: Microbiome Stacked Bar | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `Mosdepth` | Mosdepth: Fast BAM Coverage & Depth | `ComfyBIO/Genomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `MultiQC` | MultiQC: Comprehensive Analysis QC Report | `ComfyBIO/Reporting` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `OncoPrint` | Visualization: OncoPrint Mutation Landscape | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `PcoaScatter` | Visualization: Microbiome PCoA Scatter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `PhylogeneticTree` | Visualization: Phylogenetic Tree | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `PicardMarkDuplicates` | Picard: Mark Duplicates | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Preseq` | Preseq: Library Complexity Estimation | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Prokka` | Prokka: Rapid Prokaryotic Genome Annotation | `ComfyBIO/Microbiome` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `ProteinLigandInteraction` | Visualization: Protein-Ligand Interaction | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `PybedtoolsIntersect` | pybedtools: Interval Intersect | `ComfyBIO/Genomics` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `PysamStats` | Pysam: Alignment Statistics & Summary | `ComfyBIO/Genomics` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `QqPlot` | Visualization: QQ Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `QualimapRNASeq` | Qualimap: RNA-seq Quality Control | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Quast` | QUAST: Evaluate Assembly | `ComfyBIO/Genome Assembly` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `RDKitDescriptor` | RDKit: Chemical Descriptors & Lipinski Rules | `ComfyBIO/CADD` | In-Process Python/Bio Library | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `RSEMCalculateExpression` | RSEM: Calculate Expression from BAM | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `RSeQCGeneBodyCoverage` | RSeQC: Gene Body Coverage | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `RSeQCInferExperiment` | RSeQC: Infer Experiment Strandedness | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `RSeQCJunctionSaturation` | RSeQC: Junction Saturation | `ComfyBIO/Quality Control` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `RamachandranPlot` | Visualization: Ramachandran Dihedral Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `STARAlignReads` | STAR: Spliced Read Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `STARGenomeGenerate` | STAR: Generate Genome Index | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SalmonIndex` | Salmon: Transcriptome Index | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SalmonQuantAlignment` | Salmon: Transcriptome BAM Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SalmonQuantReads` | Salmon: FastQ Quasi-mapping Quant | `ComfyBIO/Quantification` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SamtoolsIndex` | samtools: Index Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SamtoolsMarkdup` | samtools: Mark Duplicates | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `SamtoolsSort` | samtools: Sort Alignment | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SankeyCellFate` | Visualization: Sankey Cell Fate Alluvial | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `ScanpyCluster` | Scanpy: PCA, UMAP & Leiden Clustering | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `ScanpyNormalize` | Scanpy: Normalize & Highly Variable Genes | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `ScanpyQC` | Scanpy: Single-cell Quality Control | `ComfyBIO/Single-Cell` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SeqKitStats` | SeqKit: Sequence File Statistics | `ComfyBIO/Genomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `SortMeRNA` | SortMeRNA: Ribosomal RNA Filter | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Spades` | SPAdes: Assemble Genome | `ComfyBIO/Genome Assembly` | Subprocess CLI Execution Wrapper | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `SpatialTissueOverlay` | Visualization: Spatial Tissue Overlay | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `StringTie` | StringTie: Transcript Assembly & Quant | `ComfyBIO/Assembly` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `SylphProfile` | Sylph: Fast Metagenomic Profiler | `ComfyBIO/Metagenomics` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `SyntenyGenome` | Visualization: Synteny Genome Alignment | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `TrimGalore` | Trim Galore!: Quality & Adapter Trimming | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Trimmomatic` | Trimmomatic: Read Trimming | `ComfyBIO/Read Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `Tximport` | tximport: Transcript to Gene Count Matrix | `ComfyBIO/Quantification` | In-Process Python/Bio Library | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `UmapScatter` | Visualization: UMAP Scatter | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
| `UmiToolsDedup` | UMI-tools: Deduplicate BAM | `ComfyBIO/Alignment` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `UmiToolsExtract` | UMI-tools: Extract UMI to Header | `ComfyBIO/Preprocessing` | Subprocess CLI Execution Wrapper | 실제 실행 / 인터페이스 확인 (Execution & Interface Verified) |  |
| `VolcanoPlot` | Visualization: Volcano Plot | `ComfyBIO/Visualization` | In-Process Scientific Visualizer (Matplotlib/PyTorch) | 원 도구 동등성 확인 (Native Concordance Verified) |  |
