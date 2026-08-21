# 생물정보학 분석 도메인 및 세부 워크플로우 명세 (8대 도메인 24개 워크플로우)

본 문서는 ComfyBioWMS 플랫폼이 지원 및 지향하는 **8대 핵심 생물정보학 도메인**과 **24개 세부 표준 워크플로우**의 입력, 단계별 도구, 소켓 타입 계약, 출력 및 출판급 시각화 노드 명세를 정의합니다.

---

## 1. 전사체학 (Transcriptomics)

### 1.1 Bulk mRNA-Seq 전사체 정량 및 차등 발현 분석
- **입력:** Raw Paired-End FASTQ, Sample Metadata CSV, Transcriptome FASTA
- **단계 및 도구:**
  1. `FastpQCNode`: 어댑터 트리밍 및 QC JSON (`FASTQ_PAIR` $\rightarrow$ `FASTP_QC_JSON`)
  2. `SalmonIndexNode` & `SalmonQuantNode`: 미끼 서열(Decoy) 기반 준정렬 정량 (`SALMON_INDEX` $\rightarrow$ `SALMON_QUANT_DIR`)
  3. `TximportNode`: 전사체-유전자 수준 카운트 행렬 변환 (`SALMON_QUANT_DIR` $\rightarrow$ `COUNT_MATRIX_CSV`)
  4. `DESeq2AnalysisNode`: 음이항 분포 기반 Wald 차등 발현 검정 (`COUNT_MATRIX_CSV` $\rightarrow$ `DEG_RESULTS_CSV`)
  5. `VolcanoPlotVisualizerNode` & `ClusteredHeatmapVisualizerNode`: 300 DPI 출판용 도판 렌더링 (`DEG_RESULTS_CSV` $\rightarrow$ `IMAGE_TENSOR`)

### 1.2 Total RNA-Seq 및 신규 전사체 조립 (rRNA Depleted)
- **입력:** Raw Stranded FASTQ, Reference Genome FASTA & GTF
- **단계 및 도구:** `FastpQCNode` $\rightarrow$ `SortMeRNANode` (rRNA 제거) $\rightarrow$ `STARAlignNode` / `HISAT2Node` (스플라이스 정렬) $\rightarrow$ `StringTieNode` (신규 전사체 어셈블리) $\rightarrow$ `GffCompareNode` $\rightarrow$ `SplicingEventVisualizerNode`.

### 1.3 Small RNA / miRNA 시퀀싱 분석
- **입력:** Raw Single-End FASTQ (Short reads), miRBase FASTA
- **단계 및 도구:** `FastpShortTrimNode` $\rightarrow$ `miRTraceQCNode` $\rightarrow$ `BowtieAlignNode` $\rightarrow$ `miRDeep2QuantNode` $\rightarrow$ `TargetScanPredictNode` $\rightarrow$ `MiRNAExpressionBarplotNode`.

### 1.4 단일세포 전사체학 (Single-Cell RNA-Seq)
- **입력:** 10x Genomics Chromium FASTQ / Feature-Barcode Matrix
- **단계 및 도구:** `CellRangerCountNode` / `AlevinFryNode` $\rightarrow$ `ScanpyQCFilterNode` (미토콘드리아 유전자, 이중체 제거) $\rightarrow$ `NormalizeAndLog1pNode` $\rightarrow$ `HighlyVariableGenesNode` $\rightarrow$ `PCANode` $\rightarrow$ `LeidenClusteringNode` $\rightarrow$ `UMAPVisualizerNode` & `MarkerGeneViolinVisualizerNode`.

### 1.5 공간 전사체학 (Spatial Transcriptomics)
- **입력:** 10x Visium FASTQ, 고해상도 조직 H&E 슬라이드 이미지
- **단계 및 도구:** `SpaceRangerCountNode` $\rightarrow$ `SquidpySpatialQCNode` $\rightarrow$ `SpatialDomainClusteringNode` $\rightarrow$ `SpatialGeneOverlayVisualizerNode` (조직 절편 위 발현 핫스팟 중첩 렌더링).

---

## 2. 유전체학 및 변이 분석 (Genomics & Variant Discovery)

### 2.1 생식세포 변이 분석 (Germline WGS/WES)
- **입력:** Raw FASTQ / CRAM, Reference Genome (GRCh38), Known Sites VCF
- **단계 및 도구:**
  1. `FastpQCNode`: 서열 전처리
  2. `BwaMem2AlignNode`: 고속 유전체 정렬 (`FASTQ_PAIR` $\rightarrow$ `SAM_BAM`)
  3. `MarkDuplicatesNode`: 광학 및 PCR 중복 리드 마킹 (`SAM_BAM` $\rightarrow$ `DEDUP_BAM`)
  4. `BcftoolsCallNode` / `GATKHaplotypeCallerNode`: SNV 및 Small InDel 변이 탐지 (`DEDUP_BAM` $\rightarrow$ `RAW_VCF`)
  5. `BcftoolsFilterNode`: 품질 점수 및 깊이(DP) 필터링 (`RAW_VCF` $\rightarrow$ `FILTERED_VCF`)
  6. `SnpEffAnnotationNode` / `VEPAnnotationNode`: 변이 기능 주석 (`FILTERED_VCF` $\rightarrow$ `ANNOTATED_VCF`)
  7. `ManhattanPlotVisualizerNode` & `VariantDensityPlotNode` (`IMAGE_TENSOR`).

### 2.2 체세포 변이 분석 (Somatic Tumor-Normal Mutation)
- **입력:** 종양-정상 페어 FASTQ, Target BED (Exome)
- **단계 및 도구:** `BwaMem2AlignNode` $\rightarrow$ `MarkDuplicatesNode` $\rightarrow$ `Mutect2PairedCallNode` $\rightarrow$ `LearnReadOrientationModelNode` $\rightarrow$ `FilterMutectCallsNode` $\rightarrow$ `OncoPrintVisualizerNode` & `MutationSignatureVisualizerNode`.

### 2.3 구조적 변이 및 복제수 변이 (SV & CNV)
- **입력:** Dedup BAM, Reference Genome
- **단계 및 도구:** `MantaSVNode` $\rightarrow$ `DellySVNode` $\rightarrow$ `CNVkitNode` $\rightarrow$ `CircosPlotVisualizerNode` & `GenomeWideCNVTrackNode`.

### 2.4 롱리드 유전체학 (PacBio HiFi / Oxford Nanopore)
- **입력:** PacBio HiFi / ONT FASTQ, Reference Genome
- **단계 및 도구:** `NanoPlotQCNode` $\rightarrow$ `Minimap2LongAlignNode` $\rightarrow$ `Clair3VariantCallNode` $\rightarrow$ `Sniffles2SVNode` $\rightarrow$ `PhasedHaplotypeReportNode`.

---

## 3. 후성유전학 및 크로마틴 동역학 (Epigenomics)

### 3.1 ATAC-Seq 크로마틴 접근성 분석
- **입력:** Paired-End ATAC-Seq FASTQ, Reference Genome
- **단계 및 도구:** `FastpTrimNode` $\rightarrow$ `Bowtie2AlignNode` (미토콘드리아 배제) $\rightarrow$ `MarkDuplicatesNode` $\rightarrow$ `AtacQualityFilterNode` $\rightarrow$ `Macs3PeakCallingNode` (Narrow Peak) $\rightarrow$ `TssEnrichmentNode` $\rightarrow$ `GenomicTrackVisualizerNode`.

### 3.2 ChIP-Seq 전사인자 및 히스톤 수식화 분석
- **입력:** ChIP FASTQ + Control Input FASTQ, Reference Genome
- **단계 및 도구:** `Bowtie2AlignNode` $\rightarrow$ `MarkDuplicatesNode` $\rightarrow$ `Macs3BroadPeakNode` $\rightarrow$ `HomerMotifDiscoveryNode` $\rightarrow$ `PeakDistributionHeatmapNode`.

### 3.3 CUT&RUN / CUT&Tag 초저투입 후성유전체 프로파일링
- **입력:** Low-input Paired-End FASTQ, Spike-in Control
- **단계 및 도구:** `Bowtie2AlignNode` $\rightarrow$ `SeacrPeakCallingNode` $\rightarrow$ `DeepToolsCoverageNode` $\rightarrow$ `ChromatinStateVisualizerNode`.

### 3.4 전장 유전체 바이설파이트 시퀀싱 (WGBS / RRBS 메틸화)
- **입력:** Bisulfite-treated FASTQ, Bisulfite-converted Reference
- **단계 및 도구:** `TrimGaloreNode` $\rightarrow$ `BismarkAlignNode` $\rightarrow$ `BismarkMethylationExtractorNode` $\rightarrow$ `DMRAnalysisNode` $\rightarrow$ `MethylationDensityHeatmapNode`.

---

## 4. 메타유전체학 및 마이크로바이옴 (Metagenomics & Microbiome)

### 4.1 샷건 메타유전체 분류학적 프로파일링
- **입력:** Shotgun Metagenome FASTQ, Kraken2 Standard DB
- **단계 및 도구:** `FastpTrimNode` $\rightarrow$ `Kraken2ClassifyNode` $\rightarrow$ `BrackenAbundanceNode` $\rightarrow$ `TaxonomicStackedBarplotNode` & `KronaSunburstNode`.

### 4.2 기능적 메타유전체 및 대사 경로 분석
- **입력:** Host-depleted FASTQ, UniRef50 / MetaCyc DB
- **단계 및 도구:** `KneadDataHostFilterNode` $\rightarrow$ `Humann3PathwayNode` $\rightarrow$ `PathwayAbundanceHeatmapNode`.

### 4.3 16S/18S/ITS 앰플리콘 시퀀싱
- **입력:** Amplicon FASTQ, Primers, SILVA/Greengenes DB
- **단계 및 도구:** `Dada2QCFilterNode` $\rightarrow$ `Dada2DenoiseNode` (ASV 도출) $\rightarrow$ `Dada2AssignTaxonomyNode` $\rightarrow$ `AlphaDiversityBoxplotNode` $\rightarrow$ `BetaDiversityPCoAPlotNode`.

### 4.4 메타유전체 유래 유전체 조립 (MAGs Binning)
- **입력:** Shotgun FASTQ
- **단계 및 도구:** `MetaSpadesAssemblyNode` $\rightarrow$ `MetaBat2BinningNode` $\rightarrow$ `CheckM2QCNode` $\rightarrow$ `GtdbTkPhylogenyNode`.

---

## 5. 유전체 조립 및 기능 주석화 (De Novo Assembly & Annotation)

### 5.1 미생물 및 바이러스 De Novo 어셈블리
- **입력:** Microbial/Viral Paired-End FASTQ
- **단계 및 도구:** `FastpQCNode` $\rightarrow$ `SpadesAssemblyNode` / `UnicyclerNode` $\rightarrow$ `QuastEvaluationNode` $\rightarrow$ `BandageGraphVisualizerNode`.

### 5.2 진핵생물 대형 유전체 하이브리드 어셈블리
- **입력:** PacBio HiFi + Illumina Short Reads
- **단계 및 도구:** `HifiasmAssemblyNode` $\rightarrow$ `RaconPolishingNode` $\rightarrow$ `BuscoCompletenessNode` $\rightarrow$ `AssemblyMetricVisualizerNode`.

### 5.3 원핵생물 / 진핵생물 게놈 구조 및 기능 주석화
- **입력:** Assembled Contigs FASTA
- **단계 및 도구:** `BaktaAnnotationNode` / `ProkkaNode` $\rightarrow$ `EggNogFunctionalMapNode` $\rightarrow$ `SyntenyGenomePlotNode`.

---

## 6. 단백체학 및 대사체학 (Proteomics & Metabolomics)

### 6.1 LC-MS/MS Bottom-Up 단백체학 (DDA/DIA)
- **입력:** Raw Mass Spectrometry Files (.raw, .mzML), UniProt FASTA
- **단계 및 도구:** `MsconvertRawNode` $\rightarrow$ `MaxQuantDDAQuantNode` / `DiaNNEngineNode` $\rightarrow$ `PerseusStatsNode` $\rightarrow$ `ProteinVolcanoPlotNode` & `ProteinCorrelationHeatmapNode`.

### 6.2 비타겟 대사체학 (Untargeted Metabolomics)
- **입력:** LC-MS mzML Files
- **단계 및 도구:** `MZmine3PeakPickNode` $\rightarrow$ `FeatureAlignmentNode` $\rightarrow$ `GNPSMolecularNetworkNode` $\rightarrow$ `MetabolitePCAPlsdaPlotNode`.

---

## 7. 구조생물학 및 CADD 신약개발 (Structural Biology & Drug Discovery)

### 7.1 단백질 3차 구조 예측 및 품질 평가
- **입력:** Amino Acid FASTA
- **단계 및 도구:** `ESMFoldStructureNode` / `AlphaFold2Node` $\rightarrow$ `MolProbityAssessmentNode` $\rightarrow$ `PyMOL3DRenderNode` (`IMAGE_TENSOR`).

### 7.2 분자 도킹 및 가상 스크리닝
- **입력:** Target Protein PDB, Small Molecule SMILES/SDF
- **단계 및 도구:** `RDKitMoleculePrepNode` $\rightarrow$ `AutoDockVinaDockingNode` $\rightarrow$ `PLIPInteractionDiagramNode` $\rightarrow$ `BindingAffinityBarplotNode`.

---

## 8. 비교유전체학 및 계통발생학 (Comparative Genomics & Phylogenetics)

### 8.1 오솔로그 추론 및 판게놈(Pangenome) 분석
- **입력:** Multiple Species Protein FASTA / GFF
- **단계 및 도구:** `OrthoFinderNode` / `RoaryPangenomeNode` $\rightarrow$ `CoreAccessoryGenomeNode` $\rightarrow$ `PangenomeCurveVisualizerNode`.

### 8.2 다중 서열 정렬 및 계통수 재구축
- **입력:** Gene/Protein FASTA Set
- **단계 및 도구:** `MafftAlignmentNode` $\rightarrow$ `TrimAlFilterNode` $\rightarrow$ `IqTree2PhylogenyNode` $\rightarrow$ `PhylogeneticTreeVisualizerNode` (Circular/Rectangular Tree 렌더링).
