# ComfyBioflow 노드 생태계 진화계통수 (Node Phylogenetic Tree)

본 문서는 ComfyBioflow 플랫폼에 구현된 **150여 개 전 도메인 노드(Bioinformatics Nodes)**의 계층적 연관성과 기능적 분화를 생물학의 **진화계통수(Phylogenetic Tree / Cladogram)** 형태로 체계화한 노드 분류 카탈로그입니다.

```
═══════════════════════════════════════════════════════════════════════════════════
                      🌿 ROOT: ComfyBioflow Node Ecosystem
═══════════════════════════════════════════════════════════════════════════════════
```

---

## 🌳 전체 노드 계통수 개요 (Cladogram Overview)

```text
[ROOT: ComfyBioflow Engine]
 │
 ├── 1. [Kingdom] Transcriptomics & Single-Cell / Spatial (전사체학 및 단일세포/공간)
 │    ├── 1.1 [Phylum] Bulk RNA-Seq Core (대용량 전사체 정량 및 차등발현)
 │    ├── 1.2 [Phylum] Single-Cell RNA-Seq (단일세포 분해능 오믹스)
 │    ├── 1.3 [Phylum] Spatial Transcriptomics (공간 전사체학)
 │    └── 1.4 [Phylum] Isoform & Long-Read Transcriptomics (롱리드 동형전사체)
 │
 ├── 2. [Kingdom] Genomics & Variant Discovery (유전체학 및 변이/롱리드)
 │    ├── 2.1 [Phylum] Short-Read Alignment & Indexing (숏리드 정렬 및 전처리)
 │    ├── 2.2 [Phylum] Variant Calling & Filtering (변이 탐지 및 필터링)
 │    ├── 2.3 [Phylum] Long-Read Genomics & SV Discovery (롱리드 및 구조변이)
 │    └── 2.4 [Phylum] High-Throughput Interval Operations (게놈 구간 및 SAM/BAM 연산)
 │
 ├── 3. [Kingdom] Epigenomics & 3D Chromatin (후성유전학 및 3차원 염색질)
 │    ├── 3.1 [Phylum] Chromatin Accessibility (ATAC-Seq 크로마틴 접근성)
 │    ├── 3.2 [Phylum] Protein-DNA Interactions (ChIP-Seq / CUT&Tag / Motif)
 │    ├── 3.3 [Phylum] 3D Genome Architecture (Hi-C / Micro-C 염색질 입체구조)
 │    └── 3.4 [Phylum] Epigenetic Modifications & Screens (메틸화 및 기능 스크린)
 │
 ├── 4. [Kingdom] Metagenomics, Microbiome & Virology (메타유전체학 및 바이러스학)
 │    ├── 4.1 [Phylum] Taxonomic Profiling (분류학적 동정 및 풍부도 추정)
 │    ├── 4.2 [Phylum] Amplicon & Diversity (16S/ITS ASV 및 알파/베타 다양성)
 │    ├── 4.3 [Phylum] Functional Metagenomics (기능적 유전자 및 대사 경로)
 │    └── 4.4 [Phylum] Virome & Phage Discovery (바이러스 및 박테리오파지 발굴)
 │
 ├── 5. [Kingdom] De Novo Genome Assembly & Annotation (유전체 조립 및 주석화)
 │    ├── 5.1 [Phylum] Microbial & Viral Assembly (미생물 서열 조립)
 │    ├── 5.2 [Phylum] Long-Read & Hybrid Assembly (롱리드/하이브리드 조립)
 │    ├── 5.3 [Phylum] Assembly QC & Graph Analytics (어셈블리 품질 및 조립 그래프)
 │    └── 5.4 [Phylum] Genome Structural & Functional Annotation (게놈 주석화)
 │
 ├── 6. [Kingdom] Proteomics & Metabolomics (단백체학 및 대사체학)
 │    ├── 6.1 [Phylum] Bottom-Up Proteomics (LC-MS/MS DDA/DIA 단백체)
 │    ├── 6.2 [Phylum] Mass Spectrometry Feature Engineering (질량분석 전처리)
 │    ├── 6.3 [Phylum] Untargeted & Lipid Metabolomics (비타겟 대사체 및 지질체)
 │    └── 6.4 [Phylum] Immunopeptidomics & Structural MS (면역 펩타이드 및 구조 MS)
 │
 ├── 7. [Kingdom] Structural Biology & CADD (구조생물학 및 컴퓨터 신약개발)
 │    ├── 7.1 [Phylum] AI 3D Structure Prediction (인공지능 단백질 3차 구조 예측)
 │    ├── 7.2 [Phylum] Molecular Docking & Pocket Finding (분자 도킹 및 바인딩 포켓)
 │    ├── 7.3 [Phylum] Molecular Dynamics & Free Energy (분자동역학 및 자유에너지)
 │    └── 7.4 [Phylum] Cheminformatics & Molecular ML (화합물 정보학 및 분자 ML)
 │
 ├── 8. [Kingdom] BioPython & Sequence Utilities (생물정보 유틸리티 및 비교유전체)
 │    ├── 8.1 [Phylum] Sequence I/O & Feature Translation (서열 조작 및 번역)
 │    ├── 8.2 [Phylum] Homology, Motif & Domain Search (상동성/모티프 검색)
 │    ├── 8.3 [Phylum] Primer & Molecular Cloning Design (프라이머 및 분자생물학 도구)
 │    └── 8.4 [Phylum] Comparative Genomics & Phylogenetics (비교유전체 및 계통수)
 │
 └── 9. [Kingdom] Tier-1 Publication Visualizers (출판 품질 전문 시각화 노드)
      ├── 9.1 [Phylum] Expression & Statistical Analytics (발현/통계 도판)
      ├── 9.2 [Phylum] Single-Cell & Spatial Cell-Fate (단일세포/공간 도판)
      ├── 9.3 [Phylum] Genomics & Epigenomic Tracks (유전체/후성유전 트랙 도판)
      ├── 9.4 [Phylum] Microbiome & Phylogenetics (마이크로바이옴/계통 도판)
      └── 9.5 [Phylum] Structural & CADD Visualizers (구조생물학/도킹 도판)
```

---

## 🔬 세부 계통 및 노드 목록 (Detailed Node Catalog by Taxon)

### 1. [Kingdom] Transcriptomics & Single-Cell / Spatial

* **1.1 Bulk RNA-Seq Core**
  * `SampleMetadataValidatorNode` — 시퀀싱 샘플 메타데이터 및 그룹 조건 검증
  * `FastpQCNode` — 초고속 어댑터 트리밍 및 QC JSON 프로파일링
  * `FastpTrimNode` — FASTQ 품질 필터링 및 리드 트리밍
  * `FastQCNode` — 고전적 FASTQ 리드 품질 진단
  * `TrimmomaticTrimNode` — 슬라이딩 윈도우 기반 리드 트리밍
  * `SalmonIndexNode` — Decoy 서열 포함 전사체 인덱스 빌드
  * `SalmonQuantNode` — 전사체 발현량 준정렬(Quasi-mapping) 고속 정량
  * `KallistoQuantNode` — 의사정렬(Pseudoalignment) 기반 발현 정량
  * `TximportNode` — 전사체 $\rightarrow$ 유전자 수준 카운트 행렬 변환
  * `DESeq2AnalysisNode` — 음이항 분포 기반 차등 발현 유전자(DEG) 검정
  * `EdgeRAnalysisNode` — 정확 검정(Exact Test) 및 일반화 선형 모델 DEG 분석
  * `LimmaVoomAnalysisNode` — 정규화 카운트 기반 선형 모델 DEG 분석
  * `GseapyEnrichmentNode` — GO / KEGG / Reactome 유전자 세트 농축 분석

* **1.2 Single-Cell RNA-Seq**
  * `TenxCountNode` — 10x Genomics 바코드/UMI 디멀티플렉싱 및 카운트
  * `ScanpyQCNode` — 단일세포 미토콘드리아 유전자율 및 UMI 필터링
  * `ScanpyNormalizeNode` — 총 카운트 정규화 및 $\log(1+x)$ 변환
  * `ScanpyClusterNode` — 고변이 유전자(HVG) 선별, PCA 및 Leiden 클러스터링
  * `ScanpyMarkerGenesNode` — 클러스터별 특이적 마커 유전자 랭킹
  * `ScviToolsLatentNode` — 심층 생성 모델 기반 단일세포 배치 효과 교정
  * `ScanviCellTypeNode` — 준지도 학습 기반 세포 타입 자동 주석화
  * `ScVeloDynamicsNode` — RNA 벨로시티(RNA Velocity) 기반 세포 분화 동역학
  * `CellRankFateNode` — 마르코프 연쇄 기반 세포 분화 운명(Cell Fate) 궤적 추적
  * `PyScenicRegulonNode` — 전사인자-표적 유전자 조절 네트워크(Regulon) 추론
  * `CiteSeqCountNode` — 단백질 항체 태그(ADT) 바코드 카운트 및 멀티모달 처리
  * `MuonMultimodalNode` — CITE-seq, scATAC+scRNA 다중 오믹스 통합

* **1.3 Spatial Transcriptomics**
  * `SquidpySpatialNode` — 공간 전사체 이웃 그래프 및 공간 자기상관(Moran's I) 분석
  * `Cell2locationAbundanceNode` — 공간 스팟 내 단일세포 유형 밀도 베이지안 디컨볼루션
  * `TangramSpatialMappingNode` — 단일세포 전사체를 공간 조직 슬라이드에 정렬 매핑
  * `CellPhoneDbInteractionNode` — 공간적 거리를 고려한 세포-세포 리간드-수용체 상호작용

* **1.4 Isoform & Long-Read Transcriptomics**
  * `StringTie2AssembleNode` — 스플라이스 정렬 기반 신규 전사체 조립 및 정량
  * `FlairIsoformNode` — 롱리드 Nanopore/PacBio 전사체 식별 및 대안적 스플라이싱
  * `IsoToolsSplicingNode` — 롱리드 전사체 스플라이싱 이벤트 및 엑손 스키핑 분석

---

### 2. [Kingdom] Genomics & Variant Discovery

* **2.1 Short-Read Alignment & Indexing**
  * `VariantInputValidatorNode` — 유전체 변이 분석 입력 자산(VCF, BAM, FASTA) 검증
  * `BwaMem2IndexNode` — 2배 빠른 BWA-MEM2 유전체 참조 서열 인덱싱
  * `BwaMem2AlignNode` — 페어드/싱글 숏리드 전장 유전체 고속 정렬 (SAM/BAM)
  * `Bowtie2AlignNode` — 후성유전체 및 미생물 서열 정렬
  * `MarkDuplicatesNode` — Picard 기반 PCR 및 광학 중복 리드 마킹/제거

* **2.2 Variant Calling & Filtering**
  * `BcftoolsCallNode` — 다배체/이배체 유전체 SNV 및 Small InDel 변이 검출
  * `BcftoolsFilterNode` — 품질 점수(QUAL), 깊이(DP), 가닥 편향(Strand Bias) 필터링
  * `DeepVariantCallNode` — 심층 합성곱 신경망(CNN) 기반 고정밀 변이 호출
  * `Cyvcf2VariantNode` — 초고속 C-바인딩 VCF 파싱 및 커스텀 필터링

* **2.3 Long-Read Genomics & Structural Variants**
  * `MappyAlignNode` — Minimap2 파이썬 바인딩 고속 롱리드 정렬
  * `CuteSvNode` — PacBio HiFi / Nanopore 긴 리드 기반 구조적 변이(SV) 검출
  * `Sniffles2SvNode` — 복합 유전체 구조 변이, 전위, 역위 탐지
  * `MedakaConsensusNode` — Oxford Nanopore 컨센서스 서열 보정 및 변이 호출

* **2.4 High-Throughput Interval Operations**
  * `PysamAnalysisNode` — SAMtools C-API 기반 BAM/SAM/CRAM 통계 산출
  * `PybedtoolsIntervalNode` — BEDTools 연동 게놈 구간 교집합, 병합, 거리 계산
  * `MosdepthCoverageNode` — 유전체 전장 및 엑솜 타겟 영역 커버리지 심도 고속 측정
  * `SeqKitToolNode` — 대용량 FASTA/FASTQ 통계, 서브샘플링, 역상보 변환

---

### 3. [Kingdom] Epigenomics & 3D Chromatin Architecture

* **3.1 Chromatin Accessibility & ATAC-Seq**
  * `AtacInputValidatorNode` — ATAC-seq 입력 데이터 및 중복 제거 유효성 검사
  * `AtacFastpTrimNode` — ATAC-seq 어댑터 트리밍
  * `AtacBwaMem2IndexNode` & `AtacBwaMem2AlignNode` — ATAC-seq 유전체 정렬
  * `AtacMarkDuplicatesNode` & `AtacQualityFilterNode` — 미토콘드리아 및 저품질 리드 필터
  * `Macs3PeakCallingNode` — 개방 크로마틴 좁은 피크(Narrow Peak) 통계적 호출
  * `PyGenrichNode` — ATAC-seq 전용 Tn5 시프트 보정 피크 검출
  * `TobiasFootprintNode` — 전사인자 결합 풋프린팅(Transcription Factor Footprinting)

* **3.2 Protein-DNA Interactions & Histone ChIP**
  * `HomerMotifNode` — 피크 영역 내 신규(De novo) 및 공인 전사인자 결합 모티프 탐색
  * `MemeSuiteNode` — DREME / CentriMo 기반 전사인자 결합 위치 농축 분석
  * `SeacrPeakNode` — CUT&RUN / CUT&Tag 초저노이즈 엄격 피크 호출
  * `DeepToolsProfileNode` — 전사 시작 부위(TSS) 및 유전자 본체 주변 커버리지 매트릭스 계산

* **3.3 3D Genome Architecture**
  * `CoolerMatrixNode` — Hi-C / Micro-C 염색질 접촉 행렬 빌드 및 다중 해상도 저장
  * `CooltoolsTadNode` — 위상 결합 도메인(TAD) 경계 및 구획(Compartment A/B) 탐지
  * `ChromosightLoopNode` — 컴퓨터 비전 기반 크로마틴 루프(Chromatin Loops) 검출

* **3.4 Epigenetic Modifications & Functional Screens**
  * `MethyldackelNode` — WGBS / RRBS 바이설파이트 시퀀싱 메틸화(5mC) 수준 산출
  * `Crispresso2ScreenNode` — CRISPR-Cas9/Cas12 유전자 편집 표적 및 비표적(Off-target) 정량
  * `MageckScreenNode` — 전장 유전체 CRISPR 기능 상실 스크리닝 유전자 랭킹

---

### 4. [Kingdom] Metagenomics, Microbiome & Virology

* **4.1 Taxonomic Profiling**
  * `MetagenomeInputValidatorNode` & `MetagenomeFastpTrimNode` — 메타게놈 입력 전처리
  * `Kraken2ClassifyNode` — $k$-mer 정확 매칭 기반 미생물 분류학적 동정
  * `BrackenAbundanceNode` — 베이지안 추론 기반 종/속 수준 미생물 상대 풍부도 재추정
  * `Metaphlan4ProfileNode` — 마커 유전자 기반 종 수준 정밀 프로파일링

* **4.2 Amplicon & Diversity**
  * `Dada2AmpliconNode` — 16S/18S/ITS 시퀀싱 노이즈 제거 및 ASV(Amplicon Sequence Variant) 도출
  * `BiomFormatNode` — BIOM 형식 미생물 카운트 테이블 변환 및 메타데이터 통합
  * `FastUniFracDistanceNode` — 계통수 기반 UniFrac (Weighted / Unweighted) 베타 다양성 거리 계산
  * `ScikitBioDiversityNode` — Shannon, Simpson, Chao1 알파 다양성 및 Bray-Curtis 거리 계산

* **4.3 Functional Metagenomics**
  * `Humann3PathwayNode` — 메타유전체/메타전사체 유전자 계열(UniRef) 및 MetaCyc 대사 경로 정량

* **4.4 Virome & Phage Discovery**
  * `VirSorter2Node` — 바이러스 및 박테리오파지 서열 탐지 및 선별
  * `CheckVQualityNode` — 바이러스 게놈 완전성, 오염도 및 프로파지 경계 평가
  * `GeNomadViromeNode` — 이동성 유전 인자 및 플라스미드/바이러스 동정

---

### 5. [Kingdom] De Novo Genome Assembly & Annotation

* **5.1 Microbial & Viral Assembly**
  * `AssemblyInputValidatorNode` — 게놈 어셈블리 입력 서열 유효성 검사
  * `SpadesAssembleNode` — 미생물 및 바이러스 isolate 서열 De Novo 조립
  * `FlyeAssembleNode` — 긴 리드(PacBio/ONT) 기반 단일 컨티그 De Novo 조립

* **5.2 Long-Read & Hybrid Assembly**
  * `HifiasmAssembleNode` — PacBio HiFi 기반 진핵생물 일배체형 분리 고품질 유전체 조립
  * `RaconPolishNode` — 서열 정렬 기반 어셈블리 컨티그 연마(Polishing)

* **5.3 Assembly QC & Graph Analytics**
  * `QuastQcNode` — 어셈블리 N50, L50, 오조립(Misassembly), Genome Fraction 평가
  * `Ete3TreeParserNode` — Newick/Nexus 계통수 조작 및 가시화 파싱

* **5.4 Genome Structural & Functional Annotation**
  * `ProkkaAnnotationNode` — 박테리아/고균 유전체 고속 구조 및 기능 주석화
  * `BaktaAnnotationNode` — 최신 NCBI/Uni保 표준 기반 미생물 게놈 정밀 주석화

---

### 6. [Kingdom] Proteomics & Metabolomics

* **6.1 Bottom-Up Proteomics**
  * `MsconvertConvertNode` — 원시 질량분석 바이너리(.raw)의 개방형 mzML 포맷 고속 변환
  * `MsfraggerSearchNode` — 초고속 오픈 질량 검색 기반 펩타이드/단백질 식별
  * `MaxQuantQuantNode` — DDA 라벨프리(LFQ) 및 TMT 다중 단백질 정량화
  * `DiaNNQuantNode` — 신경망 기반 고감도 DIA 펩타이드/단백질 정량
  * `PerseusCliNode` — 단백체 통계 검정, 결측치 대치 및 다변량 분석
  * `MhcquantNeoantigenNode` — 면역 펩티도믹스 및 종양 신생항원(Neoantigen) 예측

* **6.2 Mass Spectrometry Feature Engineering**
  * `PyOpenMSFeatureNode` — LC-MS 피크 피킹, 질량 정렬 및 크로마토그램 추출
  * `PyteomicsMSNode` — 단백질 서열 이론적 질량 산출 및 mzML 파싱
  * `MsDeisotopeNode` — 동위원소 디컨볼루션 및 전하 상태 분리
  * `MassqlQueryNode` — 질량 스펙트럼 SQL 패턴 쿼리 및 프래그먼트 검색

* **6.3 Untargeted & Lipid Metabolomics**
  * `MetaboAnalystRNode` — 대사체 경로 농축 분석 및 바이오마커 발굴
  * `MatchmsSpectrumNode` — 탠덤 질량 스펙트럼(MS/MS) 유사도 계산 및 라이브러리 매칭
  * `Spec2VecEmbeddingNode` — 단어 임베딩 기반 MS/MS 스펙트럼 구조 유사도 산출
  * `SiriusStructureNode` — 동위원소/MS2 기반 미지의 대사물질 분자식 및 화학 구조 규명
  * `MsDialLipidNode` — 지질체학(Lipidomics) 전용 식별 및 정량

---

### 7. [Kingdom] Structural Biology & CADD / Drug Discovery

* **7.1 AI 3D Structure Prediction**
  * `ESMFoldNode` — 대규모 단백질 언어 모델 기반 초고속 3차 구조 단일 서열 예측
  * `ColabFoldAlphaFoldNode` — 다중 서열 정렬(MSA) 기반 AlphaFold2 3차 구조 정밀 예측
  * `OpenFoldNode` — 파이토치 기반 완전 오픈소스 단백질 구조 생성 모델

* **7.2 Molecular Docking & Pocket Finding**
  * `AutoDockVinaDockingNode` — 수용체-리간드 결합 친화도 및 도킹 포즈 탐색
  * `SminaDockingNode` — 에너지 커스텀 스코어링 기반 분자 도킹
  * `GNINADockingNode` — 심층 신경망(Deep Learning) 기반 분자 도킹 및 가상 스크리닝
  * `DiffDockPredictNode` — 확산 모델(Diffusion Model) 기반 분자 결합 포즈 생성
  * `P2RankPocketNode` — 머신러닝 기반 단백질 리간드 바인딩 포켓 예측
  * `FPocketNode` — 보로노이 테셀레이션 기반 단백질 공동(Cavity) 탐색

* **7.3 Molecular Dynamics & Free Energy**
  * `OpenMMSimulationNode` — GPU 가속 생체분자 분자동역학(MD) 시뮬레이션
  * `MDAnalysisNode` — MD 궤적(Trajectory) 원자간 거리, 회전반경, RMSD/RMSF 계산
  * `MDTrajNode` — 2차 구조 전이 및 형태 공간 분석
  * `ProDyDynamicsNode` — 단백질 거대 구조 변형 및 정규 모드 분석(NMA)
  * `PlumedNode` — 메타다이내믹스 기반 자유에너지 지형(Free Energy Landscape) 계산
  * `PmxFEPNode` — 알케미컬 자유에너지 섭동(FEP) 돌연변이 결합 친화도 계산

* **7.4 Cheminformatics & Molecular ML**
  * `RDKitCheminformaticsNode` — 화합물 SMILES 처리, 분자량/LogP/TPSA/Lipinski 지표 산출
  * `OpenBabelConvertNode` — 100여 종의 화학 파일 포맷(SDF, MOL2, PDB, SMILES) 상호 변환
  * `TorchDrugNode` — 그래프 신경망(GNN) 기반 화합물 특성(ADMET) 및 생물활성 예측

---

### 8. [Kingdom] BioPython & Sequence Utilities

* **8.1 Sequence I/O & Feature Translation**
  * `BiopythonSeqIONode` — FASTA, GenBank, EMBL 서열 읽기/쓰기 및 필터링
  * `BiopythonSeqTranscribeNode` — DNA/RNA 전사, 번역, 개방형 판독틀(ORF) 탐색
  * `BioSeqAnalysisNode` — GC 함량, 분자량, 이소스팟 포인트 산출
  * `PyfaidxIndexNode` — 대용량 FASTA 서열의 특정 염기 구간 고속 추출
  * `PyfastxIndexNode` — 대용량 FASTQ 서열 통계 및 서열 파싱

* **8.2 Homology, Motif & Domain Search**
  * `BiopythonBlastNode` — NCBI 로컬/원격 BLAST (blastn, blastp, blastx) 검색
  * `PyhmmerSearchNode` — Pfam 단백질 도메인 HMMER 모델 초고속 검색
  * `BiopythonBioMotifsNode` — 위치 특이적 가중치 행렬(PWM) 및 서열 로고 분석
  * `PyTFBSScanNode` — 전사 조절 영역 내 전사인자 결합 부위 스캔
  * `EdlibAlignNode` — 레벤슈타인 편집 거리 기반 초고속 전역/반전역 서열 정렬

* **8.3 Primer & Molecular Cloning Design**
  * `Primer3DesignNode` — PCR 프라이머 및 시퀀싱 올리고 자동 설계
  * `BiopythonRestrictionNode` — 제한효소 절단 부위 매핑 및 클로닝 시뮬레이션
  * `BiopythonProtParamNode` — 단백질 불안정성 지수, 친수성, 흡광 계수 산출
  * `CodonWAnalysisNode` — 코돈 사용 편향성(Codon Usage Bias) 및 RSCU 계산

* **8.4 Comparative Genomics & Structural Alignment**
  * `BiopythonAlignIONode` — Clustal, Phylip, Stockholm 다중 정렬 포맷 변환
  * `BiotiteStructureAlignNode` — 단백질 3차원 입체 구조 간 중첩 및 C-alpha RMSD 정렬
  * `BiopythonBioPDBNode` — PDB 좌표 파싱 및 원자간 접촉 잔기 분석
  * `BioKEGGPathwayNode` — KEGG 대사 경로 매핑 및 유전자 네트워크 시각화
  * `BiopythonPhyloNode` — 계통수 Newick 트리 파싱 및 계통 거리 계산
  * `DnaFeaturesViewerNode` — 플라스미드 및 선형 유전자 서열 구조 주석 도판 생성
  * `HelicalWheelNode` — 양친매성 알파 나선(Helical Wheel) 2차원 투영 다이어그램
  * `LogomakerVisualizerNode` — DNA/RNA/단백질 서열 모티프 시퀀스 로고 렌더링
  * `PyCircosPlotNode` — 유전체 전장 변이/구조를 보여주는 Circos 원형 다이어그램

---

### 9. [Kingdom] Tier-1 Publication-Grade Visualizers (300+ DPI Native IMAGE)

* **9.1 Expression & Statistical Analytics**
  * `VolcanoPlotVisualizerNode` — 차등 발현 유전자 $\log_2\text{FC}$ vs $-\log_{10}(p)$ 화산 플롯 (상위 유전자 자동 라벨링)
  * `ClustermapHeatmapVisualizerNode` — 계층적 군집화(Hierarchical Clustering) 발현 히트맵 및 덴드로그램
  * `MaPlotVisualizerNode` — 발현 평균 강도(A) 대비 발현 변화율(M) MA 플롯
  * `GseaEnrichmentPlotVisualizerNode` — GSEA 시그니처 러닝 농축 점수(Enrichment Score) 곡선 도판
  * `KaplanMeierSurvivalVisualizerNode` — 환자 유전자 발현군별 생존율 Kaplan-Meier 곡선 및 Log-rank $p$

* **9.2 Single-Cell & Spatial Cell-Fate**
  * `UmapScatterVisualizerNode` — 단일세포 2차원/3차원 UMAP 차원축소 산점도
  * `SpatialTissueOverlayVisualizerNode` — H&E 조직 염색 이미지 위 유전자 발현 핫스팟 중첩 렌더링
  * `SankeyCellFateVisualizerNode` — 단일세포 분화 궤적 및 세포 운명 전환 샌키(Sankey) 다이어그램

* **9.3 Genomic Variation & Association Tracks**
  * `ManhattanPlotVisualizerNode` — GWAS 및 전장 유전체 변이 연관성 맨해튼 플롯
  * `QqPlotVisualizerNode` — 이론적 분포 대비 관측치 분위수-분위수(Q-Q) 플롯
  * `OncoPrintVisualizerNode` — 암 유전체 샘플별 돌연변이(SNV/InDel/CNV) OncoPrint 매트릭스
  * `LinkageDisequilibriumVisualizerNode` — 유전체 변이 간 연쇄불평형(LD) $r^2$ 히트맵
  * `ChipAtacCoverageProfileVisualizerNode` — 고해상도 게놈 브라우저 스타일 피크 커버리지 트랙 및 유전자 주석

* **9.4 Microbiome & Phylogenetics**
  * `MicrobiomeStackedBarVisualizerNode` — 샘플별 분류군 상대 풍부도 누적 막대 그래프
  * `PcoaScatterVisualizerNode` — Bray-Curtis / UniFrac 거리 기반 PCoA 2D/3D 좌표 서열화 산점도
  * `PhylogeneticTreeVisualizerNode` — 출판 규격 원형/직사각형 계통수(Phylogenetic Tree) 도판
  * `SyntenyGenomeVisualizerNode` — 이종 간 게놈 공선성(Synteny) 및 역위/전위 블록 리본 플롯

* **9.5 Structural Biology & CADD Visualizers**
  * `ProteinLigandInteractionVisualizerNode` — 수용체-리간드 수소결합, 소수성 상호작용 2D 다이어그램
  * `RamachandranPlotVisualizerNode` — 단백질 주쇄 이면각($\phi, \psi$) 라마찬드란 플롯 및 선호 영역 평가
  * `MdTrajectoryPlotterVisualizerNode` — 분자동역학 RMSD, RMSF, 회전반경 궤적 시계열 도판

---

## 📊 도메인별 노드 통계 요약 (Taxonomic Census)

| 계통 분류 (Kingdom) | 대표 핵심 도메인 | 구현 노드 수 | 주요 지원 기능 |
| :--- | :--- | :---: | :--- |
| **1. Transcriptomics** | Bulk, scRNA, Spatial, Isoform | **25개** | 발현 정량, DEG, 클러스터링, 공간 매핑, RNA 벨로시티 |
| **2. Genomics & Variants** | WGS/WES, 롱리드, SV, BAM 연산 | **16개** | 정렬, 중복 마킹, 변이 호출, 롱리드 SV, 커버리지 측정 |
| **3. Epigenomics & 3D** | ATAC-Seq, ChIP, 3D Hi-C, 메틸화 | **17개** | 피크 호출, 풋프린팅, TAD 탐지, 메틸화, CRISPR 스크린 |
| **4. Metagenomics** | 분류, 앰플리콘, 대사 경로, 바이러스 | **14개** | Kraken2 분류, ASV 도출, 대사 경로, 바이러스 발굴 |
| **5. Genome Assembly** | Isolate 조립, 롱리드 조립, 주석화 | **8개** | SPAdes, Hifiasm, QUAST 평가, Bakta/Prokka 주석화 |
| **6. Proteomics & Metabolomics**| DDA/DIA 질량분석, 대사체, 지질체 | **16개** | MaxQuant, DIA-NN, 분자 네트워킹, 지질체 정량 |
| **7. Structural Biology & CADD**| 3D 예측, 도킹, MD, 화합물 ML | **18개** | AlphaFold2, ESMFold, Vina 도킹, OpenMM 시뮬레이션 |
| **8. BioPython & Utilities** | 서열 조작, 상동성, 프라이머, 비교 | **23개** | BLAST, HMMER, 프라이머 설계, Circos, 플라스미드 |
| **9. Tier-1 Visualizers** | 300+ DPI 전문 출판 도판 렌더러 | **20개** | Volcano, Heatmap, UMAP, Manhattan, Coverage, OncoPrint |
| **총계 (Total Catalog)** | **8대 핵심 생명과학 분야 + 출판 시각화** | **157개** | **전주기 생물정보학 엔드투엔드 파이프라인 완결** |
