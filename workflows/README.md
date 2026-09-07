# ComfyBIOWMS Example Workflows

Ready-to-use ComfyUI workflow JSON files for various bioinformatics and computational biology pipelines.

## Available Workflows

### 5대 표준 생물학 골든 벤치마크 워크플로우 (논문 3절 검증 파이프라인)

1. **De novo 유전체 조립 (`01_genome_assembly_phix174.json`)**
   - **표준 도구 체인**: `Fastp v0.24.0` ➔ `SPAdes v4.3.0` ➔ `QUAST v5.3.0`
   - **노드 구성**: `SampleMetadataValidatorNode` ➔ `FastpTrimNode` ➔ `SpadesAssembleNode` ➔ `QuastQcNode`
   - **주요 산출물**: 트리밍 FASTQ, `contigs.fasta` 조립 서열, QUAST 품질 평가 지표 리포트 (`report.tsv`)

2. **후성유전체 ATAC-Seq 오픈 크로마틴 분석 (`02_epigenomics_atacseq_macs3.json`)**
   - **표준 도구 체인**: `Fastp v0.24.0` ➔ `BWA-MEM2 v2.2.1` ➔ `Samtools v1.24` ➔ `MACS3 v3.0.4`
   - **노드 구성**: `SampleMetadataValidatorNode` ➔ `FastpTrimNode` ➔ `BwaMem2IndexNode` ➔ `BwaMem2AlignNode` ➔ `MarkDuplicatesNode` ➔ `AtacQualityFilterNode` ➔ `Macs3PeakCallingNode`
   - **주요 산출물**: 정렬/중복제거/품질필터링 BAM, MACS3 피크 콜링 결과 (`narrowPeak`, `summits.bed`)

3. **DNA 생식세포 변이 탐지 (`03_dna_variant_calling_bcftools.json`)**
   - **표준 도구 체인**: `Fastp v0.24.0` ➔ `BWA-MEM2 v2.2.1` ➔ `Samtools v1.24` ➔ `BCFtools v1.22`
   - **노드 구성**: `SampleMetadataValidatorNode` ➔ `BwaMem2IndexNode` ➔ `BwaMem2AlignNode` ➔ `MarkDuplicatesNode` ➔ `BcftoolsCallNode` ➔ `BcftoolsFilterNode`
   - **주요 산출물**: 정렬 BAM, 원천 BCF/VCF (`raw.vcf`), 품질/깊이 필터링 변이 목록 (`filtered.vcf`)

4. **벌크 전사체 Bulk RNA-Seq 발현차이분석 (`04_bulk_rnaseq_deseq2.json`)**
   - **표준 도구 체인**: `Fastp v0.24.0` ➔ `Salmon v1.10.3` ➔ `tximport v1.34.0` ➔ `DESeq2 v1.46.0`
   - **노드 구성**: `SampleMetadataValidatorNode` ➔ `FastpTrimNode` ➔ `SalmonIndexNode` ➔ `SalmonQuantNode` ➔ `TximportNode` ➔ `DESeq2AnalysisNode`
   - **주요 산출물**: 트랜스크립톰 정량(`quant.sf`), 유전자 카운트 행렬(`count_matrix.csv`), DESeq2 DEG 통계표(`deseq2_results.csv`)

5. **메타게놈 균총 분류 (`05_metagenomics_kraken2_bracken.json`)**
   - **표준 도구 체인**: `Fastp v0.24.0` ➔ `Kraken2 v2.1.3` ➔ `Bracken v3.0`
   - **노드 구성**: `SampleMetadataValidatorNode` ➔ `FastpTrimNode` ➔ `Kraken2ClassifyNode` ➔ `BrackenAbundanceNode`
   - **주요 산출물**: 정제 FASTQ, Kraken2 분류 리포트(`kraken2_report.txt`), Bracken 종 수준 상대풍부도 재추정치(`bracken_report.txt`)

---

### 기타 시각화 워크플로우

- **Publication Visualizers (`03_publication_visualizers.json`)**:
  - 독립 시각화 노드 모음: Volcano Plot, Clustermap Heatmap, UMAP Scatter Plot.
  - 300+ DPI 출판 품질 그래픽 실시간 캔버스 렌더링.

## How to Load in ComfyUI
1. 웹 브라우저에서 ComfyUI를 실행합니다 (`http://127.0.0.1:8188`).
2. `workflows/` 폴더 내 원하는 `.json` 파일을 ComfyUI 캔버스 위로 드래그 앤 드롭하거나, 사이드 패널의 **Load** 버튼으로 불러옵니다.
3. 노드 위젯 필드의 입력 데이터 디렉터리 경로를 확인하거나 필요에 따라 조정합니다.
4. **Queue Prompt** 버튼을 클릭하여 파이프라인을 실행합니다.
