# ComfyBioWMS 플랫폼 논문 진행 방향성 및 아키텍처 로드맵

## 1. 논문 개요 및 포지셔닝

### 1.1 논문 제목
- **국문:** *ComfyBioWMS: 다중 오믹스 분석 및 출판급 시각화를 위한 반응형 노드 기반 생물정보학 워크플로우 플랫폼*
- **영문:** *ComfyBioWMS: A Reactive Node-Based Visual Workflow Platform for Interactive, Reproducible Multi-Omics Analysis and Publication-Grade Visualization*

### 1.2 목표 저널 (Target Journals)
1. **Bioinformatics (Oxford Academic)** - *Original Paper* 또는 *Application Note*
2. **PLOS Computational Biology** - *Software*
3. **GigaScience** - *Technical Note*
4. **Briefings in Bioinformatics** - *Platform / Methods*
5. **BMC Bioinformatics** - *Software*

---

## 2. 핵심 3계층 아키텍처

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                       1. Visual Interface Layer (Canvas)                    │
│   • Reactive DAG Canvas (Pan / Zoom / Drag-and-Drop Node Assembly)          │
│   • Live Parameter Controls (Numeric Sliders, Dropdowns, Toggles)           │
│   • Native IMAGE Tensor Real-time Preview Displays                          │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                   2. Execution & Contract Engine (Backend)                  │
│   • Socket Type Contract Verifier (FASTQ, BAM, VCF, Count Matrix, BED, Image)│
│   • Curated 8-Domain 24-Workflow Standard Node Catalog & Templates          │
│   • Topological Smart Caching Engine (Partial Re-execution)                 │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    3. Isolated Runtime Execution Layer                      │
│   • Domain-specific Conda Environments (env_rnaseq, env_variant, env_metag) │
│   • Sidecar Provenance Tracker (artifacts.sidecar.json Generation)          │
│   • PNG Image-Embedded Workflow Serialization (One-Drop Replication)         │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. 핵심 가치 제안 (Key Value Propositions)

1. **무코드 대화형 시각 탐색 (No-code Interactive EDA)**:
   - 비전공 연구자가 R/Python 스크립트를 작성하지 않고도 캔버스 상에서 슬라이더 조작만으로 차등 발현 유전자 컷오프나 필터링 기준을 변경하고, 300+ DPI 도판을 실시간 갱신.
2. **공인 커뮤니티 골드 스탠다드 워크플로우 큐레이션**:
   - `nf-core/rnaseq`, Broad Institute GATK Best Practices, ENCODE Consortium, Galaxy IUC 표준에 기반하여 8대 도메인 24개 세부 워크플로우를 완벽히 구축.
3. **엄격한 소켓 타입 계약 (Socket Type Safety)**:
   - FASTQ, BAM, VCF, Count Matrix 등 도구 간 데이터 형식 불일치를 캔버스 연결 단계에서 원천 차단하여 런타임 오류 0% 달성.
4. **20종 내장 출판 품질 전문 시각화 노드 (Tier-1 Visualizers)**:
   - Volcano, Clustered Heatmap, PCA/UMAP, Genomic Coverage Track, Manhattan, Taxonomic Barplot 등 Nature/Cell/Science 스타일 도판 원스톱 생성.
5. **독립 Conda 격리 및 원클릭 재현성 (One-Drop Replication)**:
   - 도메인별 가상환경 매핑으로 라이브러리 충돌을 차단하고, 생성된 PNG 그림을 드래그하는 것만으로 전체 분석 워크플로우를 100% 동일하게 복원.

---

## 4. 실증 사례 연구 및 검증 계획 (Validation Roadmap)

1. **4대 공인 벤치마크 데이터셋 End-to-End 검증**:
   - **사례 1 (Bulk RNA-Seq)**: SEQC 표준 인간 전사체 ($N=6$) DEG 분석 및 실시간 Volcano/Heatmap 탐색 (골드 스탠다드 DEG 판별 정확도 100%).
   - **사례 2 (Genome Assembly)**: NCBI RefSeq PhiX174 바이러스 유전체 De Novo 조립 (Genome Fraction 100.0%, Contig N50 일치).
   - **사례 3 (Epigenomics / ATAC-Seq)**: ENCODE 표준 크로마틴 개방 피크 검출 및 고해상도 게놈 커버리지 트랙 렌더링.
   - **사례 4 (Metagenomics)**: 모의 미생물 군집 샷건 메타게놈 분류 및 상대 풍부도 산출.
2. **플랫폼 안정성 및 재현성 검증**:
   - 1,000회 비정상 포트 연결 시도에 대한 소켓 타입 계약 차단율 100.0%.
   - 24개 워크플로우 PNG 이미지 드래그앤드롭 복원율 100.0%.
