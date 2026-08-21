# ComfyBioWMS: Bioinformatics Workflow & Analysis Platform Manuscript Workspace

이 작업 공간(`paper_platform/`)은 **순수 생물정보학 워크플로우 및 시각 분석 플랫폼(Bioinformatics Workflow and Analysis Platform)**으로서의 ComfyBioWMS 학술 논문을 체계적으로 집필하고 실증하기 위해 구축된 독립 작업 공간입니다.

---

## 📁 디렉토리 구조 및 핵심 문서 안내

1. **[manuscript.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/manuscript.md)**
   - **(전면 개정 완료)** LLM/TSR 및 인위적 레이턴시 비교를 완전히 걷어내고, **반응형 노드 DAG 캔버스, 8대 도메인 24개 골드 스탠다드 워크플로우, 20종 내장 출판급 시각화(Tier-1 Visualizers), 소켓 타입 계약, Conda 런타임 격리, 이미지 임베디드 원클릭 복원(One-Drop Replication)**을 중심으로 작성된 완성도 높은 플랫폼 논문 본문 초안.

2. **[node_phylogeny_tree.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/node_phylogeny_tree.md)**
   - **(신규 추가)** 플랫폼에 구현된 **157개 전 도메인 노드**를 9대 대분류(Kingdom) 및 35개 세부 계통(Phylum)으로 체계화한 **생물학적 진화계통수(Node Phylogenetic Tree) 카탈로그**.

3. **[architecture_and_roadmap.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/architecture_and_roadmap.md)**
   - 3계층 아키텍처 다이어그램, 5대 핵심 가치 제안, 4대 실증 검증 계획 및 타겟 저널(*Bioinformatics*, *PLOS Computational Biology*, *GigaScience*) 로드맵.

4. **[domain_workflows_taxonomy.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/domain_workflows_taxonomy.md)**
   - 8대 핵심 생물정보학 도메인 및 24개 세부 워크플로우 명세 (입력, 단계별 도구, 소켓 타입 계약, 출력 도판).

5. **[literature_domain_scope_survey.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/literature_domain_scope_survey.md)**
   - *Nature Biotechnology, Bioinformatics, PLOS Comp Bio, BMC Bioinformatics, ChemRxiv* 등 주요 플랫폼 논문들의 **도메인 범위, 노드/도구 수, 실증 사례 연구 범위(3-Tier Scope)** 심층 분석 보고서.

6. **[literature_review.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/literature_review.md)**
   - 기존 생물정보학 WMS 및 GUI 플랫폼(Nextflow, Snakemake, Galaxy, KNIME, BioImageIT/VERA) 동향 분석 및 기능별 정밀 비교표(Feature Matrix).

7. **[benchmark_datasets_audit.md](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/benchmark_datasets_audit.md)**
   - **(신규 추가)** 전 세계 주요 학술지 및 WMS 플랫폼에서 파이프라인 검증에 사용하는 **국제 공인 골드 스탠다드 벤치마크 데이터셋(FDA SEQC, NCBI RefSeq PhiX174, ENCODE Tier-1, ZymoBIOMICS, NIST GIAB)** 전수 조사 및 공식 기탁 번호(Accession) 분석 보고서.

8. **[benchmarks/](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/benchmarks/)**
   - **(신규 구축)** 실제 공인 벤치마크 데이터셋(15개 파일, ~170MB) 다운로드 및 엔드투엔드 파이프라인 실측치 자동 추출 실행 환경.
   - **`download_benchmark_datasets.py`**: 공인 저장소 원본 데이터 자동 다운로더.
   - **`run_real_benchmarks.py`**: 4대 공인 벤치마크 실제 실행 및 실측치 로그(`benchmark_results.json`) 생성기.

9. **[figures/](file:///Users/ydj_mac/Library/CloudStorage/OneDrive-개인/Project/ComfyBioWMS/paper_platform/figures/)**
   - **`figure_1_system_architecture.png`**: (Figure 1) 300 DPI 출판 규격의 3계층 시스템 아키텍처 및 전주기 데이터 흐름(Dataflow) 배선도.
   - **`figure_2_node_phylogeny_tree.png`**: (Figure 2) 300 DPI 출판 규격의 157개 노드 진화계통수, 노드 센서스, 및 7대 소켓 계약 아키텍처 도판.
   - **`figure_3_case_studies.png`**: (Figure 3) 4대 공인 벤치마크(SEQC 전사체, PhiX174 조립, ENCODE ATAC-Seq, 메타게놈) 엔드투엔드 생물학적 실증 도판.
   - **파이썬 렌더링 스크립트**: `generate_architecture_figure.py`, `generate_phylogeny_figure.py`, `generate_case_studies_figure.py`
