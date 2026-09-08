# ComfyBIOWMS

ComfyUI 기반 고신뢰도 생물정보학 워크플로우 관리 확장 시스템 (Bioinformatics Workflow Management System Extension for ComfyUI).

ComfyBIOWMS는 검증된 오픈소스 비순환 방향 그래프(DAG) 인터페이스인 ComfyUI를 생물정보학 파이프라인 제어 엔진으로 확장하여, 대용량 멀티오믹스 데이터의 전처리, 정렬, 정량, 유전체 변이 분석, 차등 발현(DEG), 메타게놈 분류 및 출판급 시각화를 시각적 캔버스 상에서 직관적으로 구성하고 재현 가능하게 실행합니다.

---

## 1. 핵심 아키텍처 및 특징

- **엄격한 과학적 무결성**: 임의의 합성 데이터(synthetic fallback) 및 난수 기반 플레이스홀더를 전면 배제하며, 모든 분석 및 시각화는 실제 생물학 데이터와 공식 도구 바이너리에 기반합니다.
- **2개 노드 실행 계층 (총 98개 활성 등록 노드)**:
  - **Class 1 (40개 노드)**: ComfyUI 인메모리 고속 과학 연산(Biopython, Scanpy, AnnData, tximport) 및 20종의 300+ DPI 출판급 시각화 노드(`IMAGE` 텐서 직결).
  - **Class 2 (58개 노드)**: 격리 Conda 환경 기반 고성능 외부 CLI 브리지 노드(fastp, Salmon, BWA-MEM2, samtools, bcftools, SPAdes, QUAST, Kraken2, Bracken, Bioconductor DESeq2 등).
- **입력 반응형 캐시 무효화 (`IS_CHANGED`)**: 파일 경로, 파일 크기(`st_size`), 최종 수정 시각(`st_mtime_ns`) 및 파라미터 해싱 기반의 결정론적 캐시 키를 생성하여, 파라미터 또는 파일 변경 시 영향받는 서브그래프만 선택적으로 재실행합니다 (캐시 적중 시 **0.15 ms** 이내 즉각 반환).
- **자동 감사 추적 (Audit Trail)**: 모든 외부 도구 실행 시 작업 디렉터리에 `run_manifest.json`, `run_manifest.sh`, `stdout.log`, `stderr.log`를 자동 기록하여 완전한 독립 재실행성을 보장합니다.
- **초경량 런타임 오버헤드**: 네이티브 바이너리 직접 실행 대비 래퍼 오버헤드는 단 **0.0023초 (1.0%)** 수준이며, 파일 경로(`STRING`) 포트 전달 아키텍처로 백엔드 파이썬 프로세스 메모리를 안전하게 보존합니다 (Peak RSS: **136 MB**).

---

## 2. 노드 분류 및 검증 현황

전체 98개 등록 노드의 분류, 실행 환경, 검증 수준 및 입출력 명세는 다음 문서에서 투명하게 공개되어 있습니다:
- [노드 검증 매트릭스 (Node Verification Matrix)](docs/node-verification-matrix.md)
- [보충표 S1: 전체 노드 카탈로그 (Supplementary Table S1)](supplementary_table_s1_node_catalog.md)
- [사례 연구 실험 프로토콜 (Case Study Protocols)](docs/case-study-protocols.md)
- [기준선 도구 및 nf-core 동등성 비교 보고서](results/validation/baseline_concordance_report.md)

---

## 3. 설치 및 환경 구성

### 3.1 기본 요구조건
- macOS (Apple Silicon / Intel) 또는 Linux (x86_64)
- Python 3.10+ 및 PyTorch
- Conda (Miniconda 또는 Anaconda)

### 3.2 Conda 환경 활성화
ComfyBIOWMS는 도메인별 의존성 충돌을 원천 차단하기 위해 독립된 전용 Conda 환경을 자동으로 감지하고 호출합니다:
- `bulk_rna_seq`: fastp 1.3.6, salmon 2.5.1, bioconductor-deseq2 1.46.0, bioconductor-tximport 1.38.2
- `variant_analysis`: bwa-mem2, samtools 1.24, bcftools 1.24
- `genome_assembly`: spades 4.3.0, quast 5.3.0, bwa, fastp
- `metagenome`: kraken2 2.17.1, bracken 2.9, fastp
- `epigenomics`: macs3, samtools, bedtools, fastp

---

## 4. 자동화 테스트 스위트 검증

저널 투고 기준의 엄격한 재현성을 입증하기 위해 계층화된 pytest 검증 스위트를 제공합니다:

```bash
# 1. 실행 계약, 환경 격리, 캐시 무효화 및 산출물 무결성 회귀 검증 (7 tests)
pytest -v tests/test_execution_contract.py

# 2. 원 도구(Native CLI / Bioconductor)와의 수치 동등성 E2E 검증 (3 tests)
pytest -v -m e2e tests/test_native_concordance.py

# 3. ComfyUI 노드 등록, 포트 서명 및 그래프 무결성 통합 검증 (3 tests)
pytest -v tests/test_comfyui_integration.py

# 4. 단독 모듈 및 배치 샘플 처리 검증 (30 tests)
pytest -v tests/test_standalone_*.py tests/test_final_case_study_workflows.py

# 5. 전체 43개 테스트 전수 실행
pytest -v
```

---

## 5. 런타임 오버헤드 벤치마크 실행

네이티브 CLI와 ComfyBIOWMS 노드 및 캐시 적중 속도를 직접 측정하려면 다음 스크립트를 실행합니다:

```bash
python3 engine/scripts/measure_overhead.py
```

결과는 `results/overhead_benchmark/overhead_summary.json` 및 마크다운 리포트로 저장됩니다.

---

## 6. 라이선스

본 소프트웨어는 [MIT License](LICENSE)에 따라 자유롭게 사용 및 수정할 수 있습니다.
