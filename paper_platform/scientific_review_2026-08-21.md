# ComfyBioflow manuscript scientific review

검토 대상: `paper_platform/manuscript.md`  
검토일: 2026-08-21  
판정: **Major revision — 실제 실행 전에는 결과 논문으로 제출 불가**

## 총평

ComfyUI 기반의 반응형 생물정보학 워크플로우 플랫폼이라는 연구 질문은 타당하고 소프트웨어 논문으로 발전할 가능성이 있다. 그러나 현재 원고는 구현 범위, 데이터 출처, 검증 방법 및 결과를 실제 증거보다 훨씬 강하게 서술한다. 특히 `benchmarks/run_real_benchmarks.py`는 Salmon, DESeq2, SPAdes, QUAST, Bowtie2, MACS3, Kraken2 또는 Bracken을 실행하지 않고 미리 정한 수치를 기록한다. Figure 3도 실측 결과가 아니라 난수와 고정 배열로 생성된다. 그러므로 제3장과 초록의 결과 문장은 현 단계에서 모두 검증 계획형 미래시제로 바꾸고, 실제 분석 후 사전 정의한 평가 지표로 교체해야 한다.

## 구조화된 리뷰 항목

아래 `anchor_text`는 모두 원고에 존재하는 정확한 구절이다. 대시보드 등록 시 각 행을 독립 리뷰로 사용해야 한다.

| # | Reviewer | Severity | Section | anchor_text | Comment |
|---:|---|---|---|---|---|
| 1 | Methods Reviewer | major | abstract | `분석 정확도 100%와 완벽한 시각 탐색 편의성을 입증하였다.` | 실제 벤치마크가 아직 수행되지 않았으므로 완료형 결과와 “완벽”이라는 표현은 허용할 수 없다. 현재는 “검증할 예정이다”로 바꾸고 실제 실행 후 사전 정의한 도메인별 지표를 보고하라. |
| 2 | Methods Reviewer | major | results | `실증 평가를 수행하였다` | 제공된 실행기는 실제 도구 체인을 실행하지 않는다. 각 사례에서 원시 입력부터 최종 결과까지 플랫폼 노드가 생성한 명령, 버전, 로그, 산출물 및 실패 상태를 보존하는 진짜 엔드투엔드 실행이 필요하다. |
| 3 | Methods Reviewer | major | results | `4대 공인 벤치마크 데이터셋 기반 엔드투엔드 생물학적 사례 연구` | Figure 3은 `generate_case_studies_figure.py`에서 난수와 고정값으로 만든 도식이며 실험 산출물이 아니다. 실행 전에는 “예상 평가 설계의 개념도”로 표시하고, 실행 후에는 원본 결과 파일에서 재현 가능하게 그린 도판으로 교체하라. |
| 4 | Methods Reviewer | major | methods | `골든 패스(Golden Path) 파이프라인을 1:1로 매핑` | 표준 파이프라인과의 1:1 동등성을 입증할 매핑표가 없다. 각 워크플로우에 대해 참조 파이프라인 릴리스, 도구·버전, 파라미터, 전처리, 참조 유전체/DB, 출력 스키마 및 허용 오차를 명시하고 동일 입력에 대한 concordance를 검증하라. |
| 5 | Methods Reviewer | major | methods | `전체 파이프라인을 100% 동일하게 복원할 수 있다.` | PNG의 그래프·파라미터 메타데이터만으로 입력 데이터, 코드 커밋, Conda 해결 결과, 외부 DB 및 OS 차이를 복원할 수 없다. “그래프 구성을 복원”으로 범위를 좁히거나 입력 해시, 잠금 파일/컨테이너 digest, 코드 버전 및 DB 버전을 함께 캡처하라. |
| 6 | Methods Reviewer | major | methods | `완벽한 계산 재현성을 보장한다.` | 명령·도구 버전·입력 해시·자원 사용량은 provenance에는 유용하지만 계산 재현성을 보장하지 않는다. 출력 해시, 난수 시드, 환경 lock/container digest, 참조 데이터 버전과 서로 다른 깨끗한 호스트에서의 재실행 검증이 필요하다. |
| 7 | Methods Reviewer | major | methods | `` `FASTQ_PAIR`, `SAM_BAM_INDEXED`, `VCF_FILTERED`, `COUNT_MATRIX_CSV`, `PEAK_BED`, `CONTIG_FASTA`, `IMAGE_TENSOR` `` | 현재 `nodes/` 구현에서 이 사용자 정의 반환 타입을 찾을 수 없고 등록 노드는 주로 `STRING`, `IMAGE`, `INT`, `FLOAT`를 사용한다. 타입 계약이 실제 구현될 때까지 이 주장과 Figure 2c를 제거하거나 현 구현의 타입 검사 범위를 정확히 기술하라. |
| 8 | Methods Reviewer | major | methods | `157개의 표준 노드` | `nodes/registry.py`의 정적 등록 키는 현재 189개이므로 원고 및 Figure 1–2의 157개와 일치하지 않는다. 커밋 고정형 자동 센서스 스크립트로 실행 가능 노드, 시각화 노드, placeholder 노드를 구분해 집계하라. |
| 9 | Methods Reviewer | major | methods | `8대 핵심 도메인 24개 세부 워크플로우` | Table 1은 8개의 대표 도구 체인만 제시하고 24개 워크플로우의 입력, 단계, 출력 및 검증 상태를 열거하지 않는다. 24개를 주장하려면 보충표에 모든 워크플로우와 구현·테스트 상태를 공개하라. |
| 10 | Methods Reviewer | major | methods | `수 밀리초 내에 재계산` | 하드웨어, 데이터 크기, warm/cold cache, 반복 횟수, 중앙값·IQR/p95 및 캐시 무효화 정확성 검증이 없다. 변경 유형별로 재실행 노드 집합과 fresh run 대비 출력 동등성을 함께 측정하라. |
| 11 | Methods Reviewer | major | discussion | `비전공자도 즉각적인 시각 피드백을 받으며` | 비전공자 접근성과 편의성은 사용자 연구 없이 주장할 수 없다. 사용자 평가를 하지 않는다면 “목표로 설계되었다”로 낮추고, 수행한다면 사전 정의 과제, 완료율, 시간, 오류, SUS 및 비교 플랫폼을 보고하라. |
| 12 | Statistics Reviewer | major | abstract | `분석 정확도 100%` | RNA-seq 분류, assembly 완성도, ATAC peak concordance, 메타게놈 조성 오차는 서로 다른 척도이므로 하나의 “정확도”로 합칠 수 없다. 도메인별 primary endpoint와 분모를 별도로 보고하라. |
| 13 | Statistics Reviewer | major | results | `100.0%의 민감도와 특이도로 완벽히 식별하였다.` | 선택한 10개 양성·10개 음성 유전자만으로 100%를 주장하면 선택 편향이 크고 불확실성도 상당하다. 독립적으로 정의된 전체 truth set, confusion matrix, ROC/PR 또는 상관·bias 지표와 95% 신뢰구간을 보고하라. |
| 14 | Statistics Reviewer | major | results | `10개 핵심 유의 차등 발현 유전자` | 이 10개 유전자 목록이 SEQC의 공인 truth set이라는 근거가 제시되지 않았다. SEQC의 TaqMan/ERCC 자료에서 분석 전에 포함 기준과 방향성을 고정하고, 원고 작성 후 선택한 유전자 패널을 검증 집합으로 사용하지 말라. |
| 15 | Statistics Reviewer | major | results | `오차 범위 1.2% 이내(일치도 98.8%)` | “일치도”의 수식, taxonomic level, 미검출/오분류 처리, reference DB build 및 반복이 정의되지 않았다. L1/Bray–Curtis 등 사전 정의한 거리, species precision/recall, false positives와 반복·subsampling 분포를 보고하라. |
| 16 | Statistics Reviewer | major | results | `Genome Fraction 100.0%, Contig N50 = 5,386 bp` | 5.4 kb 단일 유전체에서 N50만으로 assembly 정확도를 판단할 수 없고, 현재 실행기는 이 값을 계산하지 않고 상수로 기록한다. QUAST의 genome fraction, misassemblies, mismatches/100 kb, indels/100 kb, total length, circularization 및 실제 명령 로그를 보고하라. |
| 17 | Statistics Reviewer | major | results | `IDR < 0.01` | IDR은 적절한 biological replicate와 명시된 peak-ranking/IDR 절차가 필요하다. replicate 수, read depth, downsampling, QC 임계값, reference peak set 및 peak overlap/concordance 지표를 사전에 정의하라. |
| 18 | Statistics Reviewer | major | discussion | `기능별 심층 비교` | Table 2의 △/○/◎ 표시는 평가 규칙과 근거가 없어 정량 비교가 아니다. 기능 정의, 버전, 배포 방식, 동일 과제 프로토콜 및 출처를 제시하거나 표를 서술적 기능 매트릭스로 낮춰라. |
| 19 | Domain Reviewer | major | results | `` `GSE47792` / `GSE47726` `` | GSE47726은 SEQC가 아니라 류마티스관절염 anti-TNF 반응의 microarray 연구다. SEQC SuperSeries의 정확한 SubSeries와 실제 사용할 GSM/SRR run을 확정하고 샘플·lane 선택 근거를 제시하라. |
| 20 | Domain Reviewer | major | results | `` `SRR11038933` `` | NCBI 메타데이터상 이 run은 PhiX174가 아니라 야생 뒤쥐 장관의 16S amplicon 자료다. 실제 PhiX control read accession을 사용하고 reference와 read organism의 일치를 자동 검증하라. |
| 21 | Domain Reviewer | major | results | `` `ENCSR000EMS` `` | ENCSR000EMS는 GM12878 ATAC-seq가 아니라 GM12865 DNase-seq로 연결된다. 정확한 ATAC-seq experiment 및 biological replicate accession을 다시 선정하라. |
| 22 | Domain Reviewer | major | results | `` `SRR12359623` `` | NCBI 메타데이터상 이 run은 ZymoBIOMICS shotgun metagenome이 아니라 사과의 single-end genotyping-by-sequencing 자료다. 제조사 표준 조성과 정확히 연결되는 Zymo run 및 lot/composition 기준을 사용하라. |
| 23 | Domain Reviewer | major | methods | `국제 공인 표준 컨소시엄` | nf-core/IUC는 큐레이션 생태계, GATK Best Practices와 ENCODE는 도메인별 권고, CASP는 예측 평가 실험, PDB/RefSeq는 데이터 자원이다. 모두를 동일한 “공인 골드 스탠다드 파이프라인”으로 묶지 말고 각 자원의 역할과 동등성 수준을 구분하라. |
| 24 | Domain Reviewer | minor | methods | `계통분류학적 계통수(Phylogeny)` | 소프트웨어 노드 범주는 공통 조상이나 진화 거리를 나타내지 않으므로 phylogeny/cladogram/kingdom/phyla는 오해를 부른다. “hierarchical taxonomy” 또는 “node catalog hierarchy”로 바꾸는 편이 과학적으로 정확하다. |
| 25 | Domain Reviewer | minor | methods | `300+ DPI 고해상도 플로터 20종` | DPI만으로 publication-grade가 되지 않는다. 최종 물리 크기, 벡터 출력, 폰트 임베딩, 선 두께, 색각 접근성, 저널별 파일 형식과 재현 가능한 스타일 설정을 명시하라. |
| 26 | Domain Reviewer | minor | abstract | `8대 핵심 생물정보학 도메인(Bulk RNA-Seq, 단일세포/공간전사체, 생식/체세포 변이 분석, 후성유전학/ATAC-Seq, 메타유전체학, De Novo 유전체 조립, 단백체/대사체학, 구조생물학/CADD)` | 초록의 8개 범주에는 단일세포/공간전사체가 독립 범주지만 Table 1에서는 비교유전체학이 8번째 범주다. 도메인 분류를 한 체계로 통일하라. |
| 27 | Domain Reviewer | minor | references | `{doi:10.1038/nature11247}` | 본문에 인용되지만 참고문헌 목록에 없는 DOI가 4개(`10.1038/nrg.2016.49`, `10.1093/bib/bbw020`, `10.1038/nature11247`, `10.1186/s13059-017-1299-7`)다. 반대로 여러 도구 논문은 목록에만 있고 본문 연결이 없다. 인용-참고문헌 양방향 검사를 수행하라. |
| 28 | Domain Reviewer | major | availability | `GitHub 저장소 (MIT 라이선스)` | 실제 URL, 릴리스/commit, 설치 가능한 패키지, 라이선스 파일, 환경 lockfile, 예제 workflow 및 benchmark 산출물의 영구 식별자가 없다. 소프트웨어 논문 제출 전 공개·고정된 재현 패키지가 필요하다. |

요약: **28건 — major 24, minor 4, suggestion 0**.

## 실제 실행 전 제3장 권장 구성

### 3.1 Evaluation design and reproducibility protocol

- 분석 전에 데이터 accession, 포함/제외 기준, primary endpoint, 허용 오차와 실패 기준을 고정한다.
- 하드웨어, OS, ComfyBioflow commit, 각 도구 버전, Conda lock/container digest, reference genome 및 DB build를 표로 제시한다.
- 모든 실행은 원시 입력 → 노드 명령 → stdout/stderr → 산출물 → 지표 계산으로 추적되게 한다.
- reference implementation을 같은 입력·도구 버전·파라미터로 실행해 플랫폼 wrapper가 결과를 바꾸지 않는지 비교한다.

### 3.2 Analytical concordance across representative workflows

| Domain | Correct input | Primary endpoint | 필수 보조 지표 |
|---|---|---|---|
| RNA-seq | 정확한 SEQC A/B biological replicates 및 TaqMan/ERCC truth | 사전 고정 truth set에 대한 정량 상관 또는 DEG concordance | slope/bias, ROC/PR, sensitivity/specificity와 95% CI, mapping/quantification QC |
| PhiX assembly | 실제 PhiX paired reads + NC_001422.1 | QUAST genome fraction와 base-level error | misassemblies, mismatches/100 kb, indels/100 kb, total length, circularization |
| ATAC-seq | 정확한 GM12878 ATAC biological replicates | reference/replicate peak concordance | FRiP, TSS enrichment, duplicate/mitochondrial fraction, IDR, read-depth sensitivity |
| Zymo metagenome | 제조사 조성과 연결된 shotgun run | species-level abundance distance | precision/recall, false positives, L1/Bray–Curtis, DB-version sensitivity, subsampling CI |

각 subsection은 **Dataset → Workflow and parameters → Prespecified endpoint → Result with uncertainty → Reference-pipeline concordance → Failure analysis** 순서로 통일한다.

### 3.3 Platform-level validation

- **Socket safety:** 허용·거부 연결 행렬과 invalid-edge 테스트 수/통과율을 보고한다.
- **Cache correctness:** 파라미터 변경 유형별 재실행 노드, fresh run 대비 출력 hash/수치 동등성, warm/cold wall time을 최소 10회 반복해 중앙값·IQR·p95로 제시한다.
- **Reproducibility:** 깨끗한 두 환경 또는 두 호스트에서 동일 workflow를 재실행하고 deterministic output은 hash, 수치 output은 사전 허용 오차로 비교한다.
- **Robustness:** 누락 파일, 손상 FASTQ, 잘못된 metadata, tool failure 및 중단 후 재시작을 시험한다.

### 3.4 Usability evaluation 또는 주장 축소

비전공자 대상성을 핵심 기여로 유지하려면 사전 정의 과제 기반 사용자 평가가 필요하다. 최소 결과는 task completion, time-on-task, 오류/도움 요청, SUS 및 재현 성공률이다. 이 연구를 하지 않는다면 “user-friendly/편의성 입증” 대신 “visual interface를 제공한다”라는 기능 기술로 제한한다.

### 3.5 결과 요약

한 표에서 모든 primary endpoint, estimate, uncertainty, reference baseline 및 pass/fail을 보여준다. “완벽”, “100% 정확도”, “골드 스탠다드” 같은 총괄 수사는 사용하지 않는다.

## 고찰 권장 구성

1. **Principal findings:** 실제 primary endpoint가 지지하는 범위만 요약한다.
2. **Comparison with existing platforms:** Table 2의 측정 결과와 공식 문서로 확인된 기능만 비교하고, 로컬 반응형 캔버스가 제공하는 차별점을 좁게 정의한다.
3. **Trade-offs:** 로컬 Conda의 간편성 대 컨테이너/HPC 확장성, 시각 상호작용 대 대규모 분산 처리, 엄격한 타입 대 도구 호환성의 균형을 논한다.
4. **Limitations:** 제한된 대표 데이터셋, 데이터베이스 버전 의존성, OS별 환경 차이, 사용자 연구 규모, 아직 없는 클라우드/HPC 기능을 명시한다.
5. **Future work:** 실제 결과로 확인되지 않은 기능은 결론의 성취가 아니라 향후 계획으로 이동한다.

## 실행 전 즉시 조치

1. 제3장과 초록의 완료형 결과·수치·“입증” 표현을 모두 `planned evaluation`으로 전환한다.
2. Figure 3을 제거하거나 큰 글씨로 `SIMULATED / STUDY DESIGN`이라고 표시한다.
3. 잘못된 네 accession을 교정한 뒤 다운로드 파일의 organism, assay, layout, BioProject를 자동 검증한다.
4. `run_real_benchmarks.py`가 실제 CLI를 실행하고 실패 시 결과를 쓰지 않도록 교체한다.
5. 노드 수, 소켓 타입, sidecar 및 Conda 환경 주장을 현재 코드에서 자동 생성되는 표로 제한한다.

## 외부 사실 확인 출처

- NCBI GEO: [GSE47792 SEQC SuperSeries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47792)
- NCBI GEO: [GSE47726 rheumatoid arthritis microarray study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47726)
- NCBI SRA: [SRR11038933](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11038933)
- NCBI SRA: [SRR12359623](https://www.ncbi.nlm.nih.gov/sra/?term=SRR12359623)
- ENCODE: [ENCSR000EMS experiment record](https://www.encodeproject.org/experiments/ENCSR000EMS/)
- Galaxy: [local and distributed job configuration](https://docs.galaxyproject.org/en/release_25.1/admin/jobs.html)
- Galaxy: [InteractiveTools](https://docs.galaxyproject.org/en/release_26.0/admin/special_topics/interactivetools.html)

## 등록 상태

이번 세션에는 `mcp__co_scientist__add_review`와 `mcp__co_scientist__list_reviews`가 노출되지 않았고 `~/.co-scientist/config.toml`도 없어 Firestore에 위 28건을 등록하거나 기존 AI 리뷰와 중복 검사를 수행할 수 없었다. 원고가 연결된 co-scientist 프로젝트에서 MCP를 활성화하면 위 표의 각 행을 별도 `source="ai"` 리뷰로 등록하고 등록 수를 재검증해야 한다.
