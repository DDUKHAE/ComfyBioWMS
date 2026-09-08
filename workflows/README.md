# Case-study workflows

## 최종 권장 축소 세트

논문과 공개 결과에는 아래 세 워크플로만 후보로 사용합니다. JSON의 `extra.benchmark`에 주장 수준과
제약을 함께 기록했으며, 현재 입력이 부적합한 PhiX 사례는 실행 그래프를 만들지 않았습니다.

```bash
python engine/scripts/generate_case_study_workflows.py final
```

| 파일 | 용도 | 현재 상태 |
|---|---|---|
| `final_01_pancreas_scanpy_e2e.json` | 실제 Pancreas E15.5 H5AD 기반 Scanpy 기술 E2E | 실행 후보. scVelo 원 결과 재현은 아님 |
| `final_02_zymo_metaphlan_e2e.json` | SRR12324253 100,000 read pairs 기반 분류 E2E | 전체 MetaPhlAn DB를 별도로 준비해야 함 |
| `final_03_gse110004_chri_smoke.json` | GSE110004 4-run, chr-I 축소 RNA-seq 스모크 | 실행 후보. 원 논문 재현은 아님 |
| `final_case_study_status.json` | 채택·조건부·제외 상태 매니페스트 | PhiX `blocked` 사유 포함 |

`data/paper_phix174/reads_R*.fastq.gz`(SRR5458066)는 사람 피부 16S amplicon 자료이므로 PhiX174
조립 입력으로 사용하면 안 됩니다. `data/nf_core_taxprofiler`의 초소형 분류 DB도 Zymo D6300 전체
조성 정확도 평가에는 사용하지 않습니다.

## 레거시 5-케이스 생성 세트

nf-core 파이프라인의 공식 테스트 데이터, 그리고 각 원 논문의 accession에 앵커된 5종.
워크플로 JSON은 `engine/scripts/generate_case_study_workflows.py`가 레지스트리를 introspect해서
생성합니다 — 노드 시그니처를 바꿨으면 JSON을 손으로 고치지 말고 스크립트를 다시 돌리세요
(생성 시 구조 검증 포함).

```bash
python data/download_case_study_data.py                      # nf-core 공식 테스트 데이터 (기본)
python data/download_case_study_data.py --source paper       # 문헌/공개 accession 실험 세트
python engine/scripts/generate_case_study_workflows.py       # cs*.json 재생성
python engine/scripts/generate_case_study_workflows.py paper # cs*_paper.json 재생성
```

두 소스는 파일명이 겹치므로 **디렉터리를 분리**합니다 (`data/nf_core_*` vs `data/paper_*`).
`cs*_paper.json`은 입력이 nf-core 테스트 데이터와 달라지는 CS-2·3·4의 레거시 실험 그래프입니다 —
CS-1은 논문 데이터(airway)에 genome/GTF가 없어 STAR·HISAT2 경로가 성립하지 않고, CS-5는 두
소스가 동일 파일입니다.

## 데이터 출처

`data/download_case_study_data.py`는 두 소스를 지원하며, 모든 URL은 GitHub tree API로 실제 경로를
확인하고 range 요청으로 응답을 검증했습니다. 받은 파일의 sha256/크기는 `data/manifest.lock.json`에
기록되어 재실행 시 무결성 검사에 쓰입니다.

**`--source nfcore` (기본)** — `nf-core/test-datasets`, 파이프라인 브랜치별 커밋 SHA 고정
(2026-09-07 기준). nf-core CI가 실제로 돌리는 그 데이터이고 가볍습니다.

| 케이스 | 브랜치 @ SHA | 내용 |
|---|---|---|
| CS-1 | `rnaseq` @ `e07c1b1` | GSE110004 yeast RAP1, 4 샘플 paired (WT×2 / RAP1_IAA_30M×2) + genome/GTF/transcriptome |
| CS-2 | `taxprofiler` @ `c9601c9` | Kraken2 / Bracken / sylph 테스트 DB + paired FASTQ |
| CS-3 | `sarek` @ `24cdbea` | `human_g1k_v37_decoy.small.fasta` + tiny normal FASTQ |
| CS-4 | `bacass` @ `32853ba` | ERR044595 (1M reads, paired) |
| CS-5 | `theislab/scvelo_notebooks` | pancreas E15.5 h5ad (테스트 데이터 = 논문 데이터) |

**`--source paper`** — 원 논문 accession. 용량이 커서 FASTQ는 `--max-reads`(기본 100,000)로
스트리밍 서브샘플합니다.

| 케이스 | 앵커 | Accession |
|---|---|---|
| CS-1 | Himes et al. 2014 airway ASM (DESeq2 표준 데이터셋) | GSE52778 / SRR1039508·09·12·13 + GENCODE v19 |
| CS-2 | ZymoBIOMICS D6300 공개 SRA run(Nicholls et al. 2019 자료와 동일하지 않음) | SRR12324253 |
| CS-3 | GIAB HG001 / NA12878, NIST v4.2.1 benchmark | ERR194147 + truth VCF/BED + GRCh38 chr20 |
| CS-4 | Sanger 1977 PhiX174 완전서열 | NC_001422.1만 유효; SRR5458066은 PhiX 입력으로 부적합 |
| CS-5 | Bastidas-Ponce et al. 2019 | 동일 h5ad |

```bash
python data/download_case_study_data.py --source paper --max-reads 200000
python data/download_case_study_data.py cs2 cs4        # 일부만
```

## 케이스별 합격 기준

| 파일 | 앵커 | 합격 기준 |
|---|---|---|
| `cs1_nfcore_rnaseq_route_concordance.json` | nf-core/rnaseq; Patro 2017 / Dobin 2013 / Love 2014 | 3경로 카운트 Spearman ρ ≥ 0.95 |
| `cs2_nfcore_taxprofiler_multiprofiler.json` | nf-core/taxprofiler; Wood 2019 | 실 DB 준비 후 기대 10종 검출 및 조성 오차 평가 |
| `cs3_nfcore_sarek_variant_calling.json` | nf-core/sarek; GIAB | (`--source paper`) high-conf 영역 recall/precision, Ti/Tv 2.0–2.1 |
| `cs4_nfcore_bacass_assembly.json` | nf-core/bacass; Bankevich 2012 / Gurevich 2013 | 올바른 PhiX 리드 준비 후 최대 contig 5,386 bp ±1%, misassembly 0 |
| `cs5_scanpy_pancreas_endocrinogenesis.json` | Luecken & Theis 2019 / Wolf 2018 | 저자 라벨 대비 ARI ≥ 0.7 |

정량 기준은 각 사례에 맞는 정답 자료가 별도로 확인된 경우에만 평가합니다. `--source paper`라는
이름 자체가 ground truth를 보장하지 않습니다. nf-core 테스트 데이터는 실행 완주 및 경로 간
일치도(CS-1) 검증용입니다.

## 실행 전 확인

1. **`PreviewAny` 종단 노드**: 결과 경로와 요약을 캔버스에서 확인할 수 있도록 모든 분기를 ComfyUI
   코어 `PreviewAny`(IMAGE는 `PreviewImage`)로 끝냈습니다. 설치된 ComfyUI에
   `PreviewAny`가 없으면 버전을 올리거나 종단 노드를 교체하세요.
2. **Conda 환경**: Class 2 노드는 `engine/envs/*.yaml`의 격리 환경이 필요합니다.
3. **경로**: 위젯 경로는 리포지터리 루트 기준 상대경로입니다.

## 워크플로별 미비점

- **CS-1**: 디렉터리(`data/nf_core_rnaseq`) 기반 배치 처리 기능이 지원되어, 단일 `TrimGalore` 노드에서 4개 샘플(WT×2 / RAP1_IAA_30M×2) 페어를 자동 인식해 트리밍하고, 단일 `SalmonQuantReads` 노드가 일괄 정량하여 `Tximport` → `DESeq2SampleQC` / `DESeq2` → `VolcanoPlot`으로 바로 연결됩니다(이전처럼 샘플별 중복 노드 및 `JoinPaths` 노드 배치가 불필요). STAR·HISAT2 경로 역시 다중 샘플 배치를 자연스럽게 수용합니다. RSEM 경로는 `rsem-prepare-reference` 노드가 없어 제외.
- **CS-2**: host removal 단계는 nf-core 테스트 데이터에 host 레퍼런스가 없어 제외했습니다.
  분류 DB는 두 소스가 공유합니다(`data/nf_core_taxprofiler/`). **주의**: 이 테스트 DB는 극소형이라
  Zymo mock의 기대 구성원 전체를 담고 있지 않습니다. 전체 구성원 검출 기준을 실제로 평가하려면 Kraken2
  Standard(8GB+) 같은 실 DB가 필요합니다.
  `MicrobiomeStackedBar`는 Bracken TSV를 그대로 받으므로 `taxonomy_id` 같은 수치 컬럼까지
  샘플로 그립니다 — 발표용 도표에는 abundance 행렬 재구성이 필요합니다.
- **CS-3**: GIAB truth 교차(`PybedtoolsIntersect`)는 `cs3_..._paper.json`에만 있습니다.
  레퍼런스는 truth VCF와 좌표계를 맞추기 위해 GRCh38 **chr20만** 받습니다(65MB) — ERR194147
  리드 대부분은 chr20에 매핑되지 않으므로 recall은 chr20 한정으로 해석해야 합니다. truth VCF와의
  변이 단위 비교는 `bcftools isec` 노드가 없어 BED 교차까지만 그래프 안에서 처리됩니다.
  `Cyvcf2Stats`는 VCF 변이의 전이/전좌(transition/transversion) 빈도를 직접 집계하여 `ti_tv_ratio`를 계산합니다.
- **CS-4**: 기본(nfcore) 데이터에는 매칭 레퍼런스가 없어 QUAST를 reference-free로 돌립니다.
  레거시 `cs4_..._paper.json`은 NC_001422.1과 관련 없는 SRR5458066을 연결하므로 사용하지 마십시오.
  올바른 PhiX paired read set이 확보되기 전까지 최종 세트에서는 `blocked`입니다.
- **CS-5**: `ScanpyCluster`가 UMAP PNG를 직접 내므로 `UmapScatter`는 연결하지 않았습니다
  (embedding CSV를 내보내는 노드가 없음).

## 아직 막힌 케이스

**ATAC-seq (nf-core/atacseq)**: `nodes/class_2/macs3.py`의 `Macs3Callpeak`이 완성되어 있고
`engine/envs/epigenomics.yaml`에 `macs3=3.0.4`도 있으나 `nodes/class_2/__init__.py`에 import이
빠져 레지스트리에 없습니다. 배선하면 CS-6이 열립니다.
