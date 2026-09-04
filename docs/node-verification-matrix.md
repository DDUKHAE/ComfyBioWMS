# 노드별 공식 데이터 E2E 검증 매트릭스

기준일: 2026-09-04. `tests/official_data.json`의 URL과 SHA-256을 실행 전에 검증한다.

## 필수 실행 환경

| 그룹 | 요구조건 |
|---|---|
| Python | Python 3.11+, Biopython 1.88, pytest |
| QC | fastp 1.3.6; FastQC 0.12.1; Java runtime (`JAVA_HOME` 필요 가능) |
| Variant | Bioconda `bwa-mem2=2.3`, samtools 1.24, bcftools 1.24가 `PATH`에 존재 |
| Galaxy 기준 | `galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608` |

주의: 현재 Bioconda 패키지 메타데이터는 BWA-MEM2 2.3이지만 번들 실행 파일의 `bwa-mem2 version` 출력은 2.2.1이다. 기능 E2E는 통과했으나 이 불일치를 해소하기 전까지 런타임 버전 일치로 표기하지 않는다.

## 검증 결과

| 노드 | 대상 | 공식 테스트 데이터 URL | 입력 파라미터 | 기대 출력물 검증 기준 | 상태 |
|---|---|---|---|---|---|
| `BiopythonSeqIOStatsNode` | Biopython 1.88 | [Biopython `dups.fasta`](https://raw.githubusercontent.com/biopython/biopython/dc262b5c437e07a8cc1cfb8a734c0d84a4434b23/Tests/Fasta/dups.fasta) | FASTA, `fasta` | JSON 파싱, 레코드 5개, 고유 ID 4개 | 통과 |
| `BiopythonAlignmentStatsNode` | Biopython 1.88 | [Biopython `protein.aln`](https://raw.githubusercontent.com/biopython/biopython/dc262b5c437e07a8cc1cfb8a734c0d84a4434b23/Tests/Clustalw/protein.aln) | alignment, `clustal` | 20행, 411열, identity 0–100 | 통과 |
| `FastpNode` | fastp 1.3.6 | [fastp R1](https://raw.githubusercontent.com/OpenGene/fastp/dce5c40b861db8c649081e72a7b6449358a83775/testdata/R1.fq), [R2](https://raw.githubusercontent.com/OpenGene/fastp/dce5c40b861db8c649081e72a7b6449358a83775/testdata/R2.fq) | paired, threads 2, q15, u40, n5, l15 | FASTQ 2개 파싱, JSON의 입력 reads > 0, HTML 비어 있지 않음 | 통과 |
| `FastQCNode` | FastQC 0.12.1 + Java | [FastQC `minimal.fastq`](https://raw.githubusercontent.com/s-andrews/FastQC/87fb3364a2f37115833d678648926d41e184f0b1/test/data/minimal.fastq) | threads 2, kmers 7 | HTML/ZIP 비어 있지 않음, ZIP 내 `fastqc_data.txt`, `##FastQC` 헤더 | 통과 |
| `BwaMem2IndexNode` | BWA-MEM2, 패키지 2.3/버전 출력 2.2.1 | [nf-core Sarek reference](https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/reference/human_g1k_v37_decoy.small.fasta) | reference, output directory | 복사 FASTA와 5종 index sidecar 비어 있지 않음 | 기능 통과, 버전 주석 |
| `BwaMem2AlignNode` | BWA-MEM2, 패키지 2.3/버전 출력 2.2.1 | 위 reference, [Sarek R1](https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/testdata/dummy/normal/dummy_n_R1_xxx.fastq.gz), [R2](https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/testdata/dummy/normal/dummy_n_R2_xxx.fastq.gz) | paired, threads 2, Illumina, k19, w100, T30 | SAM 헤더와 alignment record 존재 | 기능 통과, 버전 주석 |
| `SamtoolsSortNode` | samtools 1.24 | 위 BWA-MEM2 SAM | coordinate, threads 2, 768M | BAM 생성 후 `samtools quickcheck` | 통과 |
| `SamtoolsIndexNode` | samtools 1.24 | 위 sorted BAM | BAI, threads 2 | index 비어 있지 않음, `samtools idxstats` 성공 | 통과 |
| `SamtoolsMarkdupNode` | samtools 1.24 | 위 sorted BAM | keep duplicates, template(`t`), distance 100 | collate→fixmate→sort→markdup, 최종 `quickcheck` | 통과 |
| `BcftoolsMpileupNode` | bcftools 1.24 | 위 reference + indexed BAM | d250, Q13, q0, threads 1 | BCF 비어 있지 않음, `bcftools view -h` 성공 | 통과 |
| `BcftoolsCallNode` | bcftools 1.24 | 위 mpileup BCF | multiallelic, variants only, default ploidy | VCF 생성, `bcftools view -h` 성공 | 통과 |
| `BcftoolsFilterNode` | bcftools 1.24 | 위 called VCF | exclude `QUAL<10` | VCF 헤더 검증, 출력 record 수 ≤ 입력 | 통과 |

## 실행 명령

```bash
pytest -q -m 'not e2e' tests/test_official_data.py tests/test_standalone_*.py tests/test_classified_nodes.py tests/test_node_mappings.py

JAVA_HOME=/opt/miniconda3/envs/bulk_rna_seq/lib/jvm \
PATH=/opt/miniconda3/envs/bulk_rna_seq/bin:$PATH \
pytest -q -m e2e tests/test_standalone_biopython.py tests/test_standalone_fastp.py tests/test_standalone_fastqc.py

PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH \
pytest -q -m e2e tests/test_standalone_bwa_mem2.py tests/test_standalone_samtools.py tests/test_standalone_bcftools.py
```

등록되지 않은 기존 `*_node.py` 파일은 이 표의 통과 범위가 아니며, 공식 데이터 E2E 전까지 사용 가능 또는 구현 완료로 간주하지 않는다.
해당 레거시 노드와 예전 워크플로우를 전제로 한 테스트도 보존만 하며 기본 pytest 수집에서 제외한다.
