# Supplementary Table S1: Verified ComfyBIOWMS Nodes

기준일: 2026-09-04. 이 표는 공식 upstream 또는 nf-core 데이터로 실제 노드 메서드를 실행하고 산출물을 재검증한 중앙 registry 항목만 포함한다.

| Node class | Tool/library | Functional verification |
|---|---|---|
| `BiopythonSeqIOStatsNode` | Biopython 1.88 | Official FASTA parsed; record/ID counts checked |
| `BiopythonAlignmentStatsNode` | Biopython 1.88 | Official Clustal alignment dimensions and identity checked |
| `FastpNode` | fastp 1.3.6 | Official paired FASTQ processed; FASTQ/JSON/HTML checked |
| `FastQCNode` | FastQC 0.12.1 | Official FASTQ processed; HTML/ZIP/data file checked |
| `BwaMem2IndexNode` | BWA-MEM2 package 2.3 | nf-core reference indexed; sidecars checked |
| `BwaMem2AlignNode` | BWA-MEM2 package 2.3 | nf-core paired FASTQ aligned; SAM parsed |
| `SamtoolsSortNode` | samtools 1.24 | BAM produced; `quickcheck` passed |
| `SamtoolsIndexNode` | samtools 1.24 | BAI produced; `idxstats` passed |
| `SamtoolsMarkdupNode` | samtools 1.24 | Real markdup chain executed; `quickcheck` passed |
| `BcftoolsMpileupNode` | bcftools 1.24 | BCF produced and parsed |
| `BcftoolsCallNode` | bcftools 1.24 | VCF produced and parsed |
| `BcftoolsFilterNode` | bcftools 1.24 | Filtered VCF parsed; record monotonicity checked |

데이터 URL, 체크섬, 입력값, 상세 판정 기준 및 BWA-MEM2 버전 표시 주석은 [`docs/node-verification-matrix.md`](docs/node-verification-matrix.md)에 있다. 그 외 소스 파일은 마이그레이션 후보이며 검증 완료로 주장하지 않는다.
