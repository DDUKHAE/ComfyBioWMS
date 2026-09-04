# ComfyBIOWMS

ComfyUI에서 실제 생물정보학 라이브러리와 CLI를 실행하는 커스텀 노드 모음입니다.

현재 중앙 registry에는 공식 upstream 또는 nf-core 데이터로 End-to-End 검증한 14개 노드만 등록됩니다. 이전 워크플로우와 클래스 이름의 호환성은 보장하지 않으며, 워크플로우는 새 인터페이스를 기준으로 다시 작성합니다.

## 검증된 노드

| 파일 | 노드 |
|---|---|
| `nodes/class_1/biopython.py` | sequence file statistics, alignment statistics |
| `nodes/class_2/fastp.py` | FASTQ filtering/trimming |
| `nodes/class_2/fastqc.py` | FASTQ/BAM/SAM quality control |
| `nodes/class_2/bwa_mem2.py` | reference index, read alignment |
| `nodes/class_2/samtools.py` | sort, index, markdup |
| `nodes/class_2/bcftools.py` | mpileup, call, filter |
| `nodes/class_2/spades.py` | de novo genome assembly |
| `nodes/class_2/quast.py` | genome assembly quality assessment |

각 파일은 단독 복사 가능한 구조이며 자체 입력 검증, 실행기, 오류 처리, `NODE_CLASS_MAPPINGS`, `NODE_DISPLAY_NAME_MAPPINGS`를 포함합니다. 대용량 생물정보학 데이터는 경로 `STRING`으로 전달합니다.

CLI 노드는 Galaxy IUC wrapper의 자주 쓰는 파라미터를 UI에 노출하고, 나머지 CLI 옵션을 위한 optional `extra_command`를 제공합니다. UI가 직접 관리하는 옵션을 `extra_command`에 다시 입력하면 해당 토큰은 제거되고 stderr에 안내됩니다.

## 설치 요구조건

- Python 3.11+ 및 Biopython 1.88
- fastp 1.3.6
- FastQC 0.12.1 및 Java runtime
- BWA-MEM2 2.3, samtools 1.24, bcftools 1.24
- SPAdes 4.3.0, QUAST 5.3.0

실행할 CLI는 ComfyUI 프로세스의 `PATH`에 있어야 합니다. 각 노드 파일 상단 docstring에 Python 및 외부 바이너리 요구조건이 기록되어 있습니다.

## 문서

- [리팩토링 및 파일 재구성 계획](docs/node-refactoring-plan.md)
- [공식 데이터 E2E 검증 매트릭스](docs/node-verification-matrix.md)
- [상세 1차 구현 계획](docs/superpowers/plans/2026-09-04-standalone-bioinformatics-nodes-phase-1.md)

## 테스트

구조 및 비-E2E 검사:

```bash
pytest -q -m 'not e2e' tests/test_official_data.py tests/test_standalone_*.py tests/test_classified_nodes.py tests/test_node_mappings.py
```

실제 CLI E2E 명령은 [검증 매트릭스](docs/node-verification-matrix.md)에 환경별로 기록되어 있습니다. 테스트 데이터는 `tests/official_data.json`의 immutable URL과 SHA-256으로 검증하며, 임의 생성 또는 fallback biological fixture를 사용하지 않습니다.

## 상태 정책

`nodes/class_1/*_node.py`와 `nodes/class_2/*_node.py`에 남은 파일은 마이그레이션 후보 소스입니다. 공식 데이터 E2E를 통과해 registry에 추가되기 전에는 검증된 노드가 아닙니다.

동일하게 기존 workflow/legacy-node 테스트 파일은 보존만 하며 기본 pytest 수집 대상이 아닙니다. 새 인터페이스로 워크플로우를 다시 만들 때 공식 데이터 E2E로 교체합니다.

## License

[MIT](LICENSE)
