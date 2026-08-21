# 최신 프론티어 LLM 기반 워크플로우 생성 스킬 및 플랫폼 연계 분석

## 1. 개요: LLM 워크플로우 생성 패러다임의 진화 (2024–2026)

초기 거대 언어 모델(LLM)을 활용한 생물정보학 자동화는 단순한 "프롬프트 투 스크립트(Prompt-to-Script)" 방식에 의존하였다. 그러나 이 방식은 도구 환각(Hallucination), 비호환 라이브러리 조합, 중간 입출력 계약 붕괴, 불투명한 런타임 오류로 인해 실제 연구 현장에서의 채택에 큰 한계를 드러냈다.

최신 프론티어 LLM(Anthropic Claude, OpenAI GPT, Google Gemini)은 단순 텍스트 생성을 넘어 **구조화된 에이전트 스킬(Scientific Agent Skills)**과 **엄격한 스키마 디코딩(Constrained Decoding / Tool Calling)**을 지원하도록 진화하였다.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                 Frontier LLM Agentic Workflow Skills Layer                  │
│                                                                             │
│  [Skill 1: Brief Extraction]       [Skill 2: Tool Selection (TSR)]          │
│   • Multi-omics modality parser     • Evidence-tier ranking (REF/ALT)       │
│   • Assay/platform identifier       • Domain rule constraints               │
│                                                                             │
│  [Skill 3: Scope Guard & Repair]   [Skill 4: DAG Synthesis & Typing]       │
│   • Unsupported assay rejection     • Socket contract enforcement           │
│   • Incompatible tool auto-fix      • Topological ordering & acyclic check  │
│                                                                             │
│  [Skill 5: Parameter Optimization] [Skill 6: Interactive Plot Refinement]   │
│   • Library/sequencer arguments     • Real-time downstream feedback loop    │
└──────────────────────────────────────┬──────────────────────────────────────┘
                                       │ Standardized Workflow Spec (JSON)
                                       ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                 ComfyBioflow Platform & Execution Substrate                 │
│                                                                             │
│  • Reactive Visual Node Canvas (ComfyUI-based DAG engine)                   │
│  • Socket Type Contract Verifier (FASTQ, BAM, VCF, Count Matrix, BED, Image)│
│  • Domain-Isolated Conda Environments (env_rnaseq, env_variant, etc.)       │
│  • Smart Topological Upstream Caching (ms-level partial re-execution)       │
│  • 20+ Tier-1 Publication-Grade Visualizers (300+ DPI IMAGE tensors)        │
│  • Provenance Sidecar & Image-Embedded Workflow Replication                 │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. 3대 프론티어 LLM의 워크플로우 생성 스킬 특성 비교

| 모델 패밀리 | 핵심 스킬 메커니즘 | 워크플로우 생성 강점 | 생물정보학 플랫폼 연계 전략 |
| :--- | :--- | :--- | :--- |
| **Anthropic Claude**<br>(Claude 3.5 Sonnet / Opus) | • 구조화된 `SKILL.md` 에이전트 스킬<br>• 고급 추론 및 Chain-of-Thought<br>• 다중 서브에이전트 조율 | • 복잡한 생물학적 제약조건 추론<br>• 도구 간 미묘한 호환성 판별<br>• 맥락 충실도 및 안전성 우수 | • 복합 다단계 파이프라인 기획자(Planner)<br>• 도구 순위 및 증거 등급 결정기 |
| **OpenAI GPT**<br>(GPT-4o, GPT-5/o1/o3) | • Strict JSON Schema Structured Outputs<br>• Python Code Interpreter 샌드박스<br>• 정밀 Function Calling | • 문법/스키마 오류 0% 보장<br>• 신속한 메타데이터 파싱<br>• 정형 JSON 워크플로우 직렬화 | • 엄격한 스키마 기반 워크플로우 컴파일러<br>• 파라미터 매핑 엔진 |
| **Google Gemini**<br>(Gemini 2.5 / 3.7 Flash & Pro) | • 초장문 컨텍스트 (1M+ 토큰)<br>• 네이티브 다중 모달(Multimodal) 그라운딩<br>• 실시간 코드 실행 샌드박스 | • 방대한 도구 문서/참조 시퀀스 검색<br>• 생성된 플롯 도판의 시각적 검증<br>• 전체 실행 로그 분석 및 디버깅 | • 시각적 결과물(Plot/QC) 분석 평가자<br>• 다기관 벤치마크 및 대용량 로그 감사자 |

---

## 3. ComfyBioflow의 6대 핵심 워크플로우 생성 스킬 (Core Agent Skills)

ComfyBioflow 플랫폼은 최신 LLM들이 공통으로 활용할 수 있는 6대 선언적 워크플로우 생성 스킬셋을 표준 인터페이스로 제공합니다:

### Skill 1: 자연어 구조화 및 분석 컨텍스트 추출 (`brief_extraction`)
- **역할**: 비정형 자연어 요청("SEQC 6개 샘플로 대조군-처리군 DEG 분석하고 Volcano 플롯 그려줘")을 표준 데이터 모델(`AnalysisBrief`)로 추출.
- **추출 엔티티**: 오믹스 모달리티(`bulk_rna_seq`), 리드 타입(`paired-end`), 참조 유전체(`GRCh38`), 실험 디자인(`condition: control vs treated`), 핵심 유전자 목록, 필터링 임계값.

### Skill 2: 도구 선택 레지스트리(TSR) 기반 후보 도구 선별 (`tool_selection`)
- **역할**: 오픈소스 생물정보학 커뮤니티(nf-core, GATK, ENCODE) 표준에 부합하는 도구 조합을 선별.
- **증거 등급(Evidence Tier)**:
  - **REF (Reference)**: 표준 권장 도구 (예: Bulk RNA-Seq 정량의 `Salmon`, 메타게놈 분류의 `Kraken2`).
  - **ALT (Alternative)**: 특정 조건에서 정당화되는 대안 도구 (예: 스플라이스 정렬의 `STAR`, 다배체 변이의 `FreeBayes`).
  - **INVALID**: 데이터 타입과 양립할 수 없는 도구 (예: 대형 진핵생물에 대한 Isolate SPAdes 적용).

### Skill 3: 지원범위 가드 및 결함 자동 복구 (`scope_guard_and_repair`)
- **Scope Guard**: 플랫폼이 지원하지 않거나 생물학적으로 모순되는 요청(예: 미생물 데이터에 단일세포 공간 전사체 알고리즘 요청)을 조기 감지하고 거부.
- **Safe Repair**: 사용자의 비표준/결함 요청(예: 단일 말단 리드에 paired 전용 트리머 지정, 폐기된 레거시 도구 요청)을 표준 실행 가능 도구로 자동 치환하고 변경 사유(Repair Diff)를 기록.

### Skill 4: 비순환 방향성 그래프(DAG) 합성 및 포트 타입 검증 (`dag_synthesis_and_typing`)
- **역할**: 선별된 도구들을 ComfyUI 커스텀 노드로 매핑하고 노드 간 입출력 슬롯을 유효한 간선(Edge)으로 연결.
- **소켓 타입 계약 검증**: `FASTQ_PAIR` $\rightarrow$ `SAM_BAM_INDEXED` $\rightarrow$ `VCF_FILTERED` $\rightarrow$ `COUNT_MATRIX_CSV` $\rightarrow$ `IMAGE_TENSOR` 계약을 확인하여 비정상 연결을 차단.

### Skill 5: 실행 파라미터 최적화 및 메타데이터 바인딩 (`parameter_binding`)
- **역할**: 시퀀싱 플랫폼 오류율, 리드 길이, 스레드 수, 어댑터 서열 등 도메인별 최적 인자를 노드 위젯 값으로 자동 바인딩.

### Skill 6: 출판급 시각화 도판 실시간 미세 조정 (`interactive_refinement`)
- **역할**: 생성된 300+ DPI 플롯(IMAGE 텐서)을 분석하고, 사용자의 후속 피드백("p-value 컷오프를 0.01로 높이고 유전자 5개 추가 라벨링해줘")을 받아 상류 재실행 없이 하류 플로터 노드만 즉각 갱신.
