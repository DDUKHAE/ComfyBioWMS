# 생물정보학 노드 리팩토링 계획

기준일: 2026-09-04

## 결정 사항

- 파일 경계는 `Python 라이브러리 1개` 또는 `CLI 실행 파일 1개`이다.
- 각 구현 파일은 형제 노드, `utils.py`, `base.py`, `execution.py`, `bioflow`를 import하지 않는다.
- 각 파일이 입력 검증, 실행, 로그, 오류 처리와 두 ComfyUI 매핑을 직접 가진다.
- FASTA/FASTQ/BAM/VCF/GFF 등은 메모리 객체가 아닌 경로 `STRING`으로 전달한다.
- 호환성을 보존하지 않는다. 기존 워크플로우는 새 인터페이스를 기준으로 다시 만든다.
- 공식 데이터 E2E가 통과한 노드만 중앙 registry에 등록한다. 현재 등록 범위는 14개다.

## 기존 파일 → 도구/라이브러리별 파일

굵은 파일은 1차 구현과 E2E가 완료됐다.

| 기존 묶음 | 변경 파일 | 단계 |
|---|---|---|
| `biopython_nodes.py`, `biopython_sequence_info.py` | **`class_1/biopython.py`**, `class_2/blast.py`, `class_1/logomaker.py`, `biotite.py`, `pyhmmer.py`, `dna_features_viewer.py`, `primer3.py`, `pyfaidx.py`, `pycircos.py`, `edlib.py`, `codonw.py`, `pytfbs.py`, `kegg.py` | Biopython 완료, 나머지 미등록 |
| `ref_nodes.py` | 입력검증, **`class_2/fastp.py`**, **`fastqc.py`**, `trimmomatic.py`, `salmon.py`, `tximport.py`, `deseq2.py`, `class_1/scanpy.py`, `cellranger.py`, 시각화, 보고서 | fastp/FastQC 완료 |
| `variant_nodes.py` | 입력검증, **`class_2/bwa_mem2.py`**, **`samtools.py`**, **`bcftools.py`**, 시각화, 보고서 | CLI 3종 완료 |
| `assembly_nodes.py` | 입력검증, **`fastp.py`**, **`class_2/spades.py`**, **`class_2/quast.py`**, 시각화, 보고서 | fastp, SPAdes, QUAST 완료 |
| `atac_nodes.py` | 입력검증, **`fastp.py`**, **`bwa_mem2.py`**, **`samtools.py`**, `macs3.py`, 시각화, 보고서 | 공통 CLI 완료 |
| `metagenome_nodes.py` | 입력검증, **`fastp.py`**, `kraken2.py`, `bracken.py`, 시각화, 보고서 | fastp 완료 |
| `genomics_longread_nodes.py` | `pysam.py`, `cyvcf2.py`, `pybedtools.py`, `mappy.py`, `pyfastx.py`, `seqkit.py`, `bowtie2.py`, `mosdepth.py`, `sniffles2.py`, `cutesv.py`, `flye.py`, `hifiasm.py`, `racon.py`, `medaka.py`, `deepvariant.py` | 미등록 |
| `transcriptomics_spatial_nodes.py` | `anndata.py`, `scvi_tools.py`, `scvelo.py`, `cellrank.py`, `squidpy.py`, `tangram.py`, `cell2location.py`, `pyscenic.py`, `cellphonedb.py`, `gseapy.py`, `muon.py`, `kallisto.py`, `alevin_fry.py`, `star.py`, `stringtie.py`, `flair.py`, `isotools.py`, `cite_seq_count.py`, `edger.py`, `limma.py` | 미등록 |
| `epigenomics_nodes.py` | `deeptools.py`, `tobias.py`, `genrich.py`, `crispresso2.py`, `mageck.py`, `scikit_fusion.py`, `seacr.py`, `homer.py`, `meme.py`, `methyldackel.py`, `cooler.py`, `cooltools.py`, `chromosight.py` | 미등록 |
| `microbiome_nodes.py` | `scikit_bio.py`, `biom_format.py`, `ete3.py`, `fastunifrac.py`, `dada2.py`, `humann.py`, `metaphlan.py`, `genomad.py`, `virsorter2.py`, `checkv.py`, `prokka.py`, `bakta.py`, `amrfinderplus.py`, `augur.py` | 미등록 |
| `proteomics_metabolomics_nodes.py` | `pyteomics.py`, `pyopenms.py`, `matchms.py`, `spec2vec.py`, `massql.py`, `ms_deisotope.py`, `diann.py`, `msfragger.py`, `msconvert.py`, `msdial.py`, `sirius.py`, `maxquant.py`, `perseus.py`, `metaboanalyst.py`, `mhcquant.py` | 미등록 |
| `cadd_structural_nodes.py` | `colabfold.py`, `esmfold.py`, `diffdock.py`, `gnina.py`, `rdkit.py`, `openmm.py`, `mdanalysis.py`, `mdtraj.py`, `torchdrug.py`, `openfold.py`, `prody.py`, `autodock_vina.py`, `p2rank.py`, `fpocket.py`, `plumed.py`, `pmx.py`, `openbabel.py`, `smina.py` | 미등록 |
| `publication_visualizer_nodes.py` | 실제 렌더러 라이브러리별 파일; 모든 노드는 실제 입력 경로 사용 | 미등록 |
| `class_1/*_node.py`, `class_2/*_node.py` | 대응 통합 파일의 공식 데이터 E2E 통과 후 도구별로 제거 | registry에서 제외 |

## Galaxy 기준 파라미터와 `extra_command`

1. Galaxy IUC wrapper의 필수 입력과 자주 쓰는 옵션을 ComfyUI `required`/`optional`에 노출한다.
2. 원본 CLI의 나머지 옵션은 모든 CLI 노드의 optional multiline `extra_command`로 전달한다.
3. 직접 노출한 옵션은 `extra_command`에서 제거한다. 긴 옵션, `--flag=value`, 짧은 옵션과 `-t8` 같은 결합형을 모두 처리하고 제거된 토큰을 stderr에 알린다.
4. `shlex.split` 후 인자 배열로 실행하며 `shell=True`를 쓰지 않는다.
5. 서브커맨드마다 같은 짧은 플래그의 의미가 다르므로 충돌 테이블도 서브커맨드별로 둔다.

| 파일/명령 | 직접 관리하여 `extra_command`에서 제거하는 옵션 |
|---|---|
| fastp | `-i/-I/-o/-O/-w/-q/-u/-n/-l/-j/-h`, 각 long alias, adapter 감지, correction |
| FastQC | `-o/-t/-c/-a/-l/--nogroup/--min_length/-k`와 long alias |
| BWA-MEM2 `mem` | `-t/-x/-k/-w/-T/-R` |
| samtools `sort` | `-@/-m/-n/-N/-o/-O` |
| samtools `index` | `-@/-b/-c/-o` |
| samtools `markdup` | `-@/-r/--mode/-d/-f` |
| bcftools `mpileup` | reference, depth, base/map quality, threads, output 형식/경로 옵션 |
| bcftools `call` | caller, variants-only, ploidy, threads, output 형식/경로 옵션 |
| bcftools `filter` | include/exclude, soft filter, SNP/indel gap, threads, output 형식/경로 옵션 |
| spades | `-1/-2/-s/--12/-o/--output-dir/-t/--threads/-m/--memory/--careful/--sc/--meta/--isolate/--only-assembler/--cov-cutoff/-k/--kmers/--phred-offset` |
| quast | `-o/--output-dir/-r/--reference/-g/--features/-m/--min-contig/-t/--threads/--large` |

## 구현 우선순위와 반복 절차

| 우선순위 | 범위 | 완료 조건 |
|---:|---|---|
| 1 | Biopython, fastp, FastQC, BWA-MEM2, samtools, bcftools | 완료: 6개 파일, 12개 노드 등록 및 공식 데이터 E2E |
| 2 | SPAdes/QUAST, MACS3, Kraken2/Bracken, Salmon/DESeq2 | SPAdes/QUAST 완료 (2개 파일, 2개 노드 추가 등록); 대표 파이프라인 후속 진행 |
| 3 | pysam/cyvcf2/pybedtools/scanpy/anndata/scikit-bio | 라이브러리별 공식 fixture 파싱과 실제 출력 검증 |
| 4 | long-read, single-cell/spatial, epigenomics | 공식 도구 example 또는 nf-core 데이터로 기능별 E2E |
| 5 | proteomics/metabolomics, 구조생물학, 시각화 | 큰 모델/DB 요구량을 명시하고 실행 가능한 CI 계층 분리 |

각 도구는 다음 순서로만 승격한다.

1. 공식 wrapper 버전과 upstream/nf-core fixture commit 및 SHA-256 고정
2. 독립 import, UI 기본값, 충돌 제거, 바이너리 부재, 실제 E2E 실패 테스트 작성
3. 한 도구 파일에 최소 구현
4. 노드 메서드를 호출해 산출물을 해당 공식 parser/CLI로 재검증
5. 검증 매트릭스 갱신 후 registry 등록
