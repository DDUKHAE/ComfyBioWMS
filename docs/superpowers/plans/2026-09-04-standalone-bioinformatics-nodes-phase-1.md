# Standalone Bioinformatics Nodes Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the first high-use library and CLI node cluster with standalone, Galaxy-aligned, genuinely E2E-verified ComfyUI modules.

**Architecture:** Keep `class_1/` and `class_2/`, but consolidate nodes by actual library or executable. Every production module contains its own validation, execution code, dependency documentation, and ComfyUI mappings; the root registry exposes only modules whose real E2E checks passed.

**Tech Stack:** Python 3.11+, ComfyUI node API, Biopython 1.88, fastp 1.3.6, FastQC 0.12.1, BWA-MEM2 2.3, samtools 1.24, bcftools 1.24, pytest, Galaxy IUC wrapper XML, upstream and nf-core official datasets.

**Spec:** `docs/superpowers/specs/2026-09-04-standalone-bioinformatics-nodes-design.md`

## Global Constraints

- Existing class names, input names, and workflow JSON compatibility do not need to be preserved.
- Preserve all unrelated uncommitted user changes; every commit uses `git commit --only` with explicit paths.
- Production node modules may import only the Python standard library and their named third-party library; no imports from `nodes`, `bioflow`, or sibling files.
- One Python library or one CLI executable owns one implementation file.
- Every implementation file defines `NODE_CLASS_MAPPINGS` and `NODE_DISPLAY_NAME_MAPPINGS` at its bottom.
- CLI executables are resolved from the active `PATH` with `shutil.which()` and executed without `shell=True`.
- Biological result streams go directly to files; human-readable logs are captured and displayed live.
- Every CLI node places `extra_command` in `optional` and removes flags managed by direct node inputs.
- No random, synthetic, handwritten, dry-run, or fallback biological output is permitted.
- Only real E2E-verified nodes enter the production `NODE_CLASS_MAPPINGS`.

---

## Repository-Wide Reorganization Map

This table is the master migration map. Phase 1 implements the bold targets; later phases follow the same contract.

| Existing bundled source | Tool/library files after migration |
|---|---|
| `nodes/assembly_nodes.py` | `input_validation.py`, **`fastp.py`**, `spades.py`, `quast.py`, `matplotlib_visualizers.py`, `reporting.py` |
| `nodes/atac_nodes.py` | `input_validation.py`, **`fastp.py`**, **`bwa_mem2.py`**, **`samtools.py`**, `macs3.py`, `matplotlib_visualizers.py`, `reporting.py` |
| `nodes/biopython_nodes.py`, `nodes/biopython_sequence_info.py` | **`biopython.py`**, `blast.py`, `logomaker.py`, `biotite.py`, `pyhmmer.py`, `dna_features_viewer.py`, `primer3.py`, `pyfaidx.py`, `pycircos.py`, `edlib.py`, `codonw.py`, `pytfbs.py`, `kegg.py`, `matplotlib_visualizers.py` |
| `nodes/cadd_structural_nodes.py` | `colabfold.py`, `esmfold.py`, `diffdock.py`, `gnina.py`, `rdkit.py`, `openmm.py`, `mdanalysis.py`, `mdtraj.py`, `torchdrug.py`, `openfold.py`, `prody.py`, `autodock_vina.py`, `p2rank.py`, `fpocket.py`, `plumed.py`, `pmx.py`, `openbabel.py`, `smina.py` |
| `nodes/epigenomics_nodes.py` | `deeptools.py`, `tobias.py`, `genrich.py`, `crispresso2.py`, `mageck.py`, `scikit_fusion.py`, `seacr.py`, `homer.py`, `meme.py`, `methyldackel.py`, `cooler.py`, `cooltools.py`, `chromosight.py` |
| `nodes/genomics_longread_nodes.py` | `pysam.py`, `cyvcf2.py`, `pybedtools.py`, `mappy.py`, `pyfastx.py`, `seqkit.py`, `bowtie2.py`, `mosdepth.py`, `sniffles2.py`, `cutesv.py`, `flye.py`, `hifiasm.py`, `racon.py`, `medaka.py`, `deepvariant.py` |
| `nodes/metagenome_nodes.py` | `input_validation.py`, **`fastp.py`**, `kraken2.py`, `bracken.py`, `matplotlib_visualizers.py`, `reporting.py` |
| `nodes/microbiome_nodes.py` | `scikit_bio.py`, `biom_format.py`, `ete3.py`, `fastunifrac.py`, `dada2.py`, `humann.py`, `metaphlan.py`, `genomad.py`, `virsorter2.py`, `checkv.py`, `prokka.py`, `bakta.py`, `amrfinderplus.py`, `augur.py` |
| `nodes/proteomics_metabolomics_nodes.py` | `pyteomics.py`, `pyopenms.py`, `matchms.py`, `spec2vec.py`, `massql.py`, `ms_deisotope.py`, `diann.py`, `msfragger.py`, `msconvert.py`, `msdial.py`, `sirius.py`, `maxquant.py`, `perseus.py`, `metaboanalyst.py`, `mhcquant.py` |
| `nodes/publication_visualizer_nodes.py` | `matplotlib_visualizers.py`, `seaborn_visualizers.py`, `biopython.py` for phylogenetic parsing; every renderer consumes an actual input path |
| `nodes/ref_nodes.py` | `input_validation.py`, **`fastp.py`**, **`fastqc.py`**, `trimmomatic.py`, `salmon.py`, `tximport.py`, `deseq2.py`, `scanpy.py`, `cellranger.py`, `matplotlib_visualizers.py`, `reporting.py` |
| `nodes/transcriptomics_spatial_nodes.py` | `anndata.py`, `scvi_tools.py`, `scvelo.py`, `cellrank.py`, `squidpy.py`, `tangram.py`, `cell2location.py`, `pyscenic.py`, `cellphonedb.py`, `gseapy.py`, `muon.py`, `kallisto.py`, `alevin_fry.py`, `star.py`, `stringtie.py`, `flair.py`, `isotools.py`, `cite_seq_count.py`, `edger.py`, `limma.py` |
| `nodes/variant_nodes.py` | `input_validation.py`, **`bwa_mem2.py`**, **`samtools.py`**, **`bcftools.py`**, `matplotlib_visualizers.py`, `reporting.py` |
| Current `nodes/class_1/*_node.py`, `nodes/class_2/*_node.py` | Removed tool-by-tool only after the consolidated replacement passes and becomes the registered implementation |

## Phase 1 Verification Matrix

All URLs are immutable. SHA-256 values are recorded in `tests/official_data.json` and checked before use.

| Node | Tool/library | Galaxy parameter source | Official test data | Inputs | Output validation |
|---|---|---|---|---|---|
| `BiopythonSeqIOStatsNode` | Biopython 1.88 | Not applicable: library node | `biopython/Tests/Fasta/dups.fasta` at `dc262b5c` | FASTA path, `fasta` | Parse returned JSON; exactly 5 records and 4 distinct IDs |
| `BiopythonAlignmentStatsNode` | Biopython 1.88 | Not applicable: library node | `biopython/Tests/Clustalw/protein.aln` at `dc262b5c` | alignment path, `clustal` | Alignment has 20 rows and 411 columns; identity is in `[0,100]` |
| `FastpNode` | fastp 1.3.6 | `tools-iuc/tools/fastp/{fastp.xml,macros.xml}` at `6a1769b` | `OpenGene/fastp/testdata/R1.fq,R2.fq` at `dce5c40b` | paired FASTQ, quality 15, unqualified 40%, N limit 5, min length 15 | Both output FASTQs parse; JSON parses and reports nonzero before-filtering reads; HTML nonempty |
| `FastQCNode` | FastQC 0.12.1 | `tools-iuc/tools/fastqc/rgFastQC.xml` at `6a1769b` | `s-andrews/FastQC/test/data/minimal.fastq` at `87fb3364` | FASTQ, threads 2, kmers 7 | HTML and ZIP nonempty; ZIP contains `fastqc_data.txt`; report begins `##FastQC` |
| `BwaMem2IndexNode` | BWA-MEM2 2.3 | `tools-iuc/tools/bwa_mem2/bwa-mem2-idx.xml` at `6a1769b` | `nf-core/test-datasets/reference/human_g1k_v37_decoy.small.fasta` at `24cdbea4` | reference path, output directory | Copied reference and all BWA-MEM2 index sidecars are nonempty |
| `BwaMem2AlignNode` | BWA-MEM2 2.3 | `tools-iuc/tools/bwa_mem2/bwa-mem2.xml` at `6a1769b` | same reference plus Sarek dummy paired FASTQ at `24cdbea4` | indexed reference, R1/R2, threads 2, Illumina preset | SAM is nonempty; `samtools view -c` succeeds; header and at least one alignment exist |
| `SamtoolsSortNode` | samtools 1.24 | `tools-iuc/tool_collections/samtools/samtools_sort/samtools_sort.xml` at `6a1769b` | BWA-MEM2 SAM above | SAM, coordinate order, threads 2 | `samtools quickcheck` passes; header sort order is `coordinate` |
| `SamtoolsIndexNode` | samtools 1.24 | native Galaxy collection at `6a1769b`; index has no dedicated current wrapper | sorted BAM | BAM, threads 2, BAI format | `.bai` nonempty; `samtools idxstats` succeeds |
| `SamtoolsMarkdupNode` | samtools 1.24 | `tools-iuc/tool_collections/samtools/samtools_markdup/samtools_markdup.xml` at `6a1769b` | sorted BAM above | BAM, threads 2, remove duplicates false, mode template | BAM passes `quickcheck`; `flagstat` parses; duplicate count is nonnegative |
| `BcftoolsMpileupNode` | bcftools 1.24 | `tools-iuc/tools/bcftools/bcftools_mpileup.xml` at `6a1769b` | indexed reference and marked BAM above | reference, BAM, max depth 250, base quality 13, map quality 0 | BCF nonempty; `bcftools view -h` succeeds |
| `BcftoolsCallNode` | bcftools 1.24 | `tools-iuc/tools/bcftools/bcftools_call.xml` at `6a1769b` | mpileup BCF above | multiallelic caller, variants only true, default diploid ploidy | VCF nonempty; header validates; each record has REF and ALT |
| `BcftoolsFilterNode` | bcftools 1.24 | `tools-iuc/tools/bcftools/bcftools_filter.xml` at `6a1769b` | called VCF above | exclude `QUAL<10`, output VCF | `bcftools view -h` succeeds; output record count does not exceed input |

Exact Phase 1 source checksums:

```json
{
  "biopython_dups.fasta": "fad2a6461314182c94d176d96f591574fd705c1e635595d67d3279a67f2c6cee",
  "biopython_protein.aln": "61cf11d97a68848248f78677228319bbede85cbd7f5a0c1f299c61369e9b6a19",
  "fastp_R1.fq": "d0a93989b8f7350efe11197f4a4783da87d1fa1280a7d7b33a953801ea28b0e3",
  "fastp_R2.fq": "eaa72e8178c8f2a915b706cc887e4cfd46920c04af6c644cc69d45de455ef7d4",
  "fastqc_minimal.fastq": "b26ccc879f1f9276d7dec09da7024a87a3881e5ad4b943480722e99cf10aed45",
  "sarek_ref.fasta": "e6da32e7a11739023a513ec4f4bf32a96d13a0cb7cbb2b64c2367bca3ca550ff",
  "sarek_R1.fastq.gz": "c44e3211b1b10eb14ce13a268b5512d7f15bfd4adbe291a0a44f8294444b263a",
  "sarek_R2.fastq.gz": "c31608e380f9b44d9c4315289144707b6b09dbdbec4b24145315f877b47292f4"
}
```

---

### Task 1: Official Dataset Fetcher and Provenance Manifest

**Files:**

- Create: `tests/official_data.json`
- Create: `tests/official_data.py`
- Create: `tests/test_official_data.py`
- Modify: `pyproject.toml`

**Interfaces:**

- Produces: `fetch_official_data(name: str, cache_dir: Path) -> Path`
- Consumes: JSON entries containing `url` and `sha256`.

- [ ] **Step 1: Write the failing checksum test**

```python
from tests.official_data import fetch_official_data

def test_fetches_biopython_fixture_with_verified_checksum(tmp_path):
    path = fetch_official_data("biopython_dups.fasta", tmp_path)
    assert path.read_text().startswith(">")
    assert path.stat().st_size == 129
```

- [ ] **Step 2: Run the test and confirm the missing module failure**

Run: `pytest tests/test_official_data.py -v`

Expected: FAIL because `tests.official_data` does not exist.

- [ ] **Step 3: Add the immutable manifest and minimum standard-library fetcher**

```python
import hashlib
import json
import urllib.request
from pathlib import Path

_MANIFEST = json.loads(Path(__file__).with_name("official_data.json").read_text())

def fetch_official_data(name: str, cache_dir: Path) -> Path:
    item = _MANIFEST[name]
    target = cache_dir / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        with urllib.request.urlopen(item["url"], timeout=60) as response, target.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    if actual != item["sha256"]:
        target.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {name}: {actual}")
    return target
```

Write this exact `official_data.json` content:

```json
{
  "biopython_dups.fasta": {
    "url": "https://raw.githubusercontent.com/biopython/biopython/dc262b5c437e07a8cc1cfb8a734c0d84a4434b23/Tests/Fasta/dups.fasta",
    "sha256": "fad2a6461314182c94d176d96f591574fd705c1e635595d67d3279a67f2c6cee"
  },
  "biopython_protein.aln": {
    "url": "https://raw.githubusercontent.com/biopython/biopython/dc262b5c437e07a8cc1cfb8a734c0d84a4434b23/Tests/Clustalw/protein.aln",
    "sha256": "61cf11d97a68848248f78677228319bbede85cbd7f5a0c1f299c61369e9b6a19"
  },
  "fastp_R1.fq": {
    "url": "https://raw.githubusercontent.com/OpenGene/fastp/dce5c40b861db8c649081e72a7b6449358a83775/testdata/R1.fq",
    "sha256": "d0a93989b8f7350efe11197f4a4783da87d1fa1280a7d7b33a953801ea28b0e3"
  },
  "fastp_R2.fq": {
    "url": "https://raw.githubusercontent.com/OpenGene/fastp/dce5c40b861db8c649081e72a7b6449358a83775/testdata/R2.fq",
    "sha256": "eaa72e8178c8f2a915b706cc887e4cfd46920c04af6c644cc69d45de455ef7d4"
  },
  "fastqc_minimal.fastq": {
    "url": "https://raw.githubusercontent.com/s-andrews/FastQC/87fb3364a2f37115833d678648926d41e184f0b1/test/data/minimal.fastq",
    "sha256": "b26ccc879f1f9276d7dec09da7024a87a3881e5ad4b943480722e99cf10aed45"
  },
  "sarek_ref.fasta": {
    "url": "https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/reference/human_g1k_v37_decoy.small.fasta",
    "sha256": "e6da32e7a11739023a513ec4f4bf32a96d13a0cb7cbb2b64c2367bca3ca550ff"
  },
  "sarek_R1.fastq.gz": {
    "url": "https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/testdata/dummy/normal/dummy_n_R1_xxx.fastq.gz",
    "sha256": "c44e3211b1b10eb14ce13a268b5512d7f15bfd4adbe291a0a44f8294444b263a"
  },
  "sarek_R2.fastq.gz": {
    "url": "https://raw.githubusercontent.com/nf-core/test-datasets/24cdbea48c4415a29f668a724ce602fef8fed813/testdata/dummy/normal/dummy_n_R2_xxx.fastq.gz",
    "sha256": "c31608e380f9b44d9c4315289144707b6b09dbdbec4b24145315f877b47292f4"
  }
}
```

Register `e2e` in `pyproject.toml` under `tool.pytest.ini_options.markers`.

- [ ] **Step 4: Run the fetcher test**

Run: `pytest tests/test_official_data.py -v`

Expected: PASS with one downloaded, checksum-verified Biopython fixture.

- [ ] **Step 5: Commit only task files**

```bash
git add tests/official_data.json tests/official_data.py tests/test_official_data.py pyproject.toml
git commit --only tests/official_data.json tests/official_data.py tests/test_official_data.py pyproject.toml -m "test: add immutable official dataset fixtures"
```

### Task 2: Standalone Biopython Module

**Files:**

- Create: `nodes/class_1/biopython.py`
- Create: `tests/test_standalone_biopython.py`

**Interfaces:**

- Produces: `BiopythonSeqIOStatsNode.run(sequence_file: str, file_format: str) -> tuple[str, int]`
- Produces: `BiopythonAlignmentStatsNode.run(alignment_file: str, file_format: str) -> tuple[str, int, int, float]`
- Produces: module-local `NODE_CLASS_MAPPINGS` and `NODE_DISPLAY_NAME_MAPPINGS`.

- [ ] **Step 1: Write failing standalone and real-data tests**

```python
import importlib.util
import json
from pathlib import Path
from tests.official_data import fetch_official_data

MODULE = Path("nodes/class_1/biopython.py").resolve()

def load_module():
    spec = importlib.util.spec_from_file_location("standalone_biopython", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def test_biopython_file_loads_standalone_and_owns_mappings():
    module = load_module()
    assert set(module.NODE_CLASS_MAPPINGS) == {
        "BiopythonSeqIOStatsNode", "BiopythonAlignmentStatsNode"
    }
    assert set(module.NODE_DISPLAY_NAME_MAPPINGS) == set(module.NODE_CLASS_MAPPINGS)

def test_seqio_stats_parse_official_fasta(tmp_path):
    module = load_module()
    fasta = fetch_official_data("biopython_dups.fasta", tmp_path)
    summary, count = module.BiopythonSeqIOStatsNode().run(str(fasta), "fasta")
    rows = json.loads(summary)
    assert count == 5
    assert len({row["id"] for row in rows}) == 4

def test_alignment_stats_parse_official_clustal(tmp_path):
    module = load_module()
    alignment = fetch_official_data("biopython_protein.aln", tmp_path)
    summary, rows, columns, identity = module.BiopythonAlignmentStatsNode().run(str(alignment), "clustal")
    assert rows == 20
    assert columns == 411
    assert 0 <= identity <= 100
    assert json.loads(summary)["rows"] == 20
```

- [ ] **Step 2: Confirm failure because the consolidated file is missing**

Run: `pytest tests/test_standalone_biopython.py -v`

Expected: FAIL opening `nodes/class_1/biopython.py`.

- [ ] **Step 3: Implement only the two verified streaming parsers**

Create `nodes/class_1/biopython.py` with this behavior:

```python
"""Biopython nodes.

Python packages: biopython==1.88
External binaries: none
"""
import json
from itertools import combinations
from pathlib import Path

def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path

class BiopythonSeqIOStatsNode:
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("summary_json", "sequence_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {
            "sequence_file": ("STRING", {"default": ""}),
            "file_format": (["fasta", "fastq", "genbank", "embl"], {"default": "fasta"}),
        }}

    def run(self, sequence_file: str, file_format: str = "fasta"):
        from Bio import SeqIO
        rows = []
        for record in SeqIO.parse(str(_file(sequence_file, "Sequence input")), file_format):
            sequence = str(record.seq).upper()
            gc = sequence.count("G") + sequence.count("C")
            rows.append({"id": record.id, "description": record.description,
                         "length": len(sequence),
                         "gc_percent": round(100 * gc / len(sequence), 6) if sequence else 0.0})
        return json.dumps(rows), len(rows)

class BiopythonAlignmentStatsNode:
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT", "FLOAT")
    RETURN_NAMES = ("summary_json", "rows", "columns", "mean_pairwise_identity")

    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {
            "alignment_file": ("STRING", {"default": ""}),
            "file_format": (["clustal", "fasta", "phylip", "stockholm"], {"default": "clustal"}),
        }}

    def run(self, alignment_file: str, file_format: str = "clustal"):
        from Bio import AlignIO
        alignment = AlignIO.read(str(_file(alignment_file, "Alignment input")), file_format)
        identities = []
        for left, right in combinations(alignment, 2):
            pairs = [(a, b) for a, b in zip(left.seq, right.seq) if a != "-" and b != "-"]
            if pairs:
                identities.append(sum(a == b for a, b in pairs) / len(pairs))
        value = round(100 * sum(identities) / len(identities), 6) if identities else 100.0
        summary = {"rows": len(alignment), "columns": alignment.get_alignment_length(),
                   "mean_pairwise_identity": value}
        return json.dumps(summary), summary["rows"], summary["columns"], value

NODE_CLASS_MAPPINGS = {
    "BiopythonSeqIOStatsNode": BiopythonSeqIOStatsNode,
    "BiopythonAlignmentStatsNode": BiopythonAlignmentStatsNode,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "BiopythonSeqIOStatsNode": "Biopython: Sequence File Statistics",
    "BiopythonAlignmentStatsNode": "Biopython: Alignment Statistics",
}
```

- [ ] **Step 4: Run standalone and official-data checks**

Run: `pytest tests/test_standalone_biopython.py -v`

Expected: 3 PASS.

- [ ] **Step 5: Commit only the Biopython slice**

```bash
git add nodes/class_1/biopython.py tests/test_standalone_biopython.py
git commit --only nodes/class_1/biopython.py tests/test_standalone_biopython.py -m "feat: add standalone verified Biopython nodes"
```

### Task 3: Standalone fastp Module with Galaxy Options

**Files:**

- Create: `nodes/class_2/fastp.py`
- Create: `tests/test_standalone_fastp.py`

**Interfaces:**

- Produces: `FastpNode.run(read1, output_dir, read2="", threads=4, qualified_quality_phred=15, unqualified_percent_limit=40, n_base_limit=5, length_required=15, detect_adapter_for_pe=False, correction=False, extra_command="") -> tuple[str, str, str, str]`.
- Managed flags: `-i/--in1`, `-I/--in2`, `-o/--out1`, `-O/--out2`, `-w/--thread`, `-q/--qualified_quality_phred`, `-u/--unqualified_percent_limit`, `-n/--n_base_limit`, `-l/--length_required`, `--detect_adapter_for_pe`, `-c/--correction`, `-j/--json`, `-h/--html`.

- [ ] **Step 1: Write failing interface and duplicate-filter tests**

```python
def test_fastp_galaxy_options_and_extra_command_are_optional(fastp_module):
    inputs = fastp_module.FastpNode.INPUT_TYPES()
    assert set(inputs["required"]) == {"read1", "output_dir"}
    assert inputs["optional"]["qualified_quality_phred"][1]["default"] == 15
    assert inputs["optional"]["unqualified_percent_limit"][1]["default"] == 40
    assert inputs["optional"]["length_required"][1]["default"] == 15
    assert "extra_command" in inputs["optional"]

def test_fastp_removes_managed_options_but_keeps_native_extra(fastp_module):
    kept, ignored = fastp_module._filter_extra(
        "--thread=99 -q 3 --overrepresentation_analysis", fastp_module._MANAGED_OPTIONS
    )
    assert kept == ["--overrepresentation_analysis"]
    assert ignored == ["--thread=99", "-q"]
```

- [ ] **Step 2: Confirm the missing module failure**

Run: `pytest tests/test_standalone_fastp.py -v -k 'not e2e'`

Expected: FAIL opening `nodes/class_2/fastp.py`.

- [ ] **Step 3: Implement the standalone module**

The module docstring pins fastp 1.3.6 and the two IUC wrapper files. Implement the local filtering and execution boundary as follows, then build the command only from validated node inputs plus `kept` tokens:

```python
def _filter_extra(text, managed):
    tokens, kept, ignored = shlex.split(text), [], []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        match = next((flag for flag in managed if token == flag or token.startswith(flag + "=")
                      or (len(flag) == 2 and managed[flag] == 1 and token.startswith(flag) and token != flag)), None)
        if match:
            ignored.append(token)
            i += 1 + (managed[match] if token == match else 0)
        else:
            kept.append(token)
            i += 1
    return kept, ignored

def _run(argv, cwd):
    executable = shutil.which(argv[0])
    if not executable:
        raise RuntimeError("fastp executable not found on PATH; install conda package fastp=1.3.6")
    argv = [executable, *argv[1:]]
    process = subprocess.Popen(argv, cwd=cwd, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, bufsize=1)
    stdout, stderr = [], []
    def drain(pipe, target, collected):
        for line in iter(pipe.readline, ""):
            collected.append(line)
            print(line, end="", file=target, flush=True)
    threads = [threading.Thread(target=drain, args=(process.stdout, sys.stdout, stdout)),
               threading.Thread(target=drain, args=(process.stderr, sys.stderr, stderr))]
    for thread in threads: thread.start()
    for thread in threads: thread.join()
    code = process.wait()
    if code:
        raise RuntimeError(f"fastp exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")
    return "".join(stdout), "".join(stderr)
```

`FastpNode.INPUT_TYPES()` contains only `read1` and `output_dir` in `required`. Put every remaining signature value in `optional`, including the multiline empty `extra_command`. Construct:

```python
argv = ["fastp", "-i", str(read1_path), "-o", str(out1), "-w", str(threads),
        "-q", str(qualified_quality_phred), "-u", str(unqualified_percent_limit),
        "-n", str(n_base_limit), "-l", str(length_required),
        "-j", str(report_json), "-h", str(report_html)]
if read2_path:
    argv += ["-I", str(read2_path), "-O", str(out2)]
if detect_adapter_for_pe: argv.append("--detect_adapter_for_pe")
if correction: argv.append("--correction")
argv += kept
```

After `_run`, parse the JSON report and require nonempty R1, paired R2 when supplied, JSON, and HTML artifacts. End with mappings containing only `FastpNode`.

- [ ] **Step 4: Add and run the real E2E check**

```python
@pytest.mark.e2e
def test_fastp_runs_on_upstream_paired_data(fastp_module, tmp_path):
    r1 = fetch_official_data("fastp_R1.fq", tmp_path)
    r2 = fetch_official_data("fastp_R2.fq", tmp_path)
    out1, out2, report_json, report_html = fastp_module.FastpNode().run(
        str(r1), str(tmp_path / "out"), str(r2), threads=2,
        extra_command="--overrepresentation_analysis --thread 99"
    )
    assert Path(out1).stat().st_size > 0 and Path(out2).stat().st_size > 0
    assert json.loads(Path(report_json).read_text())["summary"]["before_filtering"]["total_reads"] > 0
    assert Path(report_html).stat().st_size > 0
```

Run: `PATH=/opt/miniconda3/envs/bulk_rna_seq/bin:$PATH pytest tests/test_standalone_fastp.py -v`

Expected: all tests PASS and the log warns that `--thread` was ignored.

- [ ] **Step 5: Commit only the fastp slice**

```bash
git add nodes/class_2/fastp.py tests/test_standalone_fastp.py
git commit --only nodes/class_2/fastp.py tests/test_standalone_fastp.py -m "feat: add Galaxy-aligned standalone fastp node"
```

### Task 4: Standalone FastQC Module

**Files:**

- Create: `nodes/class_2/fastqc.py`
- Create: `tests/test_standalone_fastqc.py`

**Interfaces:**

- Produces: `FastQCNode.run(input_file, output_dir, threads=2, contaminants="", adapters="", limits="", nogroup=False, min_length=0, kmers=7, extra_command="") -> tuple[str, str]`.
- Managed flags: `-o/--outdir`, `-t/--threads`, `-c/--contaminants`, `-a/--adapters`, `-l/--limits`, `--nogroup`, `--min_length`, `-k/--kmers`, `-f/--format`, `--extract`.

- [ ] **Step 1: Write failing Galaxy-interface, collision, missing-binary, and E2E tests**

Assert `kmers` has Galaxy bounds 2–10 and default 7; `extra_command` is optional; `--kmers=3 -t 99 --noextract` keeps only `--noextract`; an empty controlled `PATH` raises an error naming FastQC; and the official minimal FASTQ produces a nonempty HTML and ZIP whose `fastqc_data.txt` starts with `##FastQC`.

- [ ] **Step 2: Verify the missing module failure**

Run: `pytest tests/test_standalone_fastqc.py -v -k 'not e2e'`

Expected: FAIL opening `nodes/class_2/fastqc.py`.

- [ ] **Step 3: Implement and verify the minimum module**

Define the managed table and command construction explicitly:

```python
_MANAGED_OPTIONS = {"-o": 1, "--outdir": 1, "-t": 1, "--threads": 1,
                    "-c": 1, "--contaminants": 1, "-a": 1, "--adapters": 1,
                    "-l": 1, "--limits": 1, "--nogroup": 0, "--min_length": 1,
                    "-k": 1, "--kmers": 1, "-f": 1, "--format": 1, "--extract": 0}
argv = ["fastqc", "--outdir", str(out), "--threads", str(threads),
        "--kmers", str(kmers), "--extract"]
for flag, value in (("--contaminants", contaminants), ("--adapters", adapters),
                    ("--limits", limits)):
    if value:
        argv += [flag, str(_file(value, flag))]
if nogroup: argv.append("--nogroup")
if min_length: argv += ["--min_length", str(min_length)]
argv += kept + [str(input_path)]
```

The file-local runner resolves `fastqc`, starts `subprocess.Popen` with separate text stdout/stderr pipes, drains both with two threads, prints lines live, and raises `RuntimeError` with captured stderr on a non-zero exit. Compute the output stem by removing `.gz`, `.bz2`, `.fastq`, `.fq`, `.bam`, or `.sam`, then require `<stem>_fastqc.html` and `<stem>_fastqc.zip`. End with mappings containing only `FastQCNode`.

Run: `PATH=/opt/miniconda3/envs/bulk_rna_seq/bin:$PATH pytest tests/test_standalone_fastqc.py -v`

Expected: all tests PASS.

- [ ] **Step 4: Commit only the FastQC slice**

```bash
git add nodes/class_2/fastqc.py tests/test_standalone_fastqc.py
git commit --only nodes/class_2/fastqc.py tests/test_standalone_fastqc.py -m "feat: add Galaxy-aligned standalone FastQC node"
```

### Task 5: Standalone BWA-MEM2 Module

**Files:**

- Create: `nodes/class_2/bwa_mem2.py`
- Create: `tests/test_standalone_bwa_mem2.py`

**Interfaces:**

- Produces: `BwaMem2IndexNode.run(reference_fasta, output_dir, extra_command="") -> tuple[str]`.
- Produces: `BwaMem2AlignNode.run(indexed_reference, read1, output_sam, read2="", threads=1, preset="illumina", minimum_seed_length=19, band_width=100, minimum_score=30, read_group="", extra_command="") -> tuple[str]`.
- Managed align flags: `-t`, `-x`, `-k`, `-w`, `-T`, `-R`, `-p`, `-I`.

- [ ] **Step 1: Write failing standalone, Galaxy-default, collision, and real index/alignment tests**

Use the pinned Sarek reference and paired reads. Assert the module owns both mappings; defaults match the IUC wrapper's Illumina/full settings; `-t99 -k 5 -Y` preserves only `-Y`; indexing creates real nonempty sidecars; alignment produces a nonempty SAM with at least one non-header record.

- [ ] **Step 2: Confirm the missing module failure**

Run: `pytest tests/test_standalone_bwa_mem2.py -v -k 'not e2e'`

Expected: FAIL opening `nodes/class_2/bwa_mem2.py`.

- [ ] **Step 3: Implement index and alignment without a shell pipeline**

Copy the input FASTA to the requested output directory before indexing and construct the two commands as follows:

```python
index_argv = ["bwa-mem2", "index", *index_extra, str(copied_reference)]

preset_args = [] if preset == "illumina" else ["-x", preset]
align_argv = ["bwa-mem2", "mem", "-t", str(threads), *preset_args,
              "-k", str(minimum_seed_length), "-w", str(band_width),
              "-T", str(minimum_score)]
if read_group: align_argv += ["-R", read_group]
align_argv += align_extra + [str(indexed_reference), str(read1_path)]
if read2_path: align_argv.append(str(read2_path))
```

The file-local index runner resolves `bwa-mem2`, streams both human-readable pipes with threads, and checks exit status. The alignment runner opens `output_sam` in binary mode for stdout, uses a text stderr pipe drained to `sys.stderr`, waits for completion, and deletes the partial SAM on failure. Require nonempty `.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, and `.pac` sidecars after indexing and a nonempty SAM beginning with `@` after alignment. End with mappings for the two BWA-MEM2 classes.

- [ ] **Step 4: Run against BWA-MEM2 2.3**

First verify: `PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH bwa-mem2 version`

Required output: `2.3`. If the environment still reports `2.2.1`, recreate/update it from `engine/envs/variant_analysis.yaml` before marking the matrix row verified.

Run: `PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH pytest tests/test_standalone_bwa_mem2.py -v`

Expected: all tests PASS with BWA-MEM2 2.3.

- [ ] **Step 5: Commit only the BWA-MEM2 slice**

```bash
git add nodes/class_2/bwa_mem2.py tests/test_standalone_bwa_mem2.py
git commit --only nodes/class_2/bwa_mem2.py tests/test_standalone_bwa_mem2.py -m "feat: add standalone BWA-MEM2 nodes"
```

### Task 6: Standalone samtools Module

**Files:**

- Create: `nodes/class_2/samtools.py`
- Create: `tests/test_standalone_samtools.py`

**Interfaces:**

- Produces: `SamtoolsSortNode.run(input_alignment, output_bam, threads=1, sort_order="coordinate", memory_per_thread="768M", extra_command="") -> tuple[str]`.
- Produces: `SamtoolsIndexNode.run(input_bam, output_index="", threads=1, index_format="bai", extra_command="") -> tuple[str]`.
- Produces: `SamtoolsMarkdupNode.run(input_bam, output_bam, threads=1, remove_duplicates=False, mode="template", optical_distance=100, extra_command="") -> tuple[str]`.

- [ ] **Step 1: Write failing mapping, parameter, filtering, and E2E tests**

Assert all three classes share one standalone module; `extra_command` is optional for each; sort filters `-@`, `-m`, `-n/-N`, `-o/-O`; index filters `-@`, `-b/-c`, `-o`; markdup filters `-@`, `-r`, `--mode`, `-d`, `-f`; and a full sort→index→markdup invocation on the BWA output passes `samtools quickcheck` and `idxstats`.

- [ ] **Step 2: Confirm the missing module failure**

Run: `pytest tests/test_standalone_samtools.py -v -k 'not e2e'`

Expected: FAIL opening `nodes/class_2/samtools.py`.

- [ ] **Step 3: Implement the three operations**

Use one file-local runner for the `samtools` executable and construct these commands:

```python
sort_argv = ["samtools", "sort", "-@", str(threads), "-m", memory_per_thread,
             "-O", "BAM", "-o", str(output_bam)]
if sort_order == "name": sort_argv.append("-n")
sort_argv += sort_extra + [str(input_alignment)]

index_argv = ["samtools", "index", "-@", str(threads)]
if index_format == "csi": index_argv.append("-c")
index_argv += index_extra + [str(input_bam), str(output_index)]

commands = [
    ["samtools", "collate", "-@", str(threads), "-o", str(namesort_bam), str(input_bam)],
    ["samtools", "fixmate", "-@", str(threads), "-m", str(namesort_bam), str(fixmate_bam)],
    ["samtools", "sort", "-@", str(threads), "-O", "BAM", "-o", str(coordsort_bam), str(fixmate_bam)],
    ["samtools", "markdup", "-@", str(threads), "--mode", mode,
     "-d", str(optical_distance), *( ["-r"] if remove_duplicates else [] ),
     *markdup_extra, str(coordsort_bam), str(output_bam)],
]
```

Validate each input before creating output parents. Run `samtools quickcheck` after sort and markdup, and `samtools idxstats` after index. Markdup uses `tempfile.TemporaryDirectory(dir=output_parent)` so intermediates are removed automatically; a command failure deletes its partial declared output and raises with captured stderr. End with mappings for all three classes.

- [ ] **Step 4: Run the real tests**

Run: `PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH pytest tests/test_standalone_samtools.py -v`

Expected: all tests PASS using samtools 1.24.

- [ ] **Step 5: Commit only the samtools slice**

```bash
git add nodes/class_2/samtools.py tests/test_standalone_samtools.py
git commit --only nodes/class_2/samtools.py tests/test_standalone_samtools.py -m "feat: add standalone samtools nodes"
```

### Task 7: Standalone bcftools Module

**Files:**

- Create: `nodes/class_2/bcftools.py`
- Create: `tests/test_standalone_bcftools.py`

**Interfaces:**

- Produces: `BcftoolsMpileupNode.run(reference_fasta, input_bam, output_bcf, max_depth=250, min_base_quality=13, min_mapping_quality=0, threads=1, extra_command="") -> tuple[str]`.
- Produces: `BcftoolsCallNode.run(input_bcf, output_vcf, calling_method="multiallelic", variants_only=True, ploidy="default", threads=1, extra_command="") -> tuple[str]`.
- Produces: `BcftoolsFilterNode.run(input_vcf, output_vcf, exclude="QUAL<10", include="", soft_filter="", snp_gap=0, indel_gap=0, threads=1, extra_command="") -> tuple[str]`.

- [ ] **Step 1: Write failing mapping, Galaxy-input, collision, mutual-exclusion, and E2E tests**

Assert all three classes are mapped; all `extra_command` inputs are optional; managed long, equals, and short forms are removed; setting both `include` and `exclude` raises `ValueError`; mpileup→call→filter executes on the real indexed BAM; `bcftools view -h` succeeds for every artifact; and filtered record count is no greater than called record count.

- [ ] **Step 2: Confirm the missing module failure**

Run: `pytest tests/test_standalone_bcftools.py -v -k 'not e2e'`

Expected: FAIL opening `nodes/class_2/bcftools.py`.

- [ ] **Step 3: Implement the three operations**

Use one file-local live runner for `bcftools` and explicit managed-option tables per subcommand because flag meanings differ. Construct:

```python
mpileup_argv = ["bcftools", "mpileup", "-f", str(reference),
                "-d", str(max_depth), "-Q", str(min_base_quality),
                "-q", str(min_mapping_quality), "--threads", str(threads),
                "-Ob", "-o", str(output_bcf), *mpileup_extra, str(input_bam)]

call_argv = ["bcftools", "call", "-m" if calling_method == "multiallelic" else "-c"]
if variants_only: call_argv.append("-v")
if ploidy != "default": call_argv += ["--ploidy", ploidy]
call_argv += ["--threads", str(threads), "-Ov", "-o", str(output_vcf),
              *call_extra, str(input_bcf)]

filter_argv = ["bcftools", "filter", "--threads", str(threads), "-Ov",
               "-o", str(output_vcf)]
if exclude: filter_argv += ["-e", exclude]
if include: filter_argv += ["-i", include]
if soft_filter: filter_argv += ["-s", soft_filter]
if snp_gap: filter_argv += ["-g", str(snp_gap)]
if indel_gap: filter_argv += ["-G", str(indel_gap)]
filter_argv += filter_extra + [str(input_vcf)]
```

Emit BCF/VCF directly with `-o`, validate each output using a second real `bcftools view -h` call, and remove a partial output on command failure. End with mappings for all three classes.

- [ ] **Step 4: Run the real E2E tests**

Run: `PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH pytest tests/test_standalone_bcftools.py -v`

Expected: all tests PASS using bcftools 1.24.

- [ ] **Step 5: Commit only the bcftools slice**

```bash
git add nodes/class_2/bcftools.py tests/test_standalone_bcftools.py
git commit --only nodes/class_2/bcftools.py tests/test_standalone_bcftools.py -m "feat: add standalone bcftools nodes"
```

### Task 8: Verified Registry, Documentation, and Final Audit

**Files:**

- Modify: `nodes/class_1/__init__.py`
- Modify: `nodes/class_2/__init__.py`
- Modify: `nodes/registry.py`
- Modify: `tests/test_classified_nodes.py`
- Modify: `tests/test_node_mappings.py`
- Create: `docs/node-refactoring-plan.md`
- Create: `docs/node-verification-matrix.md`
- Modify: `README.md`
- Modify: `supplementary_table_s1_node_catalog.md`

**Interfaces:**

- Produces: the root `NODE_CLASS_MAPPINGS` containing exactly the 12 Phase 1 verified nodes.
- Produces: auditable documentation linking every registered node to its passing E2E command and immutable inputs.

- [ ] **Step 1: Replace fixed-count tests with a failing verified-registry contract**

```python
EXPECTED = {
    "BiopythonSeqIOStatsNode", "BiopythonAlignmentStatsNode", "FastpNode",
    "FastQCNode", "BwaMem2IndexNode", "BwaMem2AlignNode", "SamtoolsSortNode",
    "SamtoolsIndexNode", "SamtoolsMarkdupNode", "BcftoolsMpileupNode",
    "BcftoolsCallNode", "BcftoolsFilterNode",
}

def test_registry_contains_only_phase_1_verified_nodes():
    from nodes.registry import NODE_CLASS_MAPPINGS
    assert set(NODE_CLASS_MAPPINGS) == EXPECTED
```

- [ ] **Step 2: Run and confirm the current 179-node registry fails the contract**

Run: `pytest tests/test_classified_nodes.py tests/test_node_mappings.py -v`

Expected: FAIL because legacy and fallback-backed nodes remain registered.

- [ ] **Step 3: Aggregate only module-owned verified mappings**

Replace the package aggregators with explicit verified merges:

```python
# nodes/class_1/__init__.py
from .biopython import NODE_CLASS_MAPPINGS as _biopython_classes
from .biopython import NODE_DISPLAY_NAME_MAPPINGS as _biopython_names
CLASS_1_NODE_MAPPINGS = {**_biopython_classes}
CLASS_1_NODE_DISPLAY_NAME_MAPPINGS = {**_biopython_names}

# nodes/class_2/__init__.py
from .fastp import NODE_CLASS_MAPPINGS as _fastp_classes
from .fastqc import NODE_CLASS_MAPPINGS as _fastqc_classes
from .bwa_mem2 import NODE_CLASS_MAPPINGS as _bwa_classes
from .samtools import NODE_CLASS_MAPPINGS as _samtools_classes
from .bcftools import NODE_CLASS_MAPPINGS as _bcftools_classes
from .fastp import NODE_DISPLAY_NAME_MAPPINGS as _fastp_names
from .fastqc import NODE_DISPLAY_NAME_MAPPINGS as _fastqc_names
from .bwa_mem2 import NODE_DISPLAY_NAME_MAPPINGS as _bwa_names
from .samtools import NODE_DISPLAY_NAME_MAPPINGS as _samtools_names
from .bcftools import NODE_DISPLAY_NAME_MAPPINGS as _bcftools_names
CLASS_2_NODE_MAPPINGS = {**_fastp_classes, **_fastqc_classes, **_bwa_classes,
                         **_samtools_classes, **_bcftools_classes}
CLASS_2_NODE_DISPLAY_NAME_MAPPINGS = {**_fastp_names, **_fastqc_names, **_bwa_names,
                                      **_samtools_names, **_bcftools_names}

# nodes/registry.py
from .class_1 import CLASS_1_NODE_MAPPINGS, CLASS_1_NODE_DISPLAY_NAME_MAPPINGS
from .class_2 import CLASS_2_NODE_MAPPINGS, CLASS_2_NODE_DISPLAY_NAME_MAPPINGS
NODE_CLASS_MAPPINGS = {**CLASS_1_NODE_MAPPINGS, **CLASS_2_NODE_MAPPINGS}
NODE_DISPLAY_NAME_MAPPINGS = {
    **CLASS_1_NODE_DISPLAY_NAME_MAPPINGS, **CLASS_2_NODE_DISPLAY_NAME_MAPPINGS
}
```

Keep legacy source files present but unregistered during the staged migration. Remove fixed node-count assertions and old alias mappings.

- [ ] **Step 4: Write the requested plan and matrix documents**

Copy the repository-wide reorganization map into `docs/node-refactoring-plan.md`. Copy the 12-row Phase 1 matrix into `docs/node-verification-matrix.md`, adding the exact test command, observed tool version, execution date, and `verified` status only for commands that passed. List every remaining legacy class as `unregistered` with its target tool file and migration phase; do not attach a fabricated dataset URL.

Correct README and supplementary-table claims so they distinguish `verified`, `unregistered`, and `blocked`. Remove claims of 179 active nodes, zero placeholders, and universal E2E verification.

- [ ] **Step 5: Run complete verification freshly**

```bash
pytest -q -m 'not e2e' tests/test_official_data.py tests/test_standalone_biopython.py tests/test_standalone_fastp.py tests/test_standalone_fastqc.py tests/test_standalone_bwa_mem2.py tests/test_standalone_samtools.py tests/test_standalone_bcftools.py tests/test_classified_nodes.py tests/test_node_mappings.py
PATH=/opt/miniconda3/envs/bulk_rna_seq/bin:$PATH pytest -q -m e2e tests/test_standalone_biopython.py tests/test_standalone_fastp.py tests/test_standalone_fastqc.py
PATH=/opt/miniconda3/envs/variant_analysis/bin:$PATH pytest -q -m e2e tests/test_standalone_bwa_mem2.py tests/test_standalone_samtools.py tests/test_standalone_bcftools.py
rg -n 'np\.random|random\.|placeholder|mock return' nodes/class_1/biopython.py nodes/class_2/{fastp,fastqc,bwa_mem2,samtools,bcftools}.py
rg -n '^from \.|^from nodes|^from bioflow' nodes/class_1/biopython.py nodes/class_2/{fastp,fastqc,bwa_mem2,samtools,bcftools}.py
```

Expected: pytest commands report zero failures; both `rg` commands return exit code 1 with no matches.

- [ ] **Step 6: Commit only registration and documentation files**

```bash
git add nodes/class_1/__init__.py nodes/class_2/__init__.py nodes/registry.py tests/test_classified_nodes.py tests/test_node_mappings.py docs/node-refactoring-plan.md docs/node-verification-matrix.md README.md supplementary_table_s1_node_catalog.md
git commit --only nodes/class_1/__init__.py nodes/class_2/__init__.py nodes/registry.py tests/test_classified_nodes.py tests/test_node_mappings.py docs/node-refactoring-plan.md docs/node-verification-matrix.md README.md supplementary_table_s1_node_catalog.md -m "refactor: register only standalone verified nodes"
```

## Later Phase Order

1. Assembly: SPAdes and QUAST.
2. Metagenomics: Kraken2 and Bracken with an official miniature database.
3. Single-cell: AnnData, Scanpy, then one package per spatial/trajectory library.
4. Long-read genomics and structural variants.
5. Epigenomics and functional screens.
6. Proteomics and metabolomics.
7. Structural biology and docking; proprietary or GPU-bound tools remain blocked until their real runtime is available.
8. Real-input visualization modules.
9. Delete the final legacy per-node file and shared helper only after its registered replacement has passed.
10. Rebuild workflows from the final verified registry.
