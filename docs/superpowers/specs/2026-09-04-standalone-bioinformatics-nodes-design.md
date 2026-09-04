# Standalone Bioinformatics Nodes Design

## Goal

Refactor the current bioinformatics custom nodes into independently copyable ComfyUI modules, with one Python library or one CLI tool per file, and expose only nodes that pass real end-to-end verification with an official dataset.

## Current State

- `nodes/class_1/` and `nodes/class_2/` contain 190 node implementation files created by an existing uncommitted refactor.
- 184 implementation files import another project module such as `base.py`, `execution.py`, `helpers/`, or a stage-command module.
- All 192 Python files in those directories lack their own `NODE_CLASS_MAPPINGS` and `NODE_DISPLAY_NAME_MAPPINGS` definitions.
- At least 24 nodes contain random or fallback result generation. Examples include the edgeR, limma-voom, Cell Ranger, and publication-visualizer nodes.
- The shared command runner captures subprocess output after completion but does not stream it to the terminal. Most CLI node files do not independently check executables with `shutil.which()`.
- The current structural tests pass, but they verify ComfyUI attributes and central registration rather than standalone loading or real tool execution.

The existing uncommitted work is user-owned and must be preserved while files are migrated incrementally.

## Decisions

### Compatibility

Existing node class names, input field names, and workflow JSON compatibility do not need to be retained. Workflows will be rebuilt after the node implementation is complete.

### Migration Strategy

Use a verified-registration strategy:

1. Consolidate related implementations into one file per actual library or CLI tool.
2. Make that file self-contained.
3. Run its standalone and E2E checks.
4. Add it to the package-level registry only after those checks pass.
5. Remove the superseded per-node files for that tool only after registration moves successfully.

Unverified implementations may remain temporarily in the working tree, but they must not be described as verified or newly added to the production registry.

### Target Layout

Retain the user's `class_1/` and `class_2/` grouping while changing the unit of each implementation file:

```text
nodes/
  class_1/
    biopython.py
    scanpy.py
    pysam.py
    ...
  class_2/
    fastp.py
    fastqc.py
    bwa_mem2.py
    samtools.py
    bcftools.py
    ...
```

A library or tool file may contain multiple ComfyUI node classes when the library/tool exposes multiple operations. Pipeline-specific aliases such as separate assembly, ATAC, and metagenome fastp files will not be retained.

Files for validators or reports that use only the Python standard library are grouped by their actual responsibility, for example `input_validation.py` or `reporting.py`. A visualization file is assigned to the library that performs the rendering and must consume real input data; it may not generate demonstration data.

## Standalone File Contract

Every migrated implementation file must satisfy all of the following without importing another local module:

- A module docstring lists required Python packages and required external executables, including suggested Conda and Apt package names where applicable.
- All node classes, private validation helpers, and execution helpers required by the file are defined in that file.
- `NODE_CLASS_MAPPINGS` and `NODE_DISPLAY_NAME_MAPPINGS` are defined at the bottom of the file.
- Loading the file directly with `importlib.util.spec_from_file_location()` succeeds.
- Copying the file outside this repository does not change its behavior.
- Large biological inputs and generated artifacts cross node boundaries as path `STRING` values rather than in-memory file contents.

This explicit independence requires small duplication between CLI files. The duplication is limited to path validation and subprocess execution; no shared utility layer is introduced.

## CLI Execution Contract

Each CLI tool file defines the minimum local runner needed by its nodes:

1. Resolve the executable using `shutil.which()` and raise a tool-specific `RuntimeError` when absent.
2. Validate every input path with `Path.exists()` and the expected file/directory kind before creating output directories.
3. Create the requested output directory with `Path.mkdir(parents=True, exist_ok=True)`.
4. Build an argument list and tokenize optional arguments with `shlex.split()`; never use `shell=True`.
5. Execute with `subprocess.Popen()`.
6. Drain human-readable stdout and stderr concurrently, preserve their text, and print each stream to the corresponding terminal stream while the process is running. When stdout is the biological result stream, write it directly to the declared output file instead of loading it into memory or flooding the terminal; stderr remains live and captured.
7. Raise an exception containing the executable, exit code, command, and captured stderr when the exit status is non-zero.
8. Verify required output files or directories exist and are non-empty before returning their paths.

Conda and Apt instructions belong in the module docstring. Runtime execution uses the executable found on the active process `PATH`; a hidden project Conda wrapper is not part of the standalone contract.

## Galaxy-Compatible Parameter Interface

Galaxy's maintained tool wrappers are the primary reference for deciding which CLI parameters are common enough to expose directly on a ComfyUI node. Parameter names, types, choices, bounds, defaults, and command-line spellings are taken from a pinned revision of the relevant wrapper XML.

Source priority is:

1. the matching wrapper in the [Galaxy IUC tools repository](https://github.com/galaxyproject/tools-iuc);
2. a wrapper maintained in an official Galaxy repository;
3. the tool's own version-matched CLI documentation when no maintained Galaxy wrapper exists.

The mapping follows the Galaxy tool XML `inputs/param` and `argument` semantics described in the [Galaxy tool schema](https://docs.galaxyproject.org/en/master/dev/schema.html). The implementation does not import Galaxy or parse wrapper XML at runtime.

### Required and Optional Inputs

- `required` contains biological input paths, output paths/directories, and tool parameters for which execution has no valid general default.
- `optional` contains commonly adjusted Galaxy parameters with the same default, choices, and valid numeric bounds as the pinned wrapper.
- Every CLI node provides `extra_command` as an `optional` multiline `STRING` with an empty default.
- Parameters omitted by the Galaxy wrapper are available only through `extra_command`, unless exposing one is necessary for a supported output contract.
- Wrapper-only orchestration settings that do not map to the underlying executable are not copied into the node interface.

Each CLI tool file records the wrapper URL, commit, wrapper version, and upstream tool version in its module docstring. The verification matrix records the same source so later wrapper changes do not silently alter existing node behavior.

### Preventing Duplicate Options

Every CLI tool file declares a private managed-option table for the flags already represented by node inputs. Each entry lists all accepted short/long aliases and the number of following CLI values owned by that option. Boolean switches own zero following values.

Before execution, the file:

1. tokenizes `extra_command` with `shlex.split()`;
2. removes managed flags and their owned values;
3. recognizes both `--option value` and `--option=value` forms;
4. recognizes documented short aliases, including attached values when the tool supports them;
5. prints a warning listing ignored managed options; and
6. appends only the remaining tokens to the generated argument list.

Consequently, a value exposed in `required` or `optional` can be changed only through its node input. `extra_command` cannot override it by relying on duplicate-option ordering. Unknown and non-managed CLI options remain unchanged to preserve a CLI-like experience. Shell operators, redirections, and pipelines are passed as ordinary arguments because execution never uses `shell=True`.

The managed-option table is intentionally local to each tool file. A generic cross-tool CLI parser is not introduced because option arity and short-option syntax differ between tools.

## Python Library Contract

- Import optional scientific packages inside the node operation so that a missing package produces a focused installation error without preventing unrelated standalone files from loading.
- Read large FASTA, FASTQ, BAM, VCF, GFF, and matrix inputs through streaming or library-backed file readers.
- Write large derived data to a user-selected output path and return the path.
- Small summaries may additionally return scalar or JSON `STRING` values, but not the complete source dataset.
- Never replace a missing calculation or missing output with random, constant, demonstration, or synthetic results.

## Registration

Each standalone file owns its two ComfyUI mapping dictionaries. Package `__init__.py` files aggregate only E2E-verified modules. `nodes/registry.py` combines those verified dictionaries for the repository entry point but is not required when an individual tool file is copied elsewhere.

Registration counts are derived from mappings and are not asserted as fixed numbers. Removing an invalid node is preferable to retaining it solely to meet a node-count claim.

## Verification

### Dataset Policy

Every E2E test uses an immutable revision of one of these sources:

- the tool or library's official GitHub test/example data;
- an `nf-core/test-datasets` branch or immutable commit intended for that pipeline;
- another official upstream fixture only when neither source above contains the required format.

Generated, handwritten, randomly sampled, or LLM-authored biological fixtures are prohibited. Downloaded fixtures record the upstream URL, revision, license when known, and SHA-256 checksum in a source manifest. Large fixtures remain download-on-demand rather than being committed.

### Checks Per Migrated File

1. **Standalone import:** load the file by path with the repository removed from `sys.path`; assert both mapping dictionaries and the ComfyUI class contract.
2. **Missing dependency:** run with a controlled `PATH` that omits the executable and assert the actionable error.
3. **Input boundary:** pass a nonexistent or wrong-kind path and assert failure occurs before tool execution.
4. **Real E2E:** call the node's declared `FUNCTION` directly using official data and an installed real library or binary.
5. **Artifact validation:** parse outputs independently with the relevant format parser or tool and assert content-level invariants, not merely path existence.
6. **Failure propagation:** where an upstream fixture can trigger a genuine non-zero exit, assert the captured tool error is exposed.
7. **Galaxy parameter parity:** compare every directly exposed CLI parameter with the pinned Galaxy wrapper for name, type, choice set, bounds, default, and emitted flag.
8. **Duplicate filtering:** execute argument construction with managed flags in `extra_command`, covering long, equals, and short forms, and assert the node-input value is the only emitted value.
9. **CLI passthrough:** provide a valid non-managed option through `extra_command` and verify that the real tool receives and applies it.

Tests that inject a dry-run runner, manufacture an expected output, or allow production fallback generation do not count as E2E evidence.

### Verification Matrix

Maintain `docs/node-verification-matrix.md` with one row per registered node and these columns:

| Node | Tool/library and version | Galaxy wrapper URL, revision, and defaults | Official data URL and immutable revision | Input parameters | Required environment | Output validation | Status |
|---|---|---|---|---|---|---|---|

A row may be marked `verified` only when its referenced E2E command has passed in the current environment. Missing proprietary software, reference databases, GPU resources, or licenses are recorded as blockers rather than bypassed with fallback output.

## Implementation Order

1. Add structural tests for the standalone file contract and remove fixed-count assumptions.
2. Add tests for Galaxy-derived node inputs and managed-option filtering.
3. Migrate Biopython into one library file using official Biopython test data.
4. Migrate fastp and FastQC into separate tool files using pinned Galaxy wrappers and nf-core test data.
5. Migrate BWA-MEM2, samtools, and bcftools into separate tool files and validate a small real variant-calling chain.
6. Migrate SPAdes and QUAST, then Kraken2 and Bracken, when their databases and official fixtures are available.
7. Migrate Scanpy and related single-cell Python libraries one library at a time.
8. Migrate remaining genomics, epigenomics, proteomics, microbiome, and structural-biology tools in dependency-sized batches.
9. Replace visualization fallbacks with real file inputs, then migrate their rendering-library files.
10. Remove obsolete shared modules and per-node files only after their final consumer has migrated.
11. Rebuild workflows exclusively from the verified registry.

## Documentation Corrections

README architecture listings and `supplementary_table_s1_node_catalog.md` must be generated from or reconciled with the verified registry. Current statements that all nodes are active, contain zero placeholders, or have subprocess E2E verification must be removed until supported by the verification matrix.

## Completion Criteria

The refactor is complete when:

- every registered implementation follows the one-tool/library-per-file rule;
- every registered file loads standalone and owns both mapping dictionaries;
- no registered node imports a project-local helper;
- no registered node creates fallback biological results;
- every registered CLI node performs PATH discovery, live logging, exit-code handling, and artifact validation;
- every registered CLI node exposes Galaxy-derived common parameters, keeps `extra_command` optional, and filters collisions with managed options;
- every registered node has a verified matrix row backed by an official immutable dataset;
- the complete standalone, E2E, and package-registration suites pass freshly;
- the obsolete files and inaccurate verification claims have been removed.
