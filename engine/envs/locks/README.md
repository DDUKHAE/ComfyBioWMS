# Conda Explicit Lock Files for ComfyBIOWMS (osx-arm64)

This directory contains deterministic, explicit package lock files (`conda list --explicit`) generated on macOS 14 (Apple Silicon, osx-arm64).

## Included Environments (Core Benchmark Suite)
The 7 environments required to execute the 5 golden benchmark pipelines and visualizers are fully pinned and deterministic:
1. `bulk_rna_seq_explicit_osx-arm64.txt` (Salmon, Fastp, DESeq2, tximport)
2. `variant_analysis_explicit_osx-arm64.txt` (BWA-MEM2, Samtools, BCFtools)
3. `epigenomics_explicit_osx-arm64.txt` (MACS3, Samtools, DeepTools, BWA-MEM2)
4. `genome_assembly_explicit_osx-arm64.txt` (SPAdes, QUAST, Fastp)
5. `metagenome_explicit_osx-arm64.txt` (Kraken2, Bracken, Fastp)
6. `comfybio_biopython_sequence_ops_explicit_osx-arm64.txt` (Biopython, BLAST, Primer3)
7. `comfybio_publication_visualizer_explicit_osx-arm64.txt` (Matplotlib, Seaborn, PyTorch, Visualizers)

## Extended Domain Environments
The remaining 5 extended domain environments (`cadd_structural`, `longread_genomics`, `microbiome_virome`, `proteomics_metabolomics`, `spatial_single_cell`) have their major dependency versions pinned in `engine/envs/*.yaml` and can be built deterministically across platforms using:
```bash
conda env create -f engine/envs/<env_name>.yaml
```
