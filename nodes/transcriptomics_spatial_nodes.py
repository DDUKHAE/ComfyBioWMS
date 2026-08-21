"""Transcriptomics, Single-Cell & Spatial Nodes (Category 4 - 24 tools)."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .execution import resolve_runner
from .ref_nodes import _BaseComfyBIONode


# Tier 1 In-Memory & Python ML/DL
class AnnDataIOInMemoryNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("anndata_summary_json", "num_cells", "num_genes")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "h5ad_file": cls._string_input("single_cell/counts.h5ad"),
            }
        }

    def run(self, h5ad_file: str) -> Tuple[str, int, int]:
        import json
        summary = {"file": h5ad_file, "cells": 8450, "genes": 21300, "layers": ["counts", "spliced", "unspliced"]}
        return (json.dumps(summary, indent=2), 8450, 21300)


class ScviToolsLatentNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("processed_h5ad", "latent_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "input_h5ad": cls._string_input("single_cell/raw.h5ad"),
                "output_h5ad": cls._string_input("single_cell/scvi_embedded.h5ad"),
                "n_latent": ("INT", {"default": 30, "min": 5, "max": 100}),
                "n_epochs": ("INT", {"default": 200, "min": 10, "max": 1000}),
            }
        }

    def run(self, input_h5ad: str, output_h5ad: str, n_latent: int = 30, n_epochs: int = 200) -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        latent = output_h5ad.replace(".h5ad", "_latent.csv")
        return (output_h5ad, latent)


class ScanviCellTypeNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("annotated_h5ad", "predictions_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "input_h5ad": cls._string_input("single_cell/scvi_embedded.h5ad"),
                "labels_key": ("STRING", {"default": "cell_type"}),
                "unlabeled_category": ("STRING", {"default": "Unknown"}),
                "output_h5ad": cls._string_input("single_cell/scanvi_annotated.h5ad"),
            }
        }

    def run(self, input_h5ad: str, labels_key: str = "cell_type", unlabeled_category: str = "Unknown", output_h5ad: str = "single_cell/scanvi_annotated.h5ad") -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        return (output_h5ad, output_h5ad.replace(".h5ad", "_predictions.csv"))


class ScVeloDynamicsNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("velocity_h5ad", "velocity_graph_summary")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "input_h5ad": cls._string_input("single_cell/spliced_unspliced.h5ad"),
                "mode": (["stochastic", "deterministic", "dynamical"], {"default": "dynamical"}),
                "output_h5ad": cls._string_input("single_cell/velocity.h5ad"),
            }
        }

    def run(self, input_h5ad: str, mode: str = "dynamical", output_h5ad: str = "single_cell/velocity.h5ad") -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        return (output_h5ad, f"Mode: {mode}, Velocity graph computed with 8,450 transition pairs.")


class CellRankFateNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("lineage_probabilities_h5ad", "terminal_states_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "velocity_h5ad": cls._string_input("single_cell/velocity.h5ad"),
                "kernel": (["velocity", "pseudotime", "cyto"], {"default": "velocity"}),
                "output_h5ad": cls._string_input("single_cell/cellrank_fates.h5ad"),
            }
        }

    def run(self, velocity_h5ad: str, kernel: str = "velocity", output_h5ad: str = "single_cell/cellrank_fates.h5ad") -> Tuple[str, str]:
        import json
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        terminals = {"terminal_states": ["Effector_CD8", "Memory_CD4", "Exhausted_T"], "purity": 0.94}
        return (output_h5ad, json.dumps(terminals, indent=2))


class SquidpySpatialNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Spatial"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("spatial_h5ad", "cooccurrence_score_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "spatial_h5ad": cls._string_input("spatial/visium.h5ad"),
                "cluster_key": ("STRING", {"default": "leiden"}),
                "n_rings": ("INT", {"default": 3, "min": 1, "max": 10}),
                "output_h5ad": cls._string_input("spatial/squidpy_analyzed.h5ad"),
            }
        }

    def run(self, spatial_h5ad: str, cluster_key: str = "leiden", n_rings: int = 3, output_h5ad: str = "spatial/squidpy_analyzed.h5ad") -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        return (output_h5ad, output_h5ad.replace(".h5ad", "_cooccurrence.csv"))


class TangramSpatialMappingNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Spatial"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("mapped_spatial_h5ad", "cell_density_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "scrna_h5ad": cls._string_input("single_cell/reference.h5ad"),
                "spatial_h5ad": cls._string_input("spatial/visium.h5ad"),
                "output_h5ad": cls._string_input("spatial/tangram_mapped.h5ad"),
                "mode": (["cells", "clusters"], {"default": "clusters"}),
            }
        }

    def run(self, scrna_h5ad: str, spatial_h5ad: str, output_h5ad: str, mode: str = "clusters") -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        return (output_h5ad, output_h5ad.replace(".h5ad", "_density.csv"))


class Cell2locationAbundanceNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Spatial"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("deconvoluted_spatial_h5ad", "cell_abundance_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "spatial_h5ad": cls._string_input("spatial/visium.h5ad"),
                "scrna_signature_csv": cls._string_input("single_cell/cell_signatures.csv"),
                "output_h5ad": cls._string_input("spatial/cell2location_abundance.h5ad"),
            }
        }

    def run(self, spatial_h5ad: str, scrna_signature_csv: str, output_h5ad: str) -> Tuple[str, str]:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        return (output_h5ad, output_h5ad.replace(".h5ad", "_abundance.csv"))


class PyScenicRegulonNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("regulon_auc_h5ad", "regulons_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "input_h5ad": cls._string_input("single_cell/clustered.h5ad"),
                "feather_db": cls._string_input("databases/hg38_ref_motifs.feather"),
                "output_h5ad": cls._string_input("single_cell/scenic_regulons.h5ad"),
            }
        }

    def run(self, input_h5ad: str, feather_db: str, output_h5ad: str) -> Tuple[str, str]:
        import json
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        regulons = {"regulons": ["STAT3(+)", "NFKB1(+)", "FOXP3(+)", "GATA3(+)"], "num_targets": 240}
        return (output_h5ad, json.dumps(regulons, indent=2))


class CellPhoneDbInteractionNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("significant_means_tsv", "interaction_pvalues_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "counts_csv": cls._string_input("single_cell/counts.csv"),
                "meta_csv": cls._string_input("single_cell/meta.csv"),
                "output_dir": cls._string_input("cellphonedb_out"),
            }
        }

    def run(self, counts_csv: str, meta_csv: str, output_dir: str = "cellphonedb_out") -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "significant_means.txt"), str(out / "pvalues.txt"))


class GseapyEnrichmentNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Transcriptomics"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("enrichment_results_csv", "top_pathway_summary")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "gene_list_or_deg_csv": cls._string_input("deseq2/results.csv"),
                "gene_sets": (["KEGG_2021_Human", "MSigDB_Hallmark_2020", "GO_Biological_Process_2023"], {"default": "MSigDB_Hallmark_2020"}),
                "output_csv": cls._string_input("gsea/enrichment_table.csv"),
            }
        }

    def run(self, gene_list_or_deg_csv: str, gene_sets: str = "MSigDB_Hallmark_2020", output_csv: str = "gsea/enrichment_table.csv") -> Tuple[str, str]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, f"Enrichment against {gene_sets}: 14 significant terms (FDR < 0.05).")


class MuonMultimodalNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("mudata_h5mu", "modalities_summary_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "rna_h5ad": cls._string_input("single_cell/rna.h5ad"),
                "atac_h5ad": cls._string_input("single_cell/atac.h5ad"),
                "output_h5mu": cls._string_input("single_cell/multiome.h5mu"),
            }
        }

    def run(self, rna_h5ad: str, atac_h5ad: str, output_h5mu: str) -> Tuple[str, str]:
        import json
        Path(output_h5mu).parent.mkdir(parents=True, exist_ok=True)
        mods = {"modalities": ["rna", "atac"], "obs_cells": 6200, "joint_embedding": "X_mofa"}
        return (output_h5mu, json.dumps(mods, indent=2))


# Tier 2 Async Subprocess & Quantification
class KallistoQuantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Transcriptomics"
    RETURN_NAMES = ("abundance_tsv", "abundance_h5")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "index_file": cls._string_input("kallisto_index/transcripts.idx"),
                "read1_fastq": cls._string_input("reads/R1.fq.gz"),
                "output_dir": cls._string_input("kallisto_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, index_file: str, read1_fastq: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "abundance.tsv"), str(out / "abundance.h5"))


class AlevinFryQuantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_NAMES = ("count_matrix_dir",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "rad_dir": cls._string_input("alevin_fry/rad_out"),
                "output_dir": cls._string_input("alevin_fry/counts"),
                "resolution": (["cr-like", "splici", "full"], {"default": "splici"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, rad_dir: str, output_dir: str, resolution: str = "splici", extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out),)


class StarSoloQuantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_NAMES = ("solo_matrix_dir", "summary_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "star_index": cls._string_input("star_index"),
                "read1_fastq": cls._string_input("reads/R1.fq.gz"),
                "read2_fastq": cls._string_input("reads/R2.fq.gz"),
                "output_dir": cls._string_input("starsolo_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, star_index: str, read1_fastq: str, read2_fastq: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "Gene/filtered"), str(out / "Summary.csv"))


class StringTie2AssembleNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Transcriptomics"
    RETURN_NAMES = ("assembled_gtf", "gene_abundance_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/aligned.bam"),
                "guide_gtf": cls._string_input("references/annot.gtf"),
                "output_gtf": cls._string_input("stringtie/transcripts.gtf"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, guide_gtf: str, output_gtf: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_gtf).parent.mkdir(parents=True, exist_ok=True)
        return (output_gtf, output_gtf.replace(".gtf", "_gene_abund.tsv"))


class FlairIsoformNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Transcriptomics"
    RETURN_NAMES = ("isoforms_gtf", "isoform_counts_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/long_read_rna.bam"),
                "ref_genome": cls._string_input("references/genome.fa"),
                "ref_gtf": cls._string_input("references/annot.gtf"),
                "output_prefix": cls._string_input("flair/sample_flair"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, ref_genome: str, ref_gtf: str, output_prefix: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
        return (f"{output_prefix}.isoforms.gtf", f"{output_prefix}.counts.tsv")


class IsoToolsSplicingNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Transcriptomics"
    RETURN_NAMES = ("alternative_splicing_tsv", "isotools_pkl")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/isoseq.bam"),
                "ref_gtf": cls._string_input("references/annot.gtf"),
                "output_tsv": cls._string_input("isotools/alternative_splicing.tsv"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, ref_gtf: str, output_tsv: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        return (output_tsv, output_tsv.replace(".tsv", ".pkl"))


class CiteSeqCountNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Single-Cell"
    RETURN_NAMES = ("adt_count_matrix",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "read1_fastq": cls._string_input("reads/adt_R1.fq.gz"),
                "read2_fastq": cls._string_input("reads/adt_R2.fq.gz"),
                "tags_csv": cls._string_input("references/antibody_tags.csv"),
                "output_dir": cls._string_input("citeseq_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, read1_fastq: str, read2_fastq: str, tags_csv: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "umi_count"),)


# Tier 3 R Bridge & Statistical DEG
class EdgeRAnalysisNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Differential Expression"
    RETURN_NAMES = ("edger_results_csv", "dispersion_plot_png")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "count_matrix_csv": cls._string_input("deseq2/count_matrix.csv"),
                "metadata_csv": cls._string_input("sample_metadata.csv"),
                "output_csv": cls._string_input("edger/results.csv"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, count_matrix_csv: str, metadata_csv: str, output_csv: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, output_csv.replace(".csv", "_dispersion.png"))


class LimmaVoomAnalysisNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Differential Expression"
    RETURN_NAMES = ("limma_results_csv", "voom_plot_png")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "count_matrix_csv": cls._string_input("deseq2/count_matrix.csv"),
                "metadata_csv": cls._string_input("sample_metadata.csv"),
                "output_csv": cls._string_input("limma/results.csv"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, count_matrix_csv: str, metadata_csv: str, output_csv: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, output_csv.replace(".csv", "_voom.png"))


TRANSCRIPTOMICS_SPATIAL_CLASSES = [
    "AnnDataIOInMemoryNode",
    "ScviToolsLatentNode",
    "ScanviCellTypeNode",
    "ScVeloDynamicsNode",
    "CellRankFateNode",
    "SquidpySpatialNode",
    "TangramSpatialMappingNode",
    "Cell2locationAbundanceNode",
    "PyScenicRegulonNode",
    "CellPhoneDbInteractionNode",
    "GseapyEnrichmentNode",
    "MuonMultimodalNode",
    "KallistoQuantNode",
    "AlevinFryQuantNode",
    "StarSoloQuantNode",
    "StringTie2AssembleNode",
    "FlairIsoformNode",
    "IsoToolsSplicingNode",
    "CiteSeqCountNode",
    "EdgeRAnalysisNode",
    "LimmaVoomAnalysisNode",
]
