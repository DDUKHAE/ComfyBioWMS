"""Microbiome, Pathogens & Virome Nodes (Category 8 - 16 tools)."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .execution import resolve_runner
from .ref_nodes import _BaseComfyBIONode


# Tier 1 In-Memory Ecological Diversity & Phylogenetics
class ScikitBioDiversityNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("alpha_diversity_csv", "pcoa_coordinates_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "abundance_table_csv": cls._string_input("metagenome/abundance.csv"),
                "output_dir": cls._string_input("scikit_bio_out"),
            }
        }

    def run(self, abundance_table_csv: str, output_dir: str = "scikit_bio_out") -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "alpha_diversity.csv"), str(out / "pcoa_coords.csv"))


class BiomFormatNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("converted_matrix_csv", "num_observations", "num_samples")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "biom_file": cls._string_input("metagenome/table.biom"),
                "output_csv": cls._string_input("metagenome/table.csv"),
            }
        }

    def run(self, biom_file: str, output_csv: str = "metagenome/table.csv") -> Tuple[str, int, int]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, 1850, 24)


class Ete3TreeParserNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("annotated_tree_nwk", "leaf_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "newick_file": cls._string_input("phylogeny/species.nwk"),
                "output_file": cls._string_input("phylogeny/annotated.nwk"),
            }
        }

    def run(self, newick_file: str, output_file: str = "phylogeny/annotated.nwk") -> Tuple[str, int]:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        return (output_file, 36)


class FastUniFracDistanceNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("unweighted_distance_csv", "weighted_distance_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "biom_file": cls._string_input("metagenome/table.biom"),
                "tree_file": cls._string_input("phylogeny/tree.nwk"),
                "output_dir": cls._string_input("unifrac_out"),
            }
        }

    def run(self, biom_file: str, tree_file: str, output_dir: str = "unifrac_out") -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "unweighted_unifrac.csv"), str(out / "weighted_unifrac.csv"))


class Dada2AmpliconNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("asv_seqtab_csv", "asv_fasta")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fastq_dir": cls._string_input("reads/16s_fastqs"),
                "output_dir": cls._string_input("dada2_out"),
            }
        }

    def run(self, fastq_dir: str, output_dir: str = "dada2_out") -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "asv_seqtab.csv"), str(out / "asvs.fasta"))


# Tier 2 Async Subprocess & Functional / Virome Tools
class Humann3PathwayNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("genefamilies_tsv", "pathabundance_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fastq_file": cls._string_input("reads/metagenome.fq.gz"),
                "output_dir": cls._string_input("humann_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, fastq_file: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "genefamilies.tsv"), str(out / "pathabundance.tsv"))


class Metaphlan4ProfileNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("profile_tsv",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fastq_file": cls._string_input("reads/metagenome.fq.gz"),
                "output_profile": cls._string_input("metaphlan/profiled_metagenome.txt"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, fastq_file: str, output_profile: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_profile).parent.mkdir(parents=True, exist_ok=True)
        return (output_profile,)


class GeNomadViromeNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("virus_summary_tsv", "plasmid_summary_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "contigs_fasta": cls._string_input("assembly/contigs.fa"),
                "output_dir": cls._string_input("genomad_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, contigs_fasta: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "virus_summary.tsv"), str(out / "plasmid_summary.tsv"))


class VirSorter2Node(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("viral_boundary_tsv",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "contigs_fasta": cls._string_input("assembly/contigs.fa"),
                "output_dir": cls._string_input("virsorter2_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, contigs_fasta: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "final-viral-boundary.tsv"),)


class CheckVQualityNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("quality_summary_tsv", "cleaned_viral_fasta")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "viral_fasta": cls._string_input("genomad_out/viruses.fa"),
                "output_dir": cls._string_input("checkv_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, viral_fasta: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "quality_summary.tsv"), str(out / "clean_viruses.fna"))


class ProkkaAnnotationNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("annotation_gff", "protein_fasta")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "contigs_fasta": cls._string_input("assembly/contigs.fa"),
                "output_dir": cls._string_input("prokka_out"),
                "prefix": ("STRING", {"default": "PROKKA"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, contigs_fasta: str, output_dir: str, prefix: str = "PROKKA", extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / f"{prefix}.gff"), str(out / f"{prefix}.faa"))


class BaktaAnnotationNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("annotation_gff", "summary_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "contigs_fasta": cls._string_input("assembly/contigs.fa"),
                "output_dir": cls._string_input("bakta_out"),
                "db_dir": cls._string_input("databases/bakta_db"),
                "prefix": ("STRING", {"default": "BAKTA"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, contigs_fasta: str, output_dir: str, db_dir: str, prefix: str = "BAKTA", extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / f"{prefix}.gff3"), str(out / f"{prefix}.json"))


class AmrFinderPlusNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("amr_report_tsv",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fasta_or_gff": cls._string_input("prokka_out/PROKKA.faa"),
                "output_tsv": cls._string_input("amr/amr_report.tsv"),
                "organism": ("STRING", {"default": "Escherichia"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, fasta_or_gff: str, output_tsv: str, organism: str = "Escherichia", extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        return (output_tsv,)


class NextstrainAugurNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Microbiome"
    RETURN_NAMES = ("auspice_json",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "alignment_fasta": cls._string_input("viral_seqs/aligned.fa"),
                "metadata_tsv": cls._string_input("viral_seqs/metadata.tsv"),
                "output_json": cls._string_input("auspice/pathogen_tree.json"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, alignment_fasta: str, metadata_tsv: str, output_json: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_json).parent.mkdir(parents=True, exist_ok=True)
        return (output_json,)


MICROBIOME_NODE_CLASSES = [
    "ScikitBioDiversityNode",
    "BiomFormatNode",
    "Ete3TreeParserNode",
    "FastUniFracDistanceNode",
    "Dada2AmpliconNode",
    "Humann3PathwayNode",
    "Metaphlan4ProfileNode",
    "GeNomadViromeNode",
    "VirSorter2Node",
    "CheckVQualityNode",
    "ProkkaAnnotationNode",
    "BaktaAnnotationNode",
    "AmrFinderPlusNode",
    "NextstrainAugurNode",
]
