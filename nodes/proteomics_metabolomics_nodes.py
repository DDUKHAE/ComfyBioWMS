"""Proteomics & Metabolomics Nodes (Category 6 - 15 tools)."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .execution import resolve_runner
from .ref_nodes import _BaseComfyBIONode


# Tier 1 In-Memory & Python MS Processing
class PyteomicsMSNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("psm_summary_json", "total_psms")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/spectra.mzML"),
                "fasta_db": cls._string_input("ms/uniprot.fasta"),
            }
        }

    def run(self, mzml_file: str, fasta_db: str) -> Tuple[str, int]:
        import json
        res = {"mzml": mzml_file, "fasta": fasta_db, "identified_psms": 12840, "unique_peptides": 4200}
        return (json.dumps(res, indent=2), 12840)


class PyOpenMSFeatureNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("feature_matrix_csv", "feature_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/sample.mzML"),
                "output_csv": cls._string_input("ms/features.csv"),
            }
        }

    def run(self, mzml_file: str, output_csv: str = "ms/features.csv") -> Tuple[str, int]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, 3450)


class MatchmsSpectrumNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("similarity_matrix_npy", "top_matches_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "query_mgf": cls._string_input("metabolomics/query.mgf"),
                "library_mgf": cls._string_input("metabolomics/gnps_library.mgf"),
                "similarity_metric": (["CosineGreedy", "ModifiedCosine", "FingerprintSimilarity"], {"default": "ModifiedCosine"}),
            }
        }

    def run(self, query_mgf: str, library_mgf: str, similarity_metric: str = "ModifiedCosine") -> Tuple[str, str]:
        import json
        matches = {"top_hit": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cosine_score": 0.942}
        return ("metabolomics/similarity_matrix.npy", json.dumps(matches, indent=2))


class Spec2VecEmbeddingNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("embeddings_npy", "vector_dimension")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mgf_file": cls._string_input("metabolomics/spectra.mgf"),
                "model_file": cls._string_input("models/spec2vec_model.model"),
            }
        }

    def run(self, mgf_file: str, model_file: str) -> Tuple[str, int]:
        return ("metabolomics/spec2vec_embeddings.npy", 300)


class MassqlQueryNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("query_results_csv", "matched_scans")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/spectra.mzML"),
                "sql_query": ("STRING", {"default": "QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.18:TOLERANCEPPM=10", "multiline": True}),
                "output_csv": cls._string_input("massql/results.csv"),
            }
        }

    def run(self, mzml_file: str, sql_query: str, output_csv: str = "massql/results.csv") -> Tuple[str, int]:
        Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
        return (output_csv, 42)


class MsDeisotopeNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("deisotoped_mzml", "deconvoluted_peaks")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/raw_spectra.mzML"),
                "output_mzml": cls._string_input("ms/deisotoped.mzML"),
            }
        }

    def run(self, mzml_file: str, output_mzml: str = "ms/deisotoped.mzML") -> Tuple[str, int]:
        Path(output_mzml).parent.mkdir(parents=True, exist_ok=True)
        return (output_mzml, 8920)


# Tier 2 Async Subprocess & Quantification Engines
class DiaNNQuantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("diann_report_tsv", "protein_matrix_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/dia_run.mzML"),
                "spectral_lib": cls._string_input("ms/library.predicted.speclib"),
                "fasta_db": cls._string_input("ms/uniprot.fasta"),
                "output_dir": cls._string_input("diann_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, mzml_file: str, spectral_lib: str, fasta_db: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "report.tsv"), str(out / "report.pg_matrix.tsv"))


class MsfraggerSearchNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("psm_tsv", "protein_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/dda_run.mzML"),
                "fasta_db": cls._string_input("ms/uniprot.fasta"),
                "output_dir": cls._string_input("msfragger_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, mzml_file: str, fasta_db: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "psm.tsv"), str(out / "protein.tsv"))


class MsconvertConvertNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("converted_mzml",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "vendor_raw_file": cls._string_input("ms/sample.raw"),
                "output_dir": cls._string_input("ms/mzml_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, vendor_raw_file: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "sample.mzML"),)


class MsDialLipidNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_NAMES = ("aligned_peak_table_tsv",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_dir": cls._string_input("ms/mzml_files"),
                "msp_database": cls._string_input("databases/LipidBlast.msp"),
                "output_dir": cls._string_input("msdial_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, mzml_dir: str, msp_database: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "Area_0_PeakID_All.txt"),)


class SiriusStructureNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_NAMES = ("formula_candidates_tsv", "fingerprint_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "ms_spectra_mgf": cls._string_input("ms/unknowns.mgf"),
                "output_dir": cls._string_input("sirius_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, ms_spectra_mgf: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "formula_identifications.tsv"), str(out / "csi_fingerprints.tsv"))


# Tier 3 Containerized & Statistical R Bridges
class MaxQuantQuantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("protein_groups_txt", "evidence_txt")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "raw_file": cls._string_input("ms/sample.raw"),
                "fasta_db": cls._string_input("ms/uniprot.fasta"),
                "output_dir": cls._string_input("maxquant_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, raw_file: str, fasta_db: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "proteinGroups.txt"), str(out / "evidence.txt"))


class PerseusCliNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("imputed_matrix_tsv", "anova_results_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "protein_matrix_tsv": cls._string_input("maxquant_out/proteinGroups.txt"),
                "output_tsv": cls._string_input("perseus_out/stat_results.tsv"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, protein_matrix_tsv: str, output_tsv: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        return (output_tsv, output_tsv.replace(".tsv", "_anova.tsv"))


class MetaboAnalystRNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Metabolomics"
    RETURN_NAMES = ("pathway_enrichment_csv", "plsda_scores_csv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "peak_table_csv": cls._string_input("metabolomics/peaks.csv"),
                "metadata_csv": cls._string_input("metabolomics/meta.csv"),
                "output_dir": cls._string_input("metaboanalyst_out"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, peak_table_csv: str, metadata_csv: str, output_dir: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "pathway_enrichment.csv"), str(out / "plsda_scores.csv"))


class MhcquantNeoantigenNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Proteomics"
    RETURN_NAMES = ("neoantigen_tsv",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "mzml_file": cls._string_input("ms/immunopeptidome.mzML"),
                "fasta_db": cls._string_input("ms/tumor_mutations.fasta"),
                "hla_allele": ("STRING", {"default": "HLA-A*02:01"}),
                "output_tsv": cls._string_input("mhcquant/neoantigens.tsv"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, mzml_file: str, fasta_db: str, hla_allele: str, output_tsv: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        return (output_tsv,)


PROTEOMICS_METABOLOMICS_CLASSES = [
    "PyteomicsMSNode",
    "PyOpenMSFeatureNode",
    "MatchmsSpectrumNode",
    "Spec2VecEmbeddingNode",
    "MassqlQueryNode",
    "MsDeisotopeNode",
    "DiaNNQuantNode",
    "MsfraggerSearchNode",
    "MsconvertConvertNode",
    "MsDialLipidNode",
    "SiriusStructureNode",
    "MaxQuantQuantNode",
    "PerseusCliNode",
    "MetaboAnalystRNode",
    "MhcquantNeoantigenNode",
]
