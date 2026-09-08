"""Biopython nodes.

Python packages: biopython==1.88
External binaries: none
"""

import json
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and str(custom_dir).strip():
        out = Path(custom_dir).expanduser().resolve()
    else:
        try:
            folder_paths = __import__("folder_paths")
            base = Path(folder_paths.get_output_directory())
        except Exception:
            base = Path.cwd() / "ComfyUI" / "output"
        out = base / node_name
    return out


def _get_sequence(sequence: str, sequence_file: str, label: str = "Sequence") -> str:
    if sequence_file and str(sequence_file).strip():
        path = _file(sequence_file, f"{label} file")
        from Bio import SeqIO
        records = list(SeqIO.parse(str(path), "fasta" if path.suffix.lower() in [".fa", ".fasta", ".fna"] else ("fastq" if path.suffix.lower() in [".fq", ".fastq"] else "fasta")))
        if not records:
            return path.read_text(encoding="utf-8", errors="ignore").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
        return str(records[0].seq).strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    clean = sequence.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if not clean:
        raise ValueError(f"{label} is empty. Provide a valid sequence string or sequence_file.")
    return clean


class BiopythonSeqIOStats:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("summary_json", "sequence_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "sequence_file": ("STRING", {"default": ""}),
                "file_format": (
                    ["fasta", "fastq", "genbank", "embl"],
                    {"default": "fasta"},
                ),
            }
        }

    def run(self, sequence_file: str, file_format: str = "fasta"):
        from Bio import SeqIO

        rows = []
        for record in SeqIO.parse(str(_file(sequence_file, "Sequence input")), file_format):
            sequence = str(record.seq).upper()
            gc = sequence.count("G") + sequence.count("C")
            rows.append(
                {
                    "id": record.id,
                    "description": record.description,
                    "length": len(sequence),
                    "gc_percent": round(100 * gc / len(sequence), 6) if sequence else 0.0,
                }
            )
        return json.dumps(rows), len(rows)


class BiopythonAlignmentStats:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT", "FLOAT")
    RETURN_NAMES = ("summary_json", "rows", "columns", "mean_pairwise_identity")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "alignment_file": ("STRING", {"default": ""}),
                "file_format": (
                    ["clustal", "fasta", "phylip", "stockholm"],
                    {"default": "clustal"},
                ),
            }
        }

    def run(self, alignment_file: str, file_format: str = "clustal"):
        from Bio import AlignIO

        alignment = AlignIO.read(str(_file(alignment_file, "Alignment input")), file_format)
        identities = []
        for left, right in combinations(alignment, 2):
            pairs = [(a, b) for a, b in zip(left.seq, right.seq) if a != "-" and b != "-"]
            if pairs:
                identities.append(sum(a == b for a, b in pairs) / len(pairs))
        value = round(100 * sum(identities) / len(identities), 6) if identities else 100.0
        summary = {
            "rows": len(alignment),
            "columns": alignment.get_alignment_length(),
            "mean_pairwise_identity": value,
        }
        return json.dumps(summary), summary["rows"], summary["columns"], value


class BiopythonSeqTransform:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = ("reverse_complement", "complement", "rna_transcribe", "protein_translate")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "sequence": ("STRING", {"default": "ATGCGATCGATCGATCGATAG", "multiline": True}),
            },
            "optional": {
                "sequence_file": ("STRING", {"default": ""}),
                "genetic_code_table": ("INT", {"default": 1, "min": 1, "max": 33}),
                "to_stop": ("BOOLEAN", {"default": False}),
            },
        }

    def run(
        self,
        sequence: str = "ATGCGATCGATCGATCGATAG",
        sequence_file: str = "",
        genetic_code_table: int = 1,
        to_stop: bool = False,
        **kwargs,
    ):
        from Bio.Seq import Seq

        dna_clean = _get_sequence(sequence, sequence_file, "DNA/RNA Sequence")
        s = Seq(dna_clean)

        rev_comp = str(s.reverse_complement())
        comp = str(s.complement())
        rna = str(s.transcribe())

        try:
            protein = str(s.translate(table=genetic_code_table, to_stop=to_stop))
        except Exception:
            trimmed_len = (len(dna_clean) // 3) * 3
            if trimmed_len > 0:
                protein = str(s[:trimmed_len].translate(table=genetic_code_table, to_stop=to_stop))
            else:
                protein = ""

        return rev_comp, comp, rna, protein


class BiopythonGCContent:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("FLOAT", "FLOAT", "INT", "FLOAT", "FLOAT", "STRING")
    RETURN_NAMES = ("gc_percent", "at_percent", "length", "molecular_weight", "melting_temp_c", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "sequence": ("STRING", {"default": "ATGCGATCGATCGATCGATAG", "multiline": True}),
            },
            "optional": {
                "sequence_file": ("STRING", {"default": ""}),
            },
        }

    def run(
        self,
        sequence: str = "ATGCGATCGATCGATCGATAG",
        sequence_file: str = "",
        **kwargs,
    ):
        from Bio.Seq import Seq
        from Bio.SeqUtils import gc_fraction, molecular_weight

        seq_clean = _get_sequence(sequence, sequence_file, "Sequence")
        s = Seq(seq_clean)
        seq_len = len(seq_clean)

        gc_pct = round(gc_fraction(s) * 100.0, 4) if seq_len > 0 else 0.0
        at_count = seq_clean.count("A") + seq_clean.count("T") + seq_clean.count("U")
        at_pct = round(100.0 * at_count / seq_len, 4) if seq_len > 0 else 0.0

        is_rna = "U" in seq_clean and "T" not in seq_clean
        seq_type = "RNA" if is_rna else "DNA"
        try:
            mw = round(float(molecular_weight(s, seq_type=seq_type)), 2)
        except Exception:
            mw = 0.0

        try:
            from Bio.SeqUtils import MeltingTemp as mt
            tm = round(float(mt.Tm_NN(s)), 2)
        except Exception:
            gc_count = seq_clean.count("G") + seq_clean.count("C")
            tm = round(float(2 * at_count + 4 * gc_count), 2)

        summary = {
            "length": seq_len,
            "gc_percent": gc_pct,
            "at_percent": at_pct,
            "molecular_weight_da": mw,
            "melting_temperature_c": tm,
            "counts": {
                "A": seq_clean.count("A"),
                "C": seq_clean.count("C"),
                "G": seq_clean.count("G"),
                "T": seq_clean.count("T"),
                "U": seq_clean.count("U"),
                "N": seq_clean.count("N"),
            },
            "seq_type": seq_type,
        }

        return gc_pct, at_pct, seq_len, mw, tm, json.dumps(summary, indent=2)


class BiopythonPairwiseAlign:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "FLOAT", "FLOAT")
    RETURN_NAMES = ("alignment_text", "alignment_score", "percent_identity")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "seq1": ("STRING", {"default": "ATGCGATCGATCG", "multiline": True}),
                "seq2": ("STRING", {"default": "ATGCGATAGATCG", "multiline": True}),
            },
            "optional": {
                "mode": (["global", "local"], {"default": "global"}),
                "match_score": ("FLOAT", {"default": 1.0, "step": 0.1}),
                "mismatch_score": ("FLOAT", {"default": -1.0, "step": 0.1}),
                "open_gap_score": ("FLOAT", {"default": -2.5, "step": 0.1}),
                "extend_gap_score": ("FLOAT", {"default": -0.5, "step": 0.1}),
            },
        }

    def run(
        self,
        seq1: str = "ATGCGATCGATCG",
        seq2: str = "ATGCGATAGATCG",
        mode: str = "global",
        match_score: float = 1.0,
        mismatch_score: float = -1.0,
        open_gap_score: float = -2.5,
        extend_gap_score: float = -0.5,
        **kwargs,
    ):
        from Bio.Align import PairwiseAligner

        s1 = seq1.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
        s2 = seq2.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
        if not s1 or not s2:
            raise ValueError("Both seq1 and seq2 must be non-empty for pairwise alignment.")

        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score

        alignments = aligner.align(s1, s2)
        score = round(float(alignments.score), 4)
        alignment_str = str(alignments[0]) if len(alignments) > 0 else "No alignment found"

        matches = sum(1 for a, b in zip(s1, s2) if a == b)
        pid = round(100.0 * matches / max(len(s1), len(s2), 1), 2)

        return alignment_str, score, pid


class BiopythonProtParam:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("FLOAT", "FLOAT", "FLOAT", "FLOAT", "STRING")
    RETURN_NAMES = ("molecular_weight", "isoelectric_point", "instability_index", "gravy", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "protein_sequence": ("STRING", {"default": "MSIVMGRKGAR", "multiline": True}),
            },
            "optional": {},
        }

    def run(self, protein_sequence: str = "MSIVMGRKGAR", **kwargs):
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        clean_prot = "".join(c for c in protein_sequence.strip().upper() if c.isalpha())
        if not clean_prot:
            raise ValueError("protein_sequence must contain valid amino acid characters.")

        pa = ProteinAnalysis(clean_prot)
        mw = round(float(pa.molecular_weight()), 2)
        pi = round(float(pa.isoelectric_point()), 2)
        instability = round(float(pa.instability_index()), 2)
        gravy = round(float(pa.gravy()), 4)
        aromaticity = round(float(pa.aromaticity()), 4)
        sec_struct = pa.secondary_structure_fraction()

        summary = {
            "length": len(clean_prot),
            "molecular_weight_da": mw,
            "isoelectric_point_pi": pi,
            "instability_index": instability,
            "is_stable": instability < 40.0,
            "gravy_hydropathicity": gravy,
            "aromaticity": aromaticity,
            "secondary_structure_fractions": {
                "helix": round(float(sec_struct[0]), 4),
                "turn": round(float(sec_struct[1]), 4),
                "sheet": round(float(sec_struct[2]), 4),
            },
            "amino_acid_counts": pa.count_amino_acids(),
        }

        return mw, pi, instability, gravy, json.dumps(summary, indent=2)


class BiopythonRestrictionDigest:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "INT")
    RETURN_NAMES = ("cut_sites_json", "fragment_lengths_json", "total_cut_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "dna_sequence": ("STRING", {"default": "GAATTCGGATCCAAGCTTGCGGCCGCCTCGAG", "multiline": True}),
            },
            "optional": {
                "enzymes": ("STRING", {"default": "EcoRI,BamHI,HindIII,NotI,XhoI"}),
                "linear": ("BOOLEAN", {"default": True}),
            },
        }

    def run(
        self,
        dna_sequence: str = "GAATTCGGATCCAAGCTTGCGGCCGCCTCGAG",
        enzymes: str = "EcoRI,BamHI,HindIII,NotI,XhoI",
        linear: bool = True,
        **kwargs,
    ):
        from Bio.Seq import Seq
        from Bio.Restriction import RestrictionBatch

        clean_dna = dna_sequence.strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
        if not clean_dna:
            raise ValueError("dna_sequence cannot be empty.")

        enzyme_list = [e.strip() for e in enzymes.split(",") if e.strip()]
        if not enzyme_list:
            raise ValueError("No valid enzyme names specified in enzymes parameter.")

        rb = RestrictionBatch(enzyme_list)
        results = rb.search(Seq(clean_dna), linear=linear)

        cut_sites = {str(k): list(v) for k, v in results.items()}
        all_cuts = sorted(list(set(c for cuts in cut_sites.values() for c in cuts)))

        seq_len = len(clean_dna)
        if all_cuts:
            cut_points = [0] + all_cuts + [seq_len]
            fragments = [cut_points[i + 1] - cut_points[i] for i in range(len(cut_points) - 1)]
        else:
            fragments = [seq_len]

        return json.dumps(cut_sites, indent=2), json.dumps(fragments), len(all_cuts)


class BiopythonSeqFilter:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Biopython"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("filtered_sequence_file", "passed_count", "total_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "sequence_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "min_length": ("INT", {"default": 100, "min": 0, "max": 10000000}),
                "max_length": ("INT", {"default": 1000000, "min": 0, "max": 10000000}),
                "min_gc": ("FLOAT", {"default": 0.0, "min": 0.0, "max": 100.0, "step": 0.1}),
                "max_gc": ("FLOAT", {"default": 100.0, "min": 0.0, "max": 100.0, "step": 0.1}),
                "file_format": (["fasta", "fastq"], {"default": "fasta"}),
            },
        }

    def run(
        self,
        sequence_file: str,
        min_length: int = 100,
        max_length: int = 1000000,
        min_gc: float = 0.0,
        max_gc: float = 100.0,
        file_format: str = "fasta",
        output_dir: str = "",
        **kwargs,
    ):
        from Bio import SeqIO

        seq_path = _file(sequence_file, "Sequence file to filter")
        out_base = _output_dir("BiopythonSeqFilter", output_dir)
        out_base.mkdir(parents=True, exist_ok=True)

        passed_records = []
        total_records = 0
        for record in SeqIO.parse(str(seq_path), file_format):
            total_records += 1
            seq_str = str(record.seq).upper()
            l = len(seq_str)
            if not (min_length <= l <= max_length):
                continue
            gc = 100.0 * (seq_str.count("G") + seq_str.count("C")) / l if l > 0 else 0.0
            if not (min_gc <= gc <= max_gc):
                continue
            passed_records.append(record)

        out_file = out_base / f"filtered_{seq_path.name}"
        SeqIO.write(passed_records, str(out_file), file_format)

        return str(out_file), len(passed_records), total_records


NODE_CLASS_MAPPINGS = {
    "BiopythonSeqIOStats": BiopythonSeqIOStats,
    "BiopythonAlignmentStats": BiopythonAlignmentStats,
    "BiopythonSeqTransform": BiopythonSeqTransform,
    "BiopythonGCContent": BiopythonGCContent,
    "BiopythonPairwiseAlign": BiopythonPairwiseAlign,
    "BiopythonProtParam": BiopythonProtParam,
    "BiopythonRestrictionDigest": BiopythonRestrictionDigest,
    "BiopythonSeqFilter": BiopythonSeqFilter,
}

NODE_DISPLAY_NAME_MAPPINGS = {
    "BiopythonSeqIOStats": "Biopython: Sequence File Statistics",
    "BiopythonAlignmentStats": "Biopython: Alignment Statistics",
    "BiopythonSeqTransform": "Biopython: Sequence Transformation (RevComp/Transcribe/Translate)",
    "BiopythonGCContent": "Biopython: GC Content & Sequence Metrics",
    "BiopythonPairwiseAlign": "Biopython: Pairwise Sequence Alignment",
    "BiopythonProtParam": "Biopython: Protein Physicochemical Properties",
    "BiopythonRestrictionDigest": "Biopython: Restriction Enzyme Digestion",
    "BiopythonSeqFilter": "Biopython: Sequence Length & GC Filter",
}

# Backward compatibility aliases
BiopythonSeqIOStatsNode = BiopythonSeqIOStats
BiopythonAlignmentStatsNode = BiopythonAlignmentStats
BiopythonSeqTransformNode = BiopythonSeqTransform
BiopythonSeqTranscribeNode = BiopythonSeqTransform
BiopythonGCContentNode = BiopythonGCContent
BiopythonPairwiseAlignNode = BiopythonPairwiseAlign
BiopythonAlignIONode = BiopythonPairwiseAlign
BiopythonProtParamNode = BiopythonProtParam
BiopythonRestrictionDigestNode = BiopythonRestrictionDigest
BiopythonRestrictionNode = BiopythonRestrictionDigest
BiopythonSeqFilterNode = BiopythonSeqFilter

__all__ = [
    "BiopythonSeqIOStats",
    "BiopythonAlignmentStats",
    "BiopythonSeqTransform",
    "BiopythonGCContent",
    "BiopythonPairwiseAlign",
    "BiopythonProtParam",
    "BiopythonRestrictionDigest",
    "BiopythonSeqFilter",
    "BiopythonSeqIOStatsNode",
    "BiopythonAlignmentStatsNode",
    "BiopythonSeqTransformNode",
    "BiopythonSeqTranscribeNode",
    "BiopythonGCContentNode",
    "BiopythonPairwiseAlignNode",
    "BiopythonAlignIONode",
    "BiopythonProtParamNode",
    "BiopythonRestrictionDigestNode",
    "BiopythonRestrictionNode",
    "BiopythonSeqFilterNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
