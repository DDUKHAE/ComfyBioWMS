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


class BiopythonSeqIOStats:
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


NODE_CLASS_MAPPINGS = {
    "BiopythonSeqIOStats": BiopythonSeqIOStats,
    "BiopythonAlignmentStats": BiopythonAlignmentStats,
}

NODE_DISPLAY_NAME_MAPPINGS = {
    "BiopythonSeqIOStats": "Biopython: Sequence File Statistics",
    "BiopythonAlignmentStats": "Biopython: Alignment Statistics",
}


# Backward compatibility aliases
BiopythonSeqIOStatsNode = BiopythonSeqIOStats
BiopythonAlignmentStatsNode = BiopythonAlignmentStats

__all__ = [
    "BiopythonSeqIOStats",
    "BiopythonSeqIOStatsNode",
    "BiopythonAlignmentStats",
    "BiopythonAlignmentStatsNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
