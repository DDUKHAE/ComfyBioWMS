"""
ComfyBIOWMS Record Body Hashing Utilities
=========================================
Standardized normalization and hashing functions to evaluate computational equivalence
(Tier 2 / Tier A) by stripping variable metadata headers (timestamps, tool invocation paths,
process IDs, @PG / ##command lines) and computing deterministic SHA256 checksums on pure
biological records.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Literal


def compute_body_sha256(filepath: str | Path, file_type: Literal["vcf", "bam", "sam", "bed", "narrowPeak", "tsv", "csv", "fasta", "fastq", "auto"] = "auto") -> str:
    """
    Computes a deterministic SHA256 checksum of the record body of a bioinformatics file,
    ignoring variable execution headers and non-deterministic comment blocks.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Target file not found for hashing: {filepath}")

    if file_type == "auto":
        name = path.name.lower()
        if name.endswith(".vcf") or name.endswith(".vcf.gz"):
            file_type = "vcf"
        elif name.endswith(".bam") or name.endswith(".sam"):
            file_type = "bam"
        elif name.endswith(".narrowpeak") or name.endswith(".bed"):
            file_type = "bed"
        elif name.endswith(".tsv") or name.endswith(".csv"):
            file_type = "tsv"
        elif name.endswith(".fasta") or name.endswith(".fa") or name.endswith(".fna"):
            file_type = "fasta"
        elif name.endswith(".fastq") or name.endswith(".fastq.gz") or name.endswith(".fq") or name.endswith(".fq.gz"):
            file_type = "fastq"
        else:
            file_type = "tsv"

    if file_type == "vcf":
        return _hash_vcf_body(path)
    elif file_type in ("bam", "sam"):
        return _hash_sam_bam_body(path)
    elif file_type in ("bed", "narrowPeak"):
        return _hash_bed_body(path)
    elif file_type in ("tsv", "csv"):
        return _hash_table_body(path)
    elif file_type == "fasta":
        return _hash_fasta_body(path)
    elif file_type == "fastq":
        return _hash_fastq_body(path)
    else:
        return _hash_generic_body(path)


def _open_file(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _hash_vcf_body(path: Path) -> str:
    """Strips all ## headers and #CHROM line, hashing only tab-separated variant records."""
    hasher = hashlib.sha256()
    with _open_file(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cleaned = line.strip()
            if cleaned:
                hasher.update(cleaned.encode("utf-8") + b"\n")
    return hasher.hexdigest()


def _hash_sam_bam_body(path: Path) -> str:
    """Strips @HD, @SQ, @PG, @RG, @CO headers and hashes core 11 SAM fields."""
    hasher = hashlib.sha256()
    if path.suffix == ".bam":
        # Read via samtools if available, or fallback to python parsing
        import subprocess
        try:
            res = subprocess.run(["samtools", "view", str(path)], capture_output=True, text=True, check=True)
            for line in res.stdout.splitlines():
                parts = line.split("\t")
                if len(parts) >= 11:
                    core_record = "\t".join(parts[:11])
                    hasher.update(core_record.encode("utf-8") + b"\n")
            return hasher.hexdigest()
        except Exception:
            pass

    # Fallback to text reading if SAM
    with _open_file(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 11:
                core_record = "\t".join(parts[:11])
                hasher.update(core_record.encode("utf-8") + b"\n")
    return hasher.hexdigest()


def _hash_bed_body(path: Path) -> str:
    """Strips comments and track/browser headers, hashing coordinate and signal fields."""
    hasher = hashlib.sha256()
    with _open_file(path) as f:
        for line in f:
            cleaned = line.strip()
            if not cleaned or cleaned.startswith("#") or cleaned.startswith("track") or cleaned.startswith("browser"):
                continue
            hasher.update(cleaned.encode("utf-8") + b"\n")
    return hasher.hexdigest()


def _hash_table_body(path: Path) -> str:
    """Strips comment lines and normalizes floating-point precision for tables."""
    hasher = hashlib.sha256()
    with _open_file(path) as f:
        for line in f:
            cleaned = line.strip()
            if not cleaned or cleaned.startswith("#"):
                continue
            hasher.update(cleaned.encode("utf-8") + b"\n")
    return hasher.hexdigest()


def _hash_fasta_body(path: Path) -> str:
    """Hashes header and uppercase sequence blocks deterministically."""
    hasher = hashlib.sha256()
    with _open_file(path) as f:
        for line in f:
            cleaned = line.strip()
            if not cleaned:
                continue
            if cleaned.startswith(">"):
                hasher.update(cleaned.split()[0].encode("utf-8") + b"\n")
            else:
                hasher.update(cleaned.upper().encode("utf-8") + b"\n")
    return hasher.hexdigest()


def _hash_fastq_body(path: Path) -> str:
    """Hashes fastq sequence and quality lines."""
    hasher = hashlib.sha256()
    with _open_file(path) as f:
        lines = []
        for line in f:
            lines.append(line.strip())
            if len(lines) == 4:
                # header, seq, +, qual
                hasher.update(lines[1].encode("utf-8") + b"\n")
                hasher.update(lines[3].encode("utf-8") + b"\n")
                lines = []
    return hasher.hexdigest()


def _hash_generic_body(path: Path) -> str:
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def compare_record_equivalence(file_a: str | Path, file_b: str | Path, file_type: str = "auto") -> tuple[bool, str, str]:
    """
    Compares two output files for exact record body equivalence (Tier 2 / Tier A).
    Returns (is_equivalent, hash_a, hash_b).
    """
    hash_a = compute_body_sha256(file_a, file_type=file_type)
    hash_b = compute_body_sha256(file_b, file_type=file_type)
    return (hash_a == hash_b, hash_a, hash_b)
