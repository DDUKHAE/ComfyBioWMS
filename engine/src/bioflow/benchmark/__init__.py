"""
ComfyBIOWMS Benchmark Package
=============================
Provides graph integrity checks, 3-Tier concordance benchmarks, deterministic record body
hashing, and negative control verification harnesses.
"""

from .graph_checks import check_dag
from .hashing import compute_body_sha256, compare_record_equivalence
from .negative_controls import (
    run_assembly_negative_control,
    run_rnaseq_negative_control,
    run_variant_negative_control,
    run_atacseq_negative_control,
    run_metagenome_negative_control,
)
from .concordance import (
    benchmark_phix174_assembly_e2e,
    benchmark_rnaseq_seqc_e2e,
    benchmark_metagenome_zymo_e2e,
    benchmark_atacseq_encode_e2e,
    benchmark_variant_giab_e2e,
)

__all__ = [
    "check_dag",
    "compute_body_sha256",
    "compare_record_equivalence",
    "run_assembly_negative_control",
    "run_rnaseq_negative_control",
    "run_variant_negative_control",
    "run_atacseq_negative_control",
    "run_metagenome_negative_control",
    "benchmark_phix174_assembly_e2e",
    "benchmark_rnaseq_seqc_e2e",
    "benchmark_metagenome_zymo_e2e",
    "benchmark_atacseq_encode_e2e",
    "benchmark_variant_giab_e2e",
]
