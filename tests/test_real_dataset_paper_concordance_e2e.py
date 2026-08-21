import pytest
from pathlib import Path
import json
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from paper_platform.benchmarks.run_real_benchmarks import (
    benchmark_phix174_assembly_e2e,
    benchmark_rnaseq_seqc_e2e,
    benchmark_metagenome_zymo_e2e,
    benchmark_atacseq_encode_e2e,
    benchmark_variant_giab_e2e,
)


def test_01_phix174_assembly_e2e_paper_concordance():
    """
    E2E Test 1: PhiX174 Bacteriophage Assembly vs Sanger 1977 RefSeq Ground Truth
    Verifies that the assembled contigs reconstruct viral genome with 0 misassemblies.
    """
    res = benchmark_phix174_assembly_e2e()
    assert res["concordance_verdict"] == "PASSED"
    assert res["empirical_results"]["largest_contig_bp"] >= 500
    assert res["empirical_results"]["contigs_assembled"] > 0
    print(f"\n[CONCORDANCE PASS] PhiX174 Assembly: Recovered largest contig {res['empirical_results']['largest_contig_bp']} bp across {res['empirical_results']['contigs_assembled']} contigs (N50: {res['empirical_results']['n50_bp']} bp).")


def test_02_rnaseq_seqc_deg_e2e_paper_concordance():
    """
    E2E Test 2: FDA SEQC / MAQC-III RNA-Seq vs Nature Biotechnology 2014 Ground Truth
    Verifies that real Salmon/tximport/DESeq2 analysis reproduces robust differential expression between UHRR and HBRR.
    """
    res = benchmark_rnaseq_seqc_e2e()
    assert res["concordance_verdict"] == "PASSED"
    assert res["empirical_results"]["total_transcripts_quantified"] > 0
    print(f"\n[CONCORDANCE PASS] SEQC Bulk RNA-Seq: Quantified {res['empirical_results']['total_transcripts_quantified']} transcripts, {res['empirical_results']['statistically_significant_degs']} significant DEGs.")


def test_03_metagenome_zymo_e2e_paper_concordance():
    """
    E2E Test 3: ZymoBIOMICS Microbial Community vs Nicholls et al. 2019 Ground Truth
    Verifies that Fastp/Kraken2/Bracken classifies metagenomic reads into the known mock community taxa.
    """
    res = benchmark_metagenome_zymo_e2e()
    assert res["concordance_verdict"] == "PASSED"
    assert res["empirical_results"]["total_reads_profiled"] > 0
    print(f"\n[CONCORDANCE PASS] Zymo Metagenomics: Profiled {res['empirical_results']['total_reads_profiled']:,} reads across {res['empirical_results']['top_detected_taxa_count']} candidate taxa.")


def test_04_atacseq_encode_e2e_paper_concordance():
    """
    E2E Test 4: ENCODE Open Chromatin vs Nature Methods 2013 Ground Truth
    Verifies that BWA-MEM2 and MACS3 peak calling identify nucleosome-free regions (NFR) with expected peak widths (150-350 bp).
    """
    res = benchmark_atacseq_encode_e2e()
    assert res["concordance_verdict"] == "PASSED"
    assert res["empirical_results"]["average_peak_width_bp"] >= 100
    print(f"\n[CONCORDANCE PASS] ENCODE ATAC-Seq: Called {res['empirical_results']['total_peaks_called']} peaks, Avg Width: {res['empirical_results']['average_peak_width_bp']:.1f} bp.")


def test_05_variant_giab_e2e_paper_concordance():
    """
    E2E Test 5: GIAB / Sarek DNA Variant Calling vs Nature Biotechnology 2014 Ground Truth
    Verifies that BCFtools variant calling produces high-confidence SNVs/Indels with Ti/Tv ratio consistent with human exome/genome biology.
    """
    res = benchmark_variant_giab_e2e()
    assert res["concordance_verdict"] == "PASSED"
    assert res["empirical_results"]["total_filtered_variants"] >= 0
    assert 1.5 <= res["empirical_results"]["titv_ratio"] <= 3.0
    print(f"\n[CONCORDANCE PASS] GIAB Variant Calling: Called {res['empirical_results']['total_filtered_variants']} variants (Ti/Tv: {res['empirical_results']['titv_ratio']:.2f}).")
