"""
Unit and Fault-Injection Tests for Benchmark Comparator.
Verifies all 10 fault-injection test scenarios mandated by Section 5 of docs/review-remediation-guideline.md.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from bioflow.benchmark.comparator import compare_tables, compare_metagenome_profiles


@pytest.fixture
def baseline_table():
    return pd.DataFrame({
        "gene_id": [f"GENE_{i}" for i in range(1, 11)],
        "log2FoldChange": [0.5, -1.2, 0.0, 2.3, -0.4, 1.1, -2.0, 0.8, -0.1, 3.2],
        "pvalue": [0.01, 0.005, 0.9, 0.0001, 0.4, 0.02, 0.00001, 0.08, 0.85, 0.00005],
        "sample1": [100, 250, 10, 450, 80, 120, 5, 200, 30, 600],
    })


def test_fault_01_identical_tables(baseline_table):
    """Case 1: Exactly identical tables -> must PASS."""
    res = compare_tables(baseline_table, baseline_table.copy(), id_col="gene_id", value_cols=["log2FoldChange", "pvalue"])
    assert res["status"] == "PASS"
    assert res["max_absolute_difference"] == 0.0
    assert res["exceeding_tolerance_count"] == 0


def test_fault_02_row_order_permutation(baseline_table):
    """Case 2: Reordered rows -> must PASS with identical diffs after ID sorting."""
    shuffled = baseline_table.sample(frac=1.0, random_state=42).reset_index(drop=True)
    res = compare_tables(baseline_table, shuffled, id_col="gene_id", value_cols=["log2FoldChange", "pvalue"])
    assert res["status"] == "PASS"
    assert res["max_absolute_difference"] == 0.0


def test_fault_03_missing_or_added_feature(baseline_table):
    """Case 3: Missing or added feature -> must report ID mismatch and FAIL."""
    # Drop one feature
    dropped = baseline_table.iloc[:-1].copy()
    res = compare_tables(baseline_table, dropped, id_col="gene_id", value_cols=["log2FoldChange"])
    assert res["status"] == "FAIL"
    assert "Feature ID set mismatch" in res["reason"]
    assert "GENE_10" in res["left_only_features"]

    # Add extra feature
    extra = pd.concat([baseline_table, pd.DataFrame([{"gene_id": "GENE_EXTRA", "log2FoldChange": 1.0, "pvalue": 0.05, "sample1": 50}])], ignore_index=True)
    res2 = compare_tables(baseline_table, extra, id_col="gene_id", value_cols=["log2FoldChange"])
    assert res2["status"] == "FAIL"
    assert "GENE_EXTRA" in res2["right_only_features"]


def test_fault_04_sample_column_mismatch(baseline_table):
    """Case 4: Perturbed sample column names -> sample correspondence check must FAIL."""
    renamed = baseline_table.rename(columns={"sample1": "sample_WRONG"})
    res = compare_tables(baseline_table, renamed, id_col="gene_id", value_cols=["log2FoldChange"], sample_cols=["sample1"])
    assert res["status"] == "FAIL"
    assert "Sample column mismatch" in res["reason"]


def test_fault_05_numerical_divergence_exceeding_tolerance(baseline_table):
    """Case 5: Numeric divergence > tolerance -> must FAIL."""
    perturbed = baseline_table.copy()
    perturbed.loc[0, "log2FoldChange"] += 0.5  # Large deviation
    res = compare_tables(baseline_table, perturbed, id_col="gene_id", value_cols=["log2FoldChange"], atol=1e-4, rtol=1e-4)
    assert res["status"] == "FAIL"
    assert "Numerical divergence exceeds tolerance" in res["reason"]
    assert res["exceeding_tolerance_count"] >= 1


def test_fault_06_one_sided_na(baseline_table):
    """Case 6: One table has NA while other has finite value -> must FAIL on NA mismatch."""
    with_na = baseline_table.copy()
    with_na.loc[2, "log2FoldChange"] = np.nan
    res = compare_tables(baseline_table, with_na, id_col="gene_id", value_cols=["log2FoldChange"])
    assert res["status"] == "FAIL"
    assert "NA state mismatch" in res["reason"]


def test_fault_07_all_values_na(baseline_table):
    """Case 7: All values are NA -> must return NOT_EVALUATED and strictly forbid PASS."""
    all_na1 = baseline_table.copy()
    all_na2 = baseline_table.copy()
    all_na1["log2FoldChange"] = np.nan
    all_na2["log2FoldChange"] = np.nan

    res = compare_tables(all_na1, all_na2, id_col="gene_id", value_cols=["log2FoldChange"])
    assert res["status"] == "NOT_EVALUATED"
    assert res["valid_numeric_comparisons"] == 0
    assert res["status"] != "PASS"


def test_fault_08_identical_file_independent_mode(tmp_path, baseline_table):
    """Case 8: Passing exact same file path in independent replay mode -> must be REJECTED/FAIL."""
    test_csv = tmp_path / "same_file.csv"
    baseline_table.to_csv(test_csv, index=False)

    res = compare_tables(
        test_csv, test_csv,
        id_col="gene_id",
        value_cols=["log2FoldChange"],
        independent_mode=True,
    )
    assert res["status"] == "FAIL"
    assert res["identical_file_detected"] is True
    assert "Self-comparison detected" in res["reason"]


def test_fault_09_duplicate_feature_ids():
    """Case 9: Duplicate feature IDs -> must not silently aggregate, must report error/FAIL."""
    df_dup = pd.DataFrame({
        "gene_id": ["GENE_1", "GENE_1", "GENE_2"],
        "val": [1.0, 2.0, 3.0],
    })
    df_clean = pd.DataFrame({
        "gene_id": ["GENE_1", "GENE_2"],
        "val": [1.0, 3.0],
    })
    res = compare_tables(df_dup, df_clean, id_col="gene_id", value_cols=["val"])
    assert res["status"] == "FAIL"
    assert "Duplicate feature IDs" in res["reason"]


def test_fault_10_metagenome_unexpected_taxon():
    """Case 10: Unexpected taxon outside expected list -> must be counted as FP and in TVD."""
    df_expected = pd.DataFrame({
        "name": ["Species_A", "Species_B"],
        "expected_fraction": [0.6, 0.4],
    })
    # Observed includes unexpected Species_C (false positive contamination / misclassification)
    df_observed = pd.DataFrame({
        "name": ["Species_A", "Species_B", "Species_C"],
        "fraction": [0.5, 0.3, 0.2],
    })

    res = compare_metagenome_profiles(df_observed, df_expected, taxon_col="name", abund_col="fraction")
    assert res["status"] == "COMPOSITION_DISCREPANCY"
    assert res["false_positives"] == 1
    assert "Species_C" in res["unexpected_taxa"]
    assert res["precision"] < 1.0
    assert res["tvd"] > 0.0
    # Expected fraction for Species_C should have been treated as 0.0
    merged = res["merged_profile"]
    c_row = merged[merged["name"] == "Species_C"].iloc[0]
    assert c_row["expected_fraction"] == 0.0
    assert c_row["norm_observed"] > 0.0
