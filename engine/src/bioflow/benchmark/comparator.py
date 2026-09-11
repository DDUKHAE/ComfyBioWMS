"""
ComfyBIOWMS Scientific Benchmark Comparator
============================================
Implements rigorous tabular and profile comparison with strict verification ordering:
1. File/DataFrame existence and schema
2. File path / inode identity guard (rejecting self-comparison in independent replay mode)
3. Feature ID uniqueness (detecting duplicate IDs)
4. Total ID set completeness (reporting left-only and right-only via outer join, avoiding silent inner-join drops)
5. Sample column correspondence
6. NA / Inf state consistency (strictly returning NOT_EVALUATED when zero valid entries exist, never false PASS)
7. Unrounded numerical comparison with explicit rtol / atol thresholds
8. Dynamic boolean verdict generation (PASS, FAIL, NOT_EVALUATED)
9. Taxon / metagenomic closed & open reference handling with expected_fraction=0 for false positives.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd


def compare_tables(
    df_a: pd.DataFrame | str | Path,
    df_b: pd.DataFrame | str | Path,
    id_col: str,
    value_cols: list[str],
    sample_cols: list[str] | None = None,
    rtol: float = 1e-5,
    atol: float = 1e-5,
    independent_mode: bool = False,
    file_a_path: str | Path | None = None,
    file_b_path: str | Path | None = None,
) -> dict[str, Any]:
    """
    Rigorously compare two tabular results (e.g. baseline vs replay, or comfy vs native).
    """
    # 1. Resolve paths and self-comparison check
    path_a = Path(file_a_path) if file_a_path else (Path(df_a) if isinstance(df_a, (str, Path)) else None)
    path_b = Path(file_b_path) if file_b_path else (Path(df_b) if isinstance(df_b, (str, Path)) else None)

    identical_file = False
    if path_a and path_b:
        try:
            if path_a.resolve() == path_b.resolve() or (path_a.exists() and path_b.exists() and os.path.samefile(path_a, path_b)):
                identical_file = True
        except Exception:
            pass

    if independent_mode and identical_file:
        return {
            "status": "FAIL",
            "reason": f"Self-comparison detected: input A and input B resolve to identical file ({path_a}). Independent replay must produce distinct outputs.",
            "identical_file_detected": True,
            "valid_numeric_comparisons": 0,
            "max_absolute_difference": None,
            "mean_absolute_difference": None,
        }

    # Load DataFrames if paths provided
    if isinstance(df_a, (str, Path)):
        df_a = pd.read_csv(df_a, sep=None, engine="python")
    if isinstance(df_b, (str, Path)):
        df_b = pd.read_csv(df_b, sep=None, engine="python")

    # Schema & column checks
    if id_col not in df_a.columns:
        return {"status": "FAIL", "reason": f"id_col '{id_col}' missing from table A", "valid_numeric_comparisons": 0}
    if id_col not in df_b.columns:
        return {"status": "FAIL", "reason": f"id_col '{id_col}' missing from table B", "valid_numeric_comparisons": 0}

    # Sample correspondence check
    sample_mismatches = []
    if sample_cols:
        for sc in sample_cols:
            if sc not in df_a.columns or sc not in df_b.columns:
                sample_mismatches.append(sc)
        if sample_mismatches:
            return {
                "status": "FAIL",
                "reason": f"Sample column mismatch: missing columns {sample_mismatches}",
                "sample_mismatches": sample_mismatches,
                "valid_numeric_comparisons": 0,
            }

    # 2. Check Feature ID uniqueness
    dups_a = df_a[id_col][df_a[id_col].duplicated()].tolist()
    dups_b = df_b[id_col][df_b[id_col].duplicated()].tolist()
    if dups_a or dups_b:
        return {
            "status": "FAIL",
            "reason": f"Duplicate feature IDs detected: A({len(dups_a)}), B({len(dups_b)})",
            "duplicate_ids_a": dups_a,
            "duplicate_ids_b": dups_b,
            "valid_numeric_comparisons": 0,
        }

    # 3. Total ID set check (Full Outer Join)
    set_a = set(df_a[id_col].dropna().astype(str))
    set_b = set(df_b[id_col].dropna().astype(str))
    common_ids = sorted(list(set_a & set_b))
    left_only = sorted(list(set_a - set_b))
    right_only = sorted(list(set_b - set_a))

    if left_only or right_only:
        return {
            "status": "FAIL",
            "reason": f"Feature ID set mismatch: {len(left_only)} in A only, {len(right_only)} in B only",
            "total_features_a": len(set_a),
            "total_features_b": len(set_b),
            "common_features": len(common_ids),
            "left_only_features": left_only,
            "right_only_features": right_only,
            "valid_numeric_comparisons": 0,
        }

    # Align by ID
    a_indexed = df_a.set_index(id_col).loc[common_ids]
    b_indexed = df_b.set_index(id_col).loc[common_ids]

    # 4. Numerical and NA comparisons
    all_diffs = []
    na_mismatches = []
    exceeding_count = 0
    total_valid_entries = 0

    for col in value_cols:
        if col not in a_indexed.columns or col not in b_indexed.columns:
            return {
                "status": "FAIL",
                "reason": f"Value column '{col}' missing from one of the tables",
                "valid_numeric_comparisons": 0,
            }

        va = a_indexed[col]
        vb = b_indexed[col]

        # Check NA mismatch
        na_a = va.isna()
        na_b = vb.isna()
        na_xor = na_a ^ na_b
        if na_xor.any():
            mismatched_ids = a_indexed.index[na_xor].tolist()
            na_mismatches.extend([(col, i) for i in mismatched_ids])

        # Valid numeric mask (both non-NA and finite)
        valid_mask = (~na_a) & (~na_b)
        va_valid = pd.to_numeric(va[valid_mask], errors="coerce")
        vb_valid = pd.to_numeric(vb[valid_mask], errors="coerce")
        both_num = va_valid.notna() & vb_valid.notna() & np.isfinite(va_valid) & np.isfinite(vb_valid)

        if both_num.any():
            diff = np.abs(va_valid[both_num].values - vb_valid[both_num].values)
            tol = atol + rtol * np.abs(vb_valid[both_num].values)
            exceeding = diff > tol
            exceeding_count += int(np.sum(exceeding))
            all_diffs.extend(diff.tolist())
            total_valid_entries += len(diff)

    if na_mismatches:
        return {
            "status": "FAIL",
            "reason": f"NA state mismatch across tables for {len(na_mismatches)} feature-column pairs",
            "na_mismatches": na_mismatches,
            "valid_numeric_comparisons": total_valid_entries,
        }

    # 5. Check if all values are NA or 0 valid comparisons
    if total_valid_entries == 0:
        return {
            "status": "NOT_EVALUATED",
            "reason": "All values are NA/Inf or 0 valid numerical comparisons available. PASS is strictly forbidden.",
            "total_features": len(common_ids),
            "valid_numeric_comparisons": 0,
            "max_absolute_difference": None,
            "mean_absolute_difference": None,
        }

    max_diff = float(np.max(all_diffs)) if all_diffs else 0.0
    mean_diff = float(np.mean(all_diffs)) if all_diffs else 0.0

    if exceeding_count > 0:
        return {
            "status": "FAIL",
            "reason": f"Numerical divergence exceeds tolerance (rtol={rtol}, atol={atol}) for {exceeding_count} entries. Max diff: {max_diff:.3e}",
            "total_features": len(common_ids),
            "valid_numeric_comparisons": total_valid_entries,
            "exceeding_tolerance_count": exceeding_count,
            "max_absolute_difference": max_diff,
            "mean_absolute_difference": mean_diff,
        }

    return {
        "status": "PASS",
        "reason": f"All {total_valid_entries} numerical values within tolerance (rtol={rtol}, atol={atol}). Max diff: {max_diff:.3e}",
        "total_features": len(common_ids),
        "valid_numeric_comparisons": total_valid_entries,
        "exceeding_tolerance_count": 0,
        "max_absolute_difference": max_diff,
        "mean_absolute_difference": mean_diff,
    }


def compare_metagenome_profiles(
    df_observed: pd.DataFrame,
    df_expected: pd.DataFrame,
    taxon_col: str = "name",
    abund_col: str = "fraction",
    expected_abund_col: str = "expected_fraction",
    detection_threshold: float = 0.0001,
) -> dict[str, Any]:
    """
    Rigorously compare metagenomic observed taxonomic abundance against expected composition.
    Outer join ensures unexpected taxa are NOT ignored: they are given expected_fraction=0
    and counted as False Positives and added to Total Variation Distance (TVD).
    """
    obs = df_observed[[taxon_col, abund_col]].copy()
    exp = df_expected[[taxon_col, expected_abund_col]].copy()

    merged = pd.merge(exp, obs, on=taxon_col, how="outer")
    merged[expected_abund_col] = merged[expected_abund_col].fillna(0.0)
    merged[abund_col] = merged[abund_col].fillna(0.0)

    # Re-normalize if sum > 0
    obs_sum = merged[abund_col].sum()
    if obs_sum > 0:
        merged["norm_observed"] = merged[abund_col] / obs_sum
    else:
        merged["norm_observed"] = 0.0

    exp_sum = merged[expected_abund_col].sum()
    if exp_sum > 0:
        merged["norm_expected"] = merged[expected_abund_col] / exp_sum
    else:
        merged["norm_expected"] = 0.0

    # Taxon detection metrics
    is_expected = merged["norm_expected"] >= detection_threshold
    is_detected = merged["norm_observed"] >= detection_threshold

    tp = int((is_expected & is_detected).sum())
    fp = int((~is_expected & is_detected).sum())
    fn = int((is_expected & ~is_detected).sum())

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Total Variation Distance: 0.5 * sum(|obs - exp|)
    merged["residual"] = merged["norm_observed"] - merged["norm_expected"]
    merged["abs_residual"] = merged["residual"].abs()
    tvd = float(0.5 * merged["abs_residual"].sum())

    unexpected_taxa = merged[~is_expected & is_detected][taxon_col].tolist()
    missing_taxa = merged[is_expected & ~is_detected][taxon_col].tolist()

    return {
        "status": "PASS" if fp == 0 and fn == 0 else "COMPOSITION_DISCREPANCY",
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "tvd": tvd,
        "unexpected_taxa": unexpected_taxa,
        "missing_taxa": missing_taxa,
        "merged_profile": merged,
    }
