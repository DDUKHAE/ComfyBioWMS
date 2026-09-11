"""Regression test suite for execution contract, environment isolation, audit logging, and caching."""

import json
import os
import tempfile
import time
from pathlib import Path
import pytest

from bioflow.runtime.command_runner import (
    BioCommandRunner,
    resolve_tool_environment,
    validate_extra_command_conflicts,
    DEFAULT_TOOL_ENV_MAPPINGS,
)
from bioflow.runtime import artifacts
from bioflow.runtime.artifacts import (
    compute_input_fingerprint,
    get_run_output_dir,
    validate_artifact,
)
from nodes.class_2.fastp import Fastp
from nodes.class_2.salmon import SalmonIndex, SalmonQuantReads
from nodes.class_2.deseq2 import DESeq2
from nodes.class_1.tximport import Tximport
from nodes.class_1.plots import VolcanoPlot


def test_tool_environment_resolution():
    """Verify that canonical tools are resolved to valid conda environment executables."""
    env, path = resolve_tool_environment("fastp")
    assert env == "bulk_rna_seq"
    assert Path(path).is_file()
    assert os.access(path, os.X_OK)

    env, path = resolve_tool_environment("salmon")
    assert env == "bulk_rna_seq"
    assert Path(path).is_file()

    env, path = resolve_tool_environment("Rscript")
    assert env == "bulk_rna_seq"
    assert Path(path).is_file()


def test_missing_environment_rejection_prohibits_silent_fallback():
    """Verify that requesting a non-existent conda environment raises EnvironmentError."""
    with pytest.raises(EnvironmentError, match="Silent fallback to PATH is prohibited"):
        resolve_tool_environment("fastp", requested_env="nonexistent_virtual_env_xyz")


def test_extra_command_conflict_validation():
    """Verify that extra_command flags conflicting with UI parameters are strictly rejected."""
    managed = {"-o", "--out1", "-w", "--threads"}

    # Conflicting flag should raise ValueError
    with pytest.raises(ValueError, match="Conflict detected: option '-w'"):
        validate_extra_command_conflicts("-w 8", managed)

    with pytest.raises(ValueError, match="Conflict detected: option '-o'"):
        validate_extra_command_conflicts("-o custom_output.txt", managed)

    # Non-conflicting flags should pass cleanly
    tokens = validate_extra_command_conflicts("--verbose --dont_eval_duplication", managed)
    assert tokens == ["--verbose", "--dont_eval_duplication"]


def test_audit_manifest_logging_on_success_and_failure(tmp_path):
    """Verify that BioCommandRunner creates run_manifest.json, run_manifest.sh, and logs on both success and failure."""
    # 1. Success execution
    rec = BioCommandRunner.run(
        argv=["fastp", "--version"],
        cwd=tmp_path / "success_step",
        node_type="Fastp",
    )
    assert rec.returncode == 0
    manifest_json = tmp_path / "success_step" / "run_manifest.json"
    manifest_sh = tmp_path / "success_step" / "run_manifest.sh"
    assert manifest_json.is_file()
    assert manifest_sh.is_file()

    manifest_data = json.loads(manifest_json.read_text(encoding="utf-8"))
    assert manifest_data["status"] == "SUCCESS"
    assert manifest_data["environment"] == "bulk_rna_seq"
    assert manifest_data["returncode"] == 0
    assert manifest_data["duration_seconds"] >= 0.0

    # 2. Failure execution
    fail_dir = tmp_path / "fail_step"
    with pytest.raises(RuntimeError, match="failed with returncode"):
        BioCommandRunner.run(
            argv=["fastp", "-i", "nonexistent_reads.fastq.gz"],
            cwd=fail_dir,
            node_type="Fastp",
        )

    fail_manifest = fail_dir / "run_manifest.json"
    assert fail_manifest.is_file()
    fail_data = json.loads(fail_manifest.read_text(encoding="utf-8"))
    assert fail_data["status"] == "FAILED"
    assert fail_data["returncode"] != 0
    assert Path(fail_data["stderr_log"]).is_file()


def test_input_fingerprinting_caching_contract(tmp_path):
    """Verify that compute_input_fingerprint accurately reflects input and file state changes."""
    f1 = tmp_path / "reads_R1.fq.gz"
    f1.write_text("sequence_data_v1")

    fp1 = compute_input_fingerprint(read1=str(f1), threads=4)
    fp2 = compute_input_fingerprint(read1=str(f1), threads=4)
    assert fp1 == fp2, "Fingerprint must be deterministic for identical inputs"

    # Parameter change must alter fingerprint
    fp_param = compute_input_fingerprint(read1=str(f1), threads=8)
    assert fp1 != fp_param

    # File content / mtime change must alter fingerprint
    time.sleep(0.01)
    f1.write_text("sequence_data_v2")
    os.utime(f1, (time.time() + 5, time.time() + 5))
    fp_content = compute_input_fingerprint(read1=str(f1), threads=4)
    assert fp1 != fp_content


def test_directory_fingerprint_detects_changes_beyond_scan_ceiling(tmp_path, monkeypatch):
    """A directory input must never yield a stale-identical fingerprint.

    Files past the recursive scan ceiling are not hashed individually, so reaching the
    ceiling has to force re-execution rather than report an unchanged directory.
    """
    monkeypatch.setattr(artifacts, "MAX_DIR_SCAN_FILES", 5)

    quant_dir = tmp_path / "salmon_quant"
    quant_dir.mkdir()
    for i in range(4):
        (quant_dir / f"{i:03d}" ).mkdir()
        (quant_dir / f"{i:03d}" / "quant.sf").write_text("x")

    # Below the ceiling: deterministic, and sensitive to a late file's content.
    assert compute_input_fingerprint(quant_files=str(quant_dir)) == compute_input_fingerprint(
        quant_files=str(quant_dir)
    )
    before = compute_input_fingerprint(quant_files=str(quant_dir))
    last = quant_dir / "003" / "quant.sf"
    last.write_text("xxxx")
    os.utime(last, (time.time() + 5, time.time() + 5))
    assert compute_input_fingerprint(quant_files=str(quant_dir)) != before

    # At or above the ceiling: the scan is truncated, so the cache must be invalidated.
    for i in range(4, 12):
        (quant_dir / f"{i:03d}").mkdir()
        (quant_dir / f"{i:03d}" / "quant.sf").write_text("x")
    assert compute_input_fingerprint(quant_files=str(quant_dir)) != compute_input_fingerprint(
        quant_files=str(quant_dir)
    ), "Truncated directory scan must force re-execution"


def test_is_changed_implemented_on_core_nodes():
    """Verify that IS_CHANGED classmethod is implemented and callable on core analytical nodes."""
    assert hasattr(Fastp, "IS_CHANGED")
    assert hasattr(SalmonIndex, "IS_CHANGED")
    assert hasattr(SalmonQuantReads, "IS_CHANGED")
    assert hasattr(Tximport, "IS_CHANGED")
    assert hasattr(DESeq2, "IS_CHANGED")
    assert hasattr(VolcanoPlot, "IS_CHANGED")

    res = Fastp.IS_CHANGED(read1="data/test.fq", threads=4)
    assert isinstance(res, str)
    assert len(res) == 64  # SHA-256 hex string


def test_biological_empty_results_vs_corrupt_files(tmp_path):
    """Verify that biological 0-count results (e.g. 0 significant DEGs) are validated properly."""
    # 1. Valid empty DEG table (header present, 0 data rows)
    deg_empty = tmp_path / "deg_empty.csv"
    deg_empty.write_text("gene,baseMean,log2FoldChange,pvalue,padj\n")
    is_valid, reason = validate_artifact(deg_empty, "deg_table")
    assert is_valid, f"Expected valid empty DEG table: {reason}"

    # 2. Corrupt 0-byte file without header
    corrupt_file = tmp_path / "corrupt.csv"
    corrupt_file.touch()
    is_valid, reason = validate_artifact(corrupt_file, "deg_table", allow_empty=False)
    assert not is_valid, "0-byte file without header must be rejected"

    # 3. Missing file
    is_valid, reason = validate_artifact(tmp_path / "missing.csv", "deg_table")
    assert not is_valid
