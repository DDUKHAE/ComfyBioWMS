"""Tests for standalone MACS3 callpeak node.

Pure unit tests only: no external binaries executed, no fake subprocesses,
no monkeypatched execution, and no biological records fabricated.
"""

import ast
from pathlib import Path
import pytest

MODULE_PATH = Path("nodes/class_2/macs3.py").resolve()


def test_macs3_module_ast_no_project_or_sibling_imports():
    """Verify nodes/class_2/macs3.py has no project or sibling imports and correct docstring."""
    assert MODULE_PATH.is_file(), f"{MODULE_PATH} must exist"
    tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))

    docstring = ast.get_docstring(tree)
    assert docstring is not None, "Top module docstring is required"
    assert "Python packages: none" in docstring
    assert "External binary: macs3==3.0.4" in docstring
    assert "Galaxy wrapper:" in docstring
    assert "tools/macs2/macs2_callpeak.xml" in docstring

    allowed_stdlibs = {
        "ast", "collections", "inspect", "os", "sys", "shlex", "shutil",
        "subprocess", "threading", "pathlib", "typing", "re",
    }
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root_pkg = alias.name.split(".")[0]
                assert root_pkg in allowed_stdlibs, f"Disallowed import: {alias.name}"
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0, f"Relative imports are forbidden: level {node.level}"
            assert node.module is not None
            root_pkg = node.module.split(".")[0]
            assert root_pkg in allowed_stdlibs, f"Disallowed from-import: {node.module}"


@pytest.fixture(scope="module")
def macs3_mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("standalone_macs3", MODULE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_macs3_mappings_and_node_class(macs3_mod):
    """Verify NODE_CLASS_MAPPINGS and NODE_DISPLAY_NAME_MAPPINGS."""
    assert "Macs3Callpeak" in macs3_mod.NODE_CLASS_MAPPINGS
    assert "Macs3Callpeak" in macs3_mod.NODE_DISPLAY_NAME_MAPPINGS
    node_cls = macs3_mod.NODE_CLASS_MAPPINGS["Macs3Callpeak"]
    assert node_cls == macs3_mod.Macs3Callpeak
    assert hasattr(node_cls, "run")
    assert hasattr(node_cls, "INPUT_TYPES")


def test_macs3_input_types_and_return_types(macs3_mod):
    """Verify input types, defaults, and 7 stable return signatures."""
    inputs = macs3_mod.Macs3Callpeak.INPUT_TYPES()
    assert set(inputs["required"]) == {"treatment"}
    optional = inputs["optional"]

    expected_optionals = [
        "control", "sample_name", "format", "genome_size", "tag_size", "cutoff_mode",
        "cutoff_value", "keep_dup", "nomodel", "extsize", "shift",
        "mfold_lower", "mfold_upper", "bandwidth", "broad", "broad_cutoff",
        "call_summits", "bedgraph", "spmr", "nolambda", "slocal", "llocal",
        "scale_to", "extra_command",
    ]
    for opt in expected_optionals:
        assert opt in optional, f"Missing expected optional input: {opt}"

    # Verify all 10 documented format choices are available
    expected_formats = [
        "AUTO", "BAM", "SAM", "BED", "ELAND", "ELANDMULTI", "ELANDEXPORT", "BOWTIE", "BAMPE", "BEDPE"
    ]
    assert optional["format"][0] == expected_formats

    # Verify integer mfold inputs
    assert optional["mfold_lower"][0] == "INT"
    assert optional["mfold_upper"][0] == "INT"
    assert optional["mfold_lower"][1]["default"] == 5
    assert optional["mfold_upper"][1]["default"] == 50

    # Verify tag_size input
    assert optional["tag_size"][0] == "INT"
    assert optional["tag_size"][1]["default"] == 0

    # Verify scale_to choices
    assert optional["scale_to"][0] == ["small", "large"]
    assert optional["scale_to"][1]["default"] == "small"

    # Verify return types: exactly 7 string outputs
    assert macs3_mod.Macs3Callpeak.RETURN_TYPES == (
        "STRING", "STRING", "STRING", "STRING", "STRING", "STRING", "STRING"
    )
    assert macs3_mod.Macs3Callpeak.RETURN_NAMES == (
        "primary_peaks",
        "peaks_xls",
        "summits",
        "gapped_peaks",
        "treat_bedgraph",
        "control_bedgraph",
        "output_dir",
    )


def test_build_macs3_callpeak_argv_narrow_defaults(macs3_mod):
    """Verify argv construction for single-end narrow peak calling with defaults."""
    argv = macs3_mod.build_macs3_callpeak_argv(
        executable="macs3",
        treatment="/path/to/treat.bam",
        output_dir="/path/to/out",
        sample_name="my_sample",
        format="BAM",
        genome_size="hs",
        cutoff_mode="qvalue",
        cutoff_value=0.05,
        keep_dup="1",
        nomodel=False,
        mfold_lower=5,
        mfold_upper=50,
        bandwidth=300,
        broad=False,
        call_summits=True,
        scale_to="small",
    )
    assert argv[:2] == ["macs3", "callpeak"]
    assert ["-t", "/path/to/treat.bam"] == argv[argv.index("-t"):argv.index("-t") + 2]
    assert ["--outdir", "/path/to/out"] == argv[argv.index("--outdir"):argv.index("--outdir") + 2]
    assert ["-n", "my_sample"] == argv[argv.index("-n"):argv.index("-n") + 2]
    assert ["-f", "BAM"] == argv[argv.index("-f"):argv.index("-f") + 2]
    assert ["-g", "hs"] == argv[argv.index("-g"):argv.index("-g") + 2]
    assert ["-q", "0.05"] == argv[argv.index("-q"):argv.index("-q") + 2]
    assert ["--keep-dup", "1"] == argv[argv.index("--keep-dup"):argv.index("--keep-dup") + 2]
    assert ["--scale-to", "small"] == argv[argv.index("--scale-to"):argv.index("--scale-to") + 2]
    # Single-end model mode emits -m and --bw
    assert ["-m", "5", "50"] == argv[argv.index("-m"):argv.index("-m") + 3]
    assert ["--bw", "300"] == argv[argv.index("--bw"):argv.index("--bw") + 2]
    assert "--call-summits" in argv
    # Must omit incompatible single-end model flags
    assert "--nomodel" not in argv
    assert "--shift" not in argv
    assert "--extsize" not in argv
    # Must omit uninvented flags and control-only flags
    assert "-o" not in argv
    assert "--bandwidth" not in argv
    assert "--slocal" not in argv
    assert "--llocal" not in argv


def test_build_macs3_callpeak_argv_single_end_nomodel(macs3_mod):
    """Verify argv construction for single-end nomodel mode emits extsize/shift and omits mfold/bw."""
    argv = macs3_mod.build_macs3_callpeak_argv(
        executable="macs3",
        treatment="/path/to/treat.bam",
        output_dir="/path/to/out",
        sample_name="my_sample",
        format="BED",
        nomodel=True,
        extsize=150,
        shift=-75,
        mfold_lower=10,
        mfold_upper=40,
        bandwidth=250,
    )
    assert "--nomodel" in argv
    assert ["--extsize", "150"] == argv[argv.index("--extsize"):argv.index("--extsize") + 2]
    assert ["--shift", "-75"] == argv[argv.index("--shift"):argv.index("--shift") + 2]
    # In nomodel mode, mfold and bw must be omitted
    assert "-m" not in argv
    assert "--bw" not in argv


def test_build_macs3_callpeak_argv_paired_end_omits_incompatible_arguments(macs3_mod):
    """Verify BAMPE/BEDPE omits incompatible model, shift, extsize, mfold, and bandwidth arguments."""
    for pe_format in ["BAMPE", "BEDPE"]:
        argv = macs3_mod.build_macs3_callpeak_argv(
            executable="macs3",
            treatment="/path/to/pe.bam",
            output_dir="/path/to/out",
            sample_name="pe_sample",
            format=pe_format,
            nomodel=True,
            shift=50,
            extsize=200,
            mfold_lower=10,
            mfold_upper=40,
            bandwidth=300,
        )
        assert ["-f", pe_format] == argv[argv.index("-f"):argv.index("-f") + 2]
        assert "--nomodel" not in argv
        assert "--shift" not in argv
        assert "--extsize" not in argv
        assert "-m" not in argv
        assert "--bw" not in argv


def test_build_macs3_callpeak_argv_pvalue_and_uppercase_spmr(macs3_mod):
    """Tighten pvalue and uppercase SPMR assertions."""
    argv = macs3_mod.build_macs3_callpeak_argv(
        executable="macs3",
        treatment="/path/to/treat.bam",
        output_dir="/path/to/out",
        sample_name="spmr_sample",
        cutoff_mode="pvalue",
        cutoff_value=1e-5,
        tag_size=75,
        bedgraph=True,
        spmr=True,
    )
    # Tight pvalue assertion
    assert ["-p", "1e-05"] == argv[argv.index("-p"):argv.index("-p") + 2]
    assert "-q" not in argv
    # Tag size
    assert ["-s", "75"] == argv[argv.index("-s"):argv.index("-s") + 2]
    # BedGraph and uppercase SPMR
    assert "-B" in argv
    assert "--SPMR" in argv
    assert "--spmr" not in argv


def test_build_macs3_callpeak_argv_control_and_slocal_llocal(macs3_mod):
    """Verify slocal/llocal emitted only with control and without nolambda."""
    # With control and without nolambda -> emits slocal and llocal
    argv_ctrl = macs3_mod.build_macs3_callpeak_argv(
        control="/path/to/ctrl.bam",
        slocal=2000,
        llocal=15000,
        nolambda=False,
    )
    assert ["-c", "/path/to/ctrl.bam"] == argv_ctrl[argv_ctrl.index("-c"):argv_ctrl.index("-c") + 2]
    assert ["--slocal", "2000"] == argv_ctrl[argv_ctrl.index("--slocal"):argv_ctrl.index("--slocal") + 2]
    assert ["--llocal", "15000"] == argv_ctrl[argv_ctrl.index("--llocal"):argv_ctrl.index("--llocal") + 2]

    # With control and with nolambda -> omits slocal and llocal
    argv_nolambda = macs3_mod.build_macs3_callpeak_argv(
        control="/path/to/ctrl.bam",
        slocal=2000,
        llocal=15000,
        nolambda=True,
    )
    assert "--nolambda" in argv_nolambda
    assert "--slocal" not in argv_nolambda
    assert "--llocal" not in argv_nolambda

    # Without control -> omits slocal and llocal
    argv_no_ctrl = macs3_mod.build_macs3_callpeak_argv(
        control="",
        slocal=2000,
        llocal=15000,
    )
    assert "-c" not in argv_no_ctrl
    assert "--slocal" not in argv_no_ctrl
    assert "--llocal" not in argv_no_ctrl


def test_build_macs3_callpeak_argv_broad_mode(macs3_mod):
    """Verify broad mode emits --broad and --broad-cutoff and omits --call-summits."""
    argv = macs3_mod.build_macs3_callpeak_argv(
        broad=True,
        broad_cutoff=0.05,
        call_summits=True,
    )
    assert "--broad" in argv
    assert ["--broad-cutoff", "0.05"] == argv[argv.index("--broad-cutoff"):argv.index("--broad-cutoff") + 2]
    assert "--call-summits" not in argv


def test_validation_ordering(macs3_mod, tmp_path):
    """Verify validation order: invalid output_dir and sample_name fail before input files or PATH."""
    node = macs3_mod.Macs3Callpeak()

    # 1. Existing non-directory output_dir fails with ValueError before checking treatment
    file_as_outdir = tmp_path / "file_blocking_dir"
    file_as_outdir.touch()
    with pytest.raises(ValueError, match="not a directory"):
        node.run(treatment="/nonexistent/treat.bam", output_dir=str(file_as_outdir))

    # 3. Invalid sample name fails with ValueError before checking nonexistent treatment
    with pytest.raises(ValueError, match="Sample name"):
        node.run(treatment="/nonexistent/treat.bam", output_dir=str(tmp_path / "out"), sample_name="")
    with pytest.raises(ValueError, match="path separator"):
        node.run(treatment="/nonexistent/treat.bam", output_dir=str(tmp_path / "out"), sample_name="sub/sample")

    # 4. Missing treatment fails with FileNotFoundError before checking PATH or creating outdir
    missing_treat_out = tmp_path / "out_missing_treat"
    with pytest.raises(FileNotFoundError, match="Treatment"):
        node.run(treatment=str(tmp_path / "absent.bam"), output_dir=str(missing_treat_out), sample_name="s1")
    assert not missing_treat_out.exists(), "Output directory must not be created on validation failure"

    # 5. Missing control fails with FileNotFoundError before checking PATH or creating outdir
    real_empty_treat = tmp_path / "empty_treat.bam"
    real_empty_treat.touch()
    missing_ctrl_out = tmp_path / "out_missing_ctrl"
    with pytest.raises(FileNotFoundError, match="Control"):
        node.run(
            treatment=str(real_empty_treat),
            output_dir=str(missing_ctrl_out),
            control=str(tmp_path / "absent_ctrl.bam"),
            sample_name="s1",
        )
    assert not missing_ctrl_out.exists(), "Output directory must not be created on validation failure"


def test_missing_binary_fails_before_mkdir(macs3_mod, tmp_path, monkeypatch):
    """Verify missing binary fails after input validation but before creating output directory."""
    empty_treat = tmp_path / "valid_file.bam"
    empty_treat.touch()
    out_dir = tmp_path / "no_binary_dir"
    monkeypatch.setenv("PATH", "")

    node = macs3_mod.Macs3Callpeak()
    with pytest.raises(RuntimeError, match="macs3 executable not found on PATH"):
        node.run(treatment=str(empty_treat), output_dir=str(out_dir), sample_name="s1")
    assert not out_dir.exists(), "Output directory must not be created if binary is missing"


def test_bed_peak_validation_rejects_empty_file_and_no_intervals(macs3_mod, tmp_path):
    """Verify _validate_bed_peaks rejects empty files and files with no interval rows."""
    empty_file = tmp_path / "empty.narrowPeak"
    empty_file.touch()
    with pytest.raises(RuntimeError, match="empty"):
        macs3_mod._validate_bed_peaks(empty_file)

    comments_only = tmp_path / "comments_only.narrowPeak"
    comments_only.write_text("# track description\n# second comment\n")
    with pytest.raises(RuntimeError, match="no interval rows"):
        macs3_mod._validate_bed_peaks(comments_only)


def test_output_artifacts_resolution_narrow_vs_broad(macs3_mod, tmp_path):
    """Verify artifact paths resolution for narrow vs broad modes including 7 outputs."""
    out_dir = tmp_path / "artifacts_out"

    # Narrow mode without bedgraph
    narrow_paths = macs3_mod.resolve_output_paths(
        out_dir=out_dir,
        sample_name="narrow_sample",
        broad=False,
        bedgraph=False,
    )
    assert narrow_paths["primary_peaks"] == out_dir / "narrow_sample_peaks.narrowPeak"
    assert narrow_paths["peaks_xls"] == out_dir / "narrow_sample_peaks.xls"
    assert narrow_paths["summits"] == out_dir / "narrow_sample_summits.bed"
    assert narrow_paths["gapped_peaks"] is None
    assert narrow_paths["treat_bedgraph"] is None
    assert narrow_paths["control_bedgraph"] is None
    assert narrow_paths["output_dir"] == out_dir

    # Broad mode with bedgraph
    broad_paths = macs3_mod.resolve_output_paths(
        out_dir=out_dir,
        sample_name="broad_sample",
        broad=True,
        bedgraph=True,
    )
    assert broad_paths["primary_peaks"] == out_dir / "broad_sample_peaks.broadPeak"
    assert broad_paths["peaks_xls"] == out_dir / "broad_sample_peaks.xls"
    assert broad_paths["summits"] is None
    assert broad_paths["gapped_peaks"] == out_dir / "broad_sample_peaks.gappedPeak"
    assert broad_paths["treat_bedgraph"] == out_dir / "broad_sample_treat_pileup.bdg"
    assert broad_paths["control_bedgraph"] == out_dir / "broad_sample_control_lambda.bdg"
    assert broad_paths["output_dir"] == out_dir
