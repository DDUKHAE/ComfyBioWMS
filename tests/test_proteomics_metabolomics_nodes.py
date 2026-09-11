"""Unit and functional tests for Proteomics & Metabolomics nodes."""

import json
import numpy as np
import pytest
from pathlib import Path

from bioflow.runtime.command_runner import DryRunCommandRunner
from nodes.class_1 import (
    MassqlQueryNode,
    MatchmsSpectrumNode,
    MsDeisotopeNode,
    PyOpenMSFeatureNode,
    PyteomicsMSNode,
    Spec2VecEmbeddingNode,
)
from nodes.class_2 import (
    DiaNNQuantNode,
    MaxQuantQuantNode,
    MetaboAnalystRNode,
    MhcquantNeoantigenNode,
    MsDialLipidNode,
    MsconvertConvertNode,
    MsfraggerSearchNode,
    PerseusCliNode,
    SiriusStructureNode,
)


@pytest.fixture
def fixtures_dir(tmp_path):
    """Generate realistic test fixtures for mzML, MGF, FASTA, TSV, and CSV files."""
    d = tmp_path / "fixtures"
    d.mkdir(parents=True, exist_ok=True)

    # 1. FASTA file
    fasta_path = d / "test_uniprot.fasta"
    fasta_content = (
        ">P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2\n"
        "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDK\n"
        ">P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2\n"
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKL\n"
    )
    fasta_path.write_text(fasta_content, encoding="utf-8")

    # 2. MGF file (Query)
    query_mgf = d / "test_query.mgf"
    query_mgf_content = (
        "BEGIN IONS\n"
        "TITLE=Spectrum_1\n"
        "PEPMASS=445.12 100000\n"
        "CHARGE=2+\n"
        "RTINSECONDS=65.2\n"
        "SCANS=1\n"
        "120.08 500.0\n"
        "226.18 1200.0\n"
        "305.15 800.0\n"
        "445.12 3500.0\n"
        "END IONS\n"
        "BEGIN IONS\n"
        "TITLE=Spectrum_2\n"
        "PEPMASS=580.32 85000\n"
        "CHARGE=2+\n"
        "RTINSECONDS=120.5\n"
        "SCANS=2\n"
        "150.05 300.0\n"
        "280.12 950.0\n"
        "580.32 4200.0\n"
        "END IONS\n"
    )
    query_mgf.write_text(query_mgf_content, encoding="utf-8")

    # 3. MGF file (Library)
    lib_mgf = d / "test_library.mgf"
    lib_mgf_content = (
        "BEGIN IONS\n"
        "TITLE=Caffeine_Standard\n"
        "PEPMASS=195.088\n"
        "CHARGE=1+\n"
        "RTINSECONDS=62.0\n"
        "SCANS=101\n"
        "110.05 200.0\n"
        "138.06 600.0\n"
        "195.09 5000.0\n"
        "END IONS\n"
        "BEGIN IONS\n"
        "TITLE=Serum_Albumin_Fragment\n"
        "PEPMASS=445.12\n"
        "CHARGE=2+\n"
        "RTINSECONDS=65.0\n"
        "SCANS=102\n"
        "120.08 450.0\n"
        "226.18 1150.0\n"
        "305.15 750.0\n"
        "445.12 3400.0\n"
        "END IONS\n"
    )
    lib_mgf.write_text(lib_mgf_content, encoding="utf-8")

    # 4. mzML file
    mzml_path = d / "test_spectra.mzML"
    mzml_content = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="run_test">
    <spectrumList count="2">
      <spectrum id="scan=1" index="0" defaultArrayLength="3">
        <cvParam accession="MS:1000511" name="ms level" value="1"/>
        <cvParam accession="MS:1000016" name="scan start time" value="1.0" unitName="minute"/>
      </spectrum>
      <spectrum id="scan=2" index="1" defaultArrayLength="4">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <cvParam accession="MS:1000016" name="scan start time" value="1.2" unitName="minute"/>
        <cvParam accession="MS:1000744" name="selected ion m/z" value="445.12"/>
      </spectrum>
    </spectrumList>
  </run>
</mzML>
"""
    mzml_path.write_text(mzml_content, encoding="utf-8")

    # 5. Protein Groups TSV (MaxQuant / Perseus)
    protein_tsv = d / "proteinGroups.txt"
    tsv_content = (
        "Protein IDs\tMajority protein IDs\tIntensity Sample1\tIntensity Sample2\tIntensity Sample3\n"
        "P02768\tP02768\t150000\t180000\t165000\n"
        "P68871\tP68871\t85000\t92000\t0\n"
        "P04406\tP04406\t320000\t310000\t335000\n"
    )
    protein_tsv.write_text(tsv_content, encoding="utf-8")

    # 6. Peak table CSV & Metadata CSV (MetaboAnalyst)
    peaks_csv = d / "peaks.csv"
    peaks_content = (
        "Compound,Sample_A1,Sample_A2,Sample_B1,Sample_B2\n"
        "Glucose,15400,16200,9800,10200\n"
        "Lactate,8200,8900,14200,13800\n"
        "Citrate,5100,5300,4800,5000\n"
    )
    peaks_csv.write_text(peaks_content, encoding="utf-8")

    meta_csv = d / "meta.csv"
    meta_content = (
        "Sample,Group\n"
        "Sample_A1,Control\n"
        "Sample_A2,Control\n"
        "Sample_B1,Treated\n"
        "Sample_B2,Treated\n"
    )
    meta_csv.write_text(meta_content, encoding="utf-8")

    # 7. Dummy vendor raw and speclib files
    (d / "sample.raw").write_bytes(b"\x00\x01\x02RAW_DATA")
    (d / "library.speclib").write_text("SPECLIB_DUMMY_DATA", encoding="utf-8")
    (d / "LipidBlast.msp").write_text("NAME: PC(16:0/18:1)\nPRECURSORMZ: 760.585\n", encoding="utf-8")

    return d


# =============================================================================
# Tier 1 In-Memory & Python MS Processing Tests
# =============================================================================

def test_pyteomics_ms_node(fixtures_dir):
    """Test PyteomicsMSNode parses spectra and fasta and returns real PSMs."""
    node = PyteomicsMSNode()

    # Test missing file
    with pytest.raises(FileNotFoundError):
        node.run("nonexistent.mzML", str(fixtures_dir / "test_uniprot.fasta"))
    with pytest.raises(FileNotFoundError):
        node.run(str(fixtures_dir / "test_spectra.mzML"), "nonexistent.fasta")

    # Test real execution
    summary_json, total_psms = node.run(
        str(fixtures_dir / "test_spectra.mzML"),
        str(fixtures_dir / "test_uniprot.fasta"),
    )
    assert isinstance(summary_json, str)
    assert isinstance(total_psms, int)
    data = json.loads(summary_json)
    assert data["total_spectra"] == 2
    assert data["ms1_spectra"] == 1
    assert data["ms2_spectra"] == 1
    assert "total_peaks" in data
    print(f"\n[PASS] PyteomicsMSNode parsed {data['total_spectra']} spectra, identified {total_psms} PSMs.")


def test_pyopenms_feature_node(fixtures_dir, tmp_path):
    """Test PyOpenMSFeatureNode detects LC-MS features and writes CSV."""
    node = PyOpenMSFeatureNode()
    out_csv = tmp_path / "features.csv"

    # Test missing file
    with pytest.raises(FileNotFoundError):
        node.run("nonexistent.mzML", str(out_csv))

    csv_path, feature_count = node.run(
        str(fixtures_dir / "test_spectra.mzML"),
        str(out_csv),
    )
    assert Path(csv_path).exists()
    assert feature_count >= 0
    lines = Path(csv_path).read_text(encoding="utf-8").splitlines()
    assert len(lines) >= 1
    assert "feature_id,rt,mz,intensity,charge,quality" in lines[0]
    print(f"\n[PASS] PyOpenMSFeatureNode detected {feature_count} features into {csv_path}.")


def test_matchms_spectrum_node(fixtures_dir, tmp_path):
    """Test MatchmsSpectrumNode computes spectral cosine similarity and saves matrix."""
    node = MatchmsSpectrumNode()
    out_npy = tmp_path / "sim_matrix.npy"

    # Test missing file
    with pytest.raises(FileNotFoundError):
        node.run("nonexistent_q.mgf", str(fixtures_dir / "test_library.mgf"), output_npy=str(out_npy))

    npy_path, matches_json = node.run(
        str(fixtures_dir / "test_query.mgf"),
        str(fixtures_dir / "test_library.mgf"),
        output_npy=str(out_npy),
    )
    assert Path(npy_path).exists()
    matrix = np.load(npy_path)
    assert matrix.shape == (2, 2)
    # The second query (445.12) should have high cosine similarity with the library albumin fragment (445.12)
    assert matrix[0, 1] > 0.9  # high cosine similarity for matching peaks
    data = json.loads(matches_json)
    assert data["query_spectra_count"] == 2
    assert data["library_spectra_count"] == 2
    print(f"\n[PASS] MatchmsSpectrumNode computed similarity matrix shape {matrix.shape}, top score: {matrix[0, 1]:.4f}.")


def test_spec2vec_embedding_node(fixtures_dir, tmp_path):
    """Test Spec2VecEmbeddingNode computes embeddings and returns vector dimension."""
    node = Spec2VecEmbeddingNode()
    out_npy = tmp_path / "spec2vec_embeddings.npy"

    # Test missing MGF
    with pytest.raises(FileNotFoundError):
        node.run("nonexistent.mgf", "", output_npy=str(out_npy))

    # Test missing explicit model file
    with pytest.raises(FileNotFoundError):
        node.run(str(fixtures_dir / "test_query.mgf"), "missing_model.model", output_npy=str(out_npy))

    npy_path, dim = node.run(
        str(fixtures_dir / "test_query.mgf"),
        "models/spec2vec_model.model",
        output_npy=str(out_npy),
    )
    assert Path(npy_path).exists()
    embeddings = np.load(npy_path)
    assert embeddings.shape[0] == 2
    assert embeddings.shape[1] == dim
    print(f"\n[PASS] Spec2VecEmbeddingNode generated embeddings shape {embeddings.shape}.")


def test_massql_query_node(fixtures_dir, tmp_path):
    """Test MassqlQueryNode filters scans matching MassQL query."""
    node = MassqlQueryNode()
    out_csv = tmp_path / "massql_out.csv"

    with pytest.raises(FileNotFoundError):
        node.run("nonexistent.mzML", "QUERY scaninfo(MS2DATA)", str(out_csv))

    # Query with MGF query file
    csv_path, matched_scans = node.run(
        str(fixtures_dir / "test_query.mgf"),
        "QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.18:TOLERANCEPPM=20",
        str(out_csv),
    )
    assert Path(csv_path).exists()
    assert matched_scans == 1  # Spectrum_1 contains peak 226.18
    print(f"\n[PASS] MassqlQueryNode found {matched_scans} matching scans.")


def test_ms_deisotope_node(fixtures_dir, tmp_path):
    """Test MsDeisotopeNode deisotopes spectra and writes valid XML mzML."""
    node = MsDeisotopeNode()
    out_mzml = tmp_path / "deisotoped.mzML"

    with pytest.raises(FileNotFoundError):
        node.run("nonexistent.mzML", str(out_mzml))

    result_mzml, peak_count = node.run(
        str(fixtures_dir / "test_query.mgf"),
        str(out_mzml),
    )
    assert Path(result_mzml).exists()
    content = Path(result_mzml).read_text(encoding="utf-8")
    assert "<mzML" in content and "</mzML>" in content
    assert peak_count > 0
    print(f"\n[PASS] MsDeisotopeNode deconvoluted {peak_count} peaks into {result_mzml}.")


# =============================================================================
# Tier 2 & Tier 3 CLI Dispatches & Manifest Tests
# =============================================================================

def test_diann_quant_node(fixtures_dir, tmp_path):
    """Test DiaNNQuantNode dispatches CLI and creates execution manifest."""
    node = DiaNNQuantNode()
    out_dir = tmp_path / "diann_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.mzML", str(fixtures_dir / "library.speclib"), str(fixtures_dir / "test_uniprot.fasta"), str(out_dir), runner=runner)

    report, matrix = node.run(
        str(fixtures_dir / "test_spectra.mzML"),
        str(fixtures_dir / "library.speclib"),
        str(fixtures_dir / "test_uniprot.fasta"),
        str(out_dir),
        runner=runner,
    )
    assert Path(report).name == "report.tsv"
    assert Path(matrix).name == "report.pg_matrix.tsv"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "diann" for cmd in runner.commands)
    print("\n[PASS] DiaNNQuantNode dispatched diann CLI and generated manifest.")


def test_msfragger_search_node(fixtures_dir, tmp_path):
    """Test MsfraggerSearchNode dispatches CLI and creates execution manifest."""
    node = MsfraggerSearchNode()
    out_dir = tmp_path / "msfragger_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.mzML", str(fixtures_dir / "test_uniprot.fasta"), str(out_dir), runner=runner)

    psm, protein = node.run(
        str(fixtures_dir / "test_spectra.mzML"),
        str(fixtures_dir / "test_uniprot.fasta"),
        str(out_dir),
        runner=runner,
    )
    assert Path(psm).name == "psm.tsv"
    assert Path(protein).name == "protein.tsv"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "msfragger" for cmd in runner.commands)
    print("\n[PASS] MsfraggerSearchNode dispatched msfragger CLI and generated manifest.")


def test_msconvert_convert_node(fixtures_dir, tmp_path):
    """Test MsconvertConvertNode dispatches CLI and creates execution manifest."""
    node = MsconvertConvertNode()
    out_dir = tmp_path / "msconvert_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.raw", str(out_dir), runner=runner)

    (converted,) = node.run(
        str(fixtures_dir / "sample.raw"),
        str(out_dir),
        runner=runner,
    )
    assert Path(converted).name == "sample.mzML"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "msconvert" for cmd in runner.commands)
    print("\n[PASS] MsconvertConvertNode dispatched msconvert CLI and generated manifest.")


def test_msdial_lipid_node(fixtures_dir, tmp_path):
    """Test MsDialLipidNode dispatches CLI and creates execution manifest."""
    node = MsDialLipidNode()
    out_dir = tmp_path / "msdial_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing_dir", str(fixtures_dir / "LipidBlast.msp"), str(out_dir), runner=runner)

    (peak_table,) = node.run(
        str(fixtures_dir),
        str(fixtures_dir / "LipidBlast.msp"),
        str(out_dir),
        runner=runner,
    )
    assert Path(peak_table).name == "Area_0_PeakID_All.txt"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "msdial" for cmd in runner.commands)
    print("\n[PASS] MsDialLipidNode dispatched msdial CLI and generated manifest.")


def test_sirius_structure_node(fixtures_dir, tmp_path):
    """Test SiriusStructureNode dispatches CLI and creates execution manifest."""
    node = SiriusStructureNode()
    out_dir = tmp_path / "sirius_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.mgf", str(out_dir), runner=runner)

    formula, fp = node.run(
        str(fixtures_dir / "test_query.mgf"),
        str(out_dir),
        runner=runner,
    )
    assert Path(formula).name == "formula_identifications.tsv"
    assert Path(fp).name == "csi_fingerprints.tsv"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "sirius" for cmd in runner.commands)
    print("\n[PASS] SiriusStructureNode dispatched sirius CLI and generated manifest.")


def test_maxquant_quant_node(fixtures_dir, tmp_path):
    """Test MaxQuantQuantNode dispatches CLI and creates execution manifest."""
    node = MaxQuantQuantNode()
    out_dir = tmp_path / "maxquant_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.raw", str(fixtures_dir / "test_uniprot.fasta"), str(out_dir), runner=runner)

    pg, evidence = node.run(
        str(fixtures_dir / "sample.raw"),
        str(fixtures_dir / "test_uniprot.fasta"),
        str(out_dir),
        runner=runner,
    )
    assert Path(pg).name == "proteinGroups.txt"
    assert Path(evidence).name == "evidence.txt"
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "maxquant" for cmd in runner.commands)
    print("\n[PASS] MaxQuantQuantNode dispatched maxquant CLI and generated manifest.")


def test_perseus_cli_node(fixtures_dir, tmp_path):
    """Test PerseusCliNode executes real statistical ANOVA and imputation."""
    node = PerseusCliNode()
    out_tsv = tmp_path / "perseus_out" / "stat_results.tsv"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.txt", str(out_tsv), runner=runner)

    imputed, anova = node.run(
        str(fixtures_dir / "proteinGroups.txt"),
        str(out_tsv),
        runner=runner,
    )
    assert Path(imputed).exists()
    assert Path(anova).exists()
    anova_lines = Path(anova).read_text(encoding="utf-8").splitlines()
    assert len(anova_lines) >= 2  # header + proteins
    assert "Protein\tF_statistic\tp_value\tsignificant" in anova_lines[0]
    assert (out_tsv.parent / "run_manifest.sh").exists()
    assert (out_tsv.parent / "run_manifest.json").exists()
    print(f"\n[PASS] PerseusCliNode generated imputed matrix and ANOVA results ({len(anova_lines)-1} proteins).")


def test_metaboanalyst_r_node(fixtures_dir, tmp_path):
    """Test MetaboAnalystRNode generates reproducible R script and computes PLS-DA."""
    node = MetaboAnalystRNode()
    out_dir = tmp_path / "metaboanalyst_out"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.csv", str(fixtures_dir / "meta.csv"), str(out_dir), runner=runner)

    enrichment, plsda = node.run(
        str(fixtures_dir / "peaks.csv"),
        str(fixtures_dir / "meta.csv"),
        str(out_dir),
        runner=runner,
    )
    assert Path(enrichment).exists()
    assert Path(plsda).exists()
    assert (out_dir / "run_metaboanalyst.R").exists()
    assert (out_dir / "run_manifest.sh").exists()
    assert (out_dir / "run_manifest.json").exists()
    plsda_lines = Path(plsda).read_text(encoding="utf-8").splitlines()
    assert len(plsda_lines) >= 4  # header + 4 samples
    assert "Sample,PLS_Component1,PLS_Component2" in plsda_lines[0]
    print(f"\n[PASS] MetaboAnalystRNode computed PLS-DA scores for {len(plsda_lines)-1} samples.")


def test_mhcquant_neoantigen_node(fixtures_dir, tmp_path):
    """Test MhcquantNeoantigenNode dispatches CLI and creates execution manifest."""
    node = MhcquantNeoantigenNode()
    out_tsv = tmp_path / "mhcquant" / "neoantigens.tsv"
    runner = DryRunCommandRunner()

    with pytest.raises(FileNotFoundError):
        node.run("missing.mzML", str(fixtures_dir / "test_uniprot.fasta"), "HLA-A*02:01", str(out_tsv), runner=runner)

    (res_tsv,) = node.run(
        str(fixtures_dir / "test_spectra.mzML"),
        str(fixtures_dir / "test_uniprot.fasta"),
        "HLA-A*02:01",
        str(out_tsv),
        runner=runner,
    )
    assert Path(res_tsv).name == "neoantigens.tsv"
    assert (out_tsv.parent / "run_manifest.sh").exists()
    assert (out_tsv.parent / "run_manifest.json").exists()
    assert any(cmd.argv[0] == "mhcquant" for cmd in runner.commands)
    print("\n[PASS] MhcquantNeoantigenNode dispatched mhcquant CLI and generated manifest.")
