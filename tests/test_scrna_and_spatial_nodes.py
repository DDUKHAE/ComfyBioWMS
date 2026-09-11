"""Tests for scRNA-seq and Transcriptomics/Spatial nodes.

Verifies real functional execution for:
- 7 scRNA nodes in nodes/ref_nodes.py
- In-memory analysis & CLI tools in nodes/transcriptomics_spatial_nodes.py
- Input file validation and FileNotFoundError raising
- Reproducible execution manifest generation
"""

import json
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import anndata as ad

from bioflow.runtime.command_runner import DryRunCommandRunner
from nodes.class_1 import (
    ScRNAReportNode,
    ScRNAVisualizationNode,
    ScanpyClusterNode,
    ScanpyMarkerGenesNode,
    ScanpyNormalizeNode,
    ScanpyQCNode,
)
from nodes.class_2 import TenxCountNode
from nodes.class_1 import (
    AnnDataIOInMemoryNode,
    Cell2locationAbundanceNode,
    CellPhoneDbInteractionNode,
    CellRankFateNode,
    GseapyEnrichmentNode,
    MuonMultimodalNode,
    PyScenicRegulonNode,
    ScVeloDynamicsNode,
    ScanviCellTypeNode,
    ScviToolsLatentNode,
    SquidpySpatialNode,
    TangramSpatialMappingNode,
)
from nodes.class_2 import (
    AlevinFryQuantNode,
    CiteSeqCountNode,
    EdgeRAnalysisNode,
    FlairIsoformNode,
    IsoToolsSplicingNode,
    KallistoQuantNode,
    LimmaVoomAnalysisNode,
    StarSoloQuantNode,
    StringTie2AssembleNode,
)


@pytest.fixture
def test_dir(tmp_path):
    d = tmp_path / "scrna_test"
    d.mkdir(parents=True, exist_ok=True)
    return d


def test_01_ref_nodes_scrna_pipeline_e2e(test_dir):
    """Test full e2e scRNA-seq pipeline across all 7 ref_nodes."""
    runner = DryRunCommandRunner()

    # Step 1: TenxCountNode validation & execution
    tenx = TenxCountNode()
    with pytest.raises(FileNotFoundError):
        tenx.run(str(test_dir / "nonexistent_fastqs"), "sample1", str(test_dir / "ref"), str(test_dir / "10x_out"), runner=runner)

    fastq_dir = test_dir / "fastqs"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    (fastq_dir / "sample_R1.fastq.gz").touch()
    (fastq_dir / "sample_R2.fastq.gz").touch()

    ref_dir = test_dir / "cellranger_ref"
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "genome.fa").touch()

    matrix_out_dir = test_dir / "cellranger_count" / "filtered_feature_bc_matrix"
    (res_matrix,) = tenx.run(str(fastq_dir), "sample1", str(ref_dir), str(matrix_out_dir), threads=4, runner=runner)
    assert Path(res_matrix).exists()
    assert (Path(res_matrix) / "matrix.mtx.gz").exists()
    assert (Path(res_matrix) / "barcodes.tsv.gz").exists()
    assert (Path(res_matrix) / "features.tsv.gz").exists()
    assert (matrix_out_dir.parent / "run_manifest.sh").exists()

    # Step 2: ScanpyQCNode
    qc_node = ScanpyQCNode()
    with pytest.raises(FileNotFoundError):
        qc_node.run("", str(test_dir / "nonexistent_matrix"), str(test_dir / "qc.h5ad"), runner=runner)

    qc_h5ad = test_dir / "scanpy" / "qc.h5ad"
    (res_qc,) = qc_node.run(res_matrix, res_matrix, str(qc_h5ad), min_genes=10, max_mito_pct=80.0, runner=runner)
    assert Path(res_qc).exists()
    adata_qc = ad.read_h5ad(res_qc)
    assert adata_qc.n_obs > 0
    assert "n_genes_by_counts" in adata_qc.obs
    assert "total_counts" in adata_qc.obs

    # Step 3: ScanpyNormalizeNode
    norm_node = ScanpyNormalizeNode()
    with pytest.raises(FileNotFoundError):
        norm_node.run("", str(test_dir / "nonexistent.h5ad"), str(test_dir / "norm.h5ad"), runner=runner)

    norm_h5ad = test_dir / "scanpy" / "normalized.h5ad"
    (res_norm,) = norm_node.run(res_qc, res_qc, str(norm_h5ad), target_sum=10000, runner=runner)
    assert Path(res_norm).exists()
    adata_norm = ad.read_h5ad(res_norm)
    assert adata_norm.n_obs == adata_qc.n_obs

    # Step 4: ScanpyClusterNode
    clust_node = ScanpyClusterNode()
    with pytest.raises(FileNotFoundError):
        clust_node.run("", str(test_dir / "nonexistent.h5ad"), str(test_dir / "clust.h5ad"), runner=runner)

    clust_h5ad = test_dir / "scanpy" / "clustered.h5ad"
    (res_clust,) = clust_node.run(res_norm, res_norm, str(clust_h5ad), n_pcs=10, resolution=0.8, runner=runner)
    assert Path(res_clust).exists()
    adata_clust = ad.read_h5ad(res_clust)
    assert "leiden" in adata_clust.obs
    assert "X_pca" in adata_clust.obsm
    assert "X_umap" in adata_clust.obsm

    # Step 5: ScanpyMarkerGenesNode
    marker_node = ScanpyMarkerGenesNode()
    with pytest.raises(FileNotFoundError):
        marker_node.run("", str(test_dir / "nonexistent.h5ad"), str(test_dir / "markers.csv"), runner=runner)

    markers_csv = test_dir / "scanpy" / "markers.csv"
    (res_markers,) = marker_node.run(res_clust, res_clust, str(markers_csv), groupby="leiden", method="wilcoxon", runner=runner)
    assert Path(res_markers).exists()
    m_df = pd.read_csv(res_markers)
    assert len(m_df) > 0

    # Step 6: ScRNAVisualizationNode
    viz_node = ScRNAVisualizationNode()
    with pytest.raises(FileNotFoundError):
        viz_node.run("", str(test_dir / "nonexistent.h5ad"), str(markers_csv), str(test_dir / "plots"), runner=runner)

    plot_dir = test_dir / "scanpy" / "plots"
    res_plot_dir, preview_tensor = viz_node.run(res_markers, res_clust, res_markers, str(plot_dir), runner=runner)
    assert Path(res_plot_dir).exists()
    assert (Path(res_plot_dir) / "umap.png").exists()
    assert (Path(res_plot_dir) / "pca.png").exists()
    assert preview_tensor is not None

    # Step 7: ScRNAReportNode
    rep_node = ScRNAReportNode()
    with pytest.raises(FileNotFoundError):
        rep_node.run(res_plot_dir, str(test_dir / "nonexistent.h5ad"), res_markers, res_plot_dir, str(test_dir / "report.md"), runner=runner)

    report_path = test_dir / "report" / "scrna_report.md"
    (res_report,) = rep_node.run(res_plot_dir, res_clust, res_markers, res_plot_dir, str(report_path), runner=runner)
    assert Path(res_report).exists()
    content = Path(res_report).read_text(encoding="utf-8")
    assert "# Single-Cell RNA-seq Analysis Report" in content
    assert "Total Filtered Cells" in content
    assert "Cluster Distribution" in content


def test_02_transcriptomics_spatial_in_memory_nodes(test_dir):
    """Test in-memory and python ML/DL nodes in transcriptomics_spatial_nodes.py."""
    # Create sample AnnData
    n_obs, n_vars = 40, 60
    X = np.random.poisson(2.0, (n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame({"cell_type": ["B_cell", "T_cell"] * (n_obs // 2), "leiden": ["0", "1"] * (n_obs // 2)}, index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame({"gene_symbol": [f"Gene_{j}" for j in range(n_vars)]}, index=[f"Gene_{j}" for j in range(n_vars)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["spliced"] = X
    adata.layers["unspliced"] = X * 0.3
    h5ad_file = test_dir / "toy.h5ad"
    adata.write_h5ad(h5ad_file)

    # 1. AnnDataIOInMemoryNode
    io_node = AnnDataIOInMemoryNode()
    with pytest.raises(FileNotFoundError):
        io_node.run(str(test_dir / "missing.h5ad"))

    summary_json, c_count, g_count = io_node.run(str(h5ad_file))
    assert c_count == n_obs
    assert g_count == n_vars
    summary_data = json.loads(summary_json)
    assert summary_data["cells"] == n_obs
    assert summary_data["genes"] == n_vars

    # 2. ScVeloDynamicsNode
    scv_node = ScVeloDynamicsNode()
    with pytest.raises(FileNotFoundError):
        scv_node.run(str(test_dir / "missing.h5ad"))

    velo_out = test_dir / "velocity.h5ad"
    res_velo, velo_summary = scv_node.run(str(h5ad_file), mode="stochastic", output_h5ad=str(velo_out))
    assert Path(res_velo).exists()
    assert "transition pairs" in velo_summary

    # 3. CellRankFateNode
    cr_node = CellRankFateNode()
    with pytest.raises(FileNotFoundError):
        cr_node.run(str(test_dir / "missing.h5ad"))

    fate_out = test_dir / "cellrank_fates.h5ad"
    res_fate, term_json = cr_node.run(str(res_velo), kernel="velocity", output_h5ad=str(fate_out))
    assert Path(res_fate).exists()
    terms = json.loads(term_json)
    assert "terminal_states" in terms

    # 4. SquidpySpatialNode
    sq_node = SquidpySpatialNode()
    with pytest.raises(FileNotFoundError):
        sq_node.run(str(test_dir / "missing.h5ad"))

    sq_out = test_dir / "squidpy_out.h5ad"
    res_sq, co_csv = sq_node.run(str(h5ad_file), cluster_key="cell_type", n_rings=2, output_h5ad=str(sq_out))
    assert Path(res_sq).exists()
    assert Path(co_csv).exists()

    # 5. TangramSpatialMappingNode
    tangram_node = TangramSpatialMappingNode()
    with pytest.raises(FileNotFoundError):
        tangram_node.run(str(test_dir / "missing_sc.h5ad"), str(h5ad_file), str(test_dir / "tangram.h5ad"))

    tangram_out = test_dir / "tangram_out.h5ad"
    res_tg, density_csv = tangram_node.run(str(h5ad_file), str(h5ad_file), str(tangram_out))
    assert Path(res_tg).exists()
    assert Path(density_csv).exists()

    # 6. MuonMultimodalNode
    muon_node = MuonMultimodalNode()
    with pytest.raises(FileNotFoundError):
        muon_node.run(str(test_dir / "missing_rna.h5ad"), str(h5ad_file), str(test_dir / "mu.h5mu"))

    mu_out = test_dir / "multimodal.h5mu"
    res_mu, mu_summary = muon_node.run(str(h5ad_file), str(h5ad_file), str(mu_out))
    assert Path(res_mu).exists()
    mu_dict = json.loads(mu_summary)
    assert "rna" in mu_dict["modalities"]

    # 7. GseapyEnrichmentNode
    gsea_node = GseapyEnrichmentNode()
    with pytest.raises(FileNotFoundError):
        gsea_node.run(str(test_dir / "missing_deg.csv"))

    deg_csv = test_dir / "deg_genes.csv"
    deg_csv.write_text("gene_symbol,logFC,pvalue\nTP53,2.5,0.001\nBRCA1,1.8,0.01\nCDK4,2.1,0.005\nMYC,3.0,0.0001\n", encoding="utf-8")
    enr_csv = test_dir / "enrichment.csv"
    res_enr, gsea_summary = gsea_node.run(str(deg_csv), gene_sets="MSigDB_Hallmark_2020", output_csv=str(enr_csv))
    assert Path(res_enr).exists()
    assert "significant terms" in gsea_summary

    # 8. PyScenicRegulonNode
    scenic_node = PyScenicRegulonNode()
    feather_file = test_dir / "motifs.feather"
    feather_file.touch()
    with pytest.raises(FileNotFoundError):
        scenic_node.run(str(test_dir / "missing.h5ad"), str(feather_file), str(test_dir / "scenic.h5ad"))

    scenic_out = test_dir / "scenic.h5ad"
    res_scenic, regulons_json = scenic_node.run(str(h5ad_file), str(feather_file), str(scenic_out))
    assert Path(res_scenic).exists()
    reg_data = json.loads(regulons_json)
    assert "regulons" in reg_data


def test_03_transcriptomics_cli_dispatch_and_manifests(test_dir):
    """Test CLI tools in transcriptomics_spatial_nodes.py dispatching commands with runner."""
    runner = DryRunCommandRunner()

    # Create dummy files for CLI testing
    r1 = test_dir / "R1.fq.gz"
    r2 = test_dir / "R2.fq.gz"
    bam = test_dir / "test.bam"
    gtf = test_dir / "test.gtf"
    fa = test_dir / "genome.fa"
    counts_csv = test_dir / "counts.csv"
    meta_csv = test_dir / "meta.csv"
    tags_csv = test_dir / "tags.csv"
    idx_file = test_dir / "transcripts.idx"
    rad_dir = test_dir / "rad_dir"
    star_idx = test_dir / "star_index"

    for f in [r1, r2, bam, gtf, fa, idx_file, tags_csv]:
        f.touch()
    counts_csv.write_text("gene,S1,S2\nGeneA,10,20\nGeneB,50,60\n", encoding="utf-8")
    meta_csv.write_text("sample,condition\nS1,control\nS2,treatment\n", encoding="utf-8")
    rad_dir.mkdir(parents=True, exist_ok=True)
    star_idx.mkdir(parents=True, exist_ok=True)

    # 1. KallistoQuantNode
    kallisto = KallistoQuantNode()
    with pytest.raises(FileNotFoundError):
        kallisto.run(str(test_dir / "missing.idx"), str(r1), str(test_dir / "kallisto_out"), runner=runner)
    abund_tsv, abund_h5 = kallisto.run(str(idx_file), str(r1), str(test_dir / "kallisto_out"), runner=runner)
    assert Path(abund_tsv).exists()
    assert (Path(test_dir / "kallisto_out") / "run_manifest.sh").exists()

    # 2. AlevinFryQuantNode
    alevin = AlevinFryQuantNode()
    with pytest.raises(FileNotFoundError):
        alevin.run(str(test_dir / "missing_rad"), str(test_dir / "alevin_out"), runner=runner)
    (alevin_out,) = alevin.run(str(rad_dir), str(test_dir / "alevin_out"), runner=runner)
    assert Path(alevin_out).exists()
    assert (Path(test_dir / "alevin_out") / "run_manifest.sh").exists()

    # 3. StarSoloQuantNode
    starsolo = StarSoloQuantNode()
    with pytest.raises(FileNotFoundError):
        starsolo.run(str(test_dir / "missing_star"), str(r1), str(r2), str(test_dir / "star_out"), runner=runner)
    solo_mat, solo_summary = starsolo.run(str(star_idx), str(r1), str(r2), str(test_dir / "star_out"), runner=runner)
    assert Path(solo_mat).exists()
    assert Path(solo_summary).exists()
    assert (Path(test_dir / "star_out") / "run_manifest.sh").exists()

    # 4. StringTie2AssembleNode
    stringtie = StringTie2AssembleNode()
    with pytest.raises(FileNotFoundError):
        stringtie.run(str(test_dir / "missing.bam"), str(gtf), str(test_dir / "stringtie.gtf"), runner=runner)
    out_gtf, abund_tsv = stringtie.run(str(bam), str(gtf), str(test_dir / "stringtie" / "transcripts.gtf"), runner=runner)
    assert Path(out_gtf).exists()
    assert Path(abund_tsv).exists()
    assert (Path(test_dir / "stringtie") / "run_manifest.sh").exists()

    # 5. FlairIsoformNode
    flair = FlairIsoformNode()
    with pytest.raises(FileNotFoundError):
        flair.run(str(test_dir / "missing.bam"), str(fa), str(gtf), str(test_dir / "flair" / "prefix"), runner=runner)
    flair_gtf, flair_counts = flair.run(str(bam), str(fa), str(gtf), str(test_dir / "flair" / "prefix"), runner=runner)
    assert Path(flair_gtf).exists()
    assert Path(flair_counts).exists()
    assert (Path(test_dir / "flair") / "run_manifest.sh").exists()

    # 6. IsoToolsSplicingNode
    isotools = IsoToolsSplicingNode()
    with pytest.raises(FileNotFoundError):
        isotools.run(str(test_dir / "missing.bam"), str(gtf), str(test_dir / "iso" / "out.tsv"), runner=runner)
    iso_tsv, iso_pkl = isotools.run(str(bam), str(gtf), str(test_dir / "iso" / "out.tsv"), runner=runner)
    assert Path(iso_tsv).exists()
    assert Path(iso_pkl).exists()
    assert (Path(test_dir / "iso") / "run_manifest.sh").exists()

    # 7. CiteSeqCountNode
    citeseq = CiteSeqCountNode()
    with pytest.raises(FileNotFoundError):
        citeseq.run(str(test_dir / "missing_r1"), str(r2), str(tags_csv), str(test_dir / "citeseq_out"), runner=runner)
    (cite_out,) = citeseq.run(str(r1), str(r2), str(tags_csv), str(test_dir / "citeseq_out"), runner=runner)
    assert Path(cite_out).exists()
    assert (Path(test_dir / "citeseq_out") / "run_manifest.sh").exists()

    # 8. EdgeRAnalysisNode
    edger = EdgeRAnalysisNode()
    with pytest.raises(FileNotFoundError):
        edger.run(str(test_dir / "missing_counts.csv"), str(meta_csv), str(test_dir / "edger" / "res.csv"), runner=runner)
    edger_res, edger_plot = edger.run(str(counts_csv), str(meta_csv), str(test_dir / "edger" / "res.csv"), runner=runner)
    assert Path(edger_res).exists()
    assert Path(edger_plot).exists()
    assert (Path(test_dir / "edger") / "run_manifest.sh").exists()

    # 9. LimmaVoomAnalysisNode
    limma = LimmaVoomAnalysisNode()
    with pytest.raises(FileNotFoundError):
        limma.run(str(test_dir / "missing_counts.csv"), str(meta_csv), str(test_dir / "limma" / "res.csv"), runner=runner)
    limma_res, limma_plot = limma.run(str(counts_csv), str(meta_csv), str(test_dir / "limma" / "res.csv"), runner=runner)
    assert Path(limma_res).exists()
    assert Path(limma_plot).exists()
    assert (Path(test_dir / "limma") / "run_manifest.sh").exists()
