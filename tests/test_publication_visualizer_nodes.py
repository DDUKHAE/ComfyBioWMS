import pytest
from pathlib import Path
import sys

root_dir = Path(__file__).resolve().parents[1]
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))
if str(root_dir / "engine" / "src") not in sys.path:
    sys.path.insert(0, str(root_dir / "engine" / "src"))

from nodes.publication_visualizer_nodes import (
    VolcanoPlotVisualizerNode,
    ManhattanPlotVisualizerNode,
    UmapScatterVisualizerNode,
    ClustermapHeatmapVisualizerNode,
    GseaEnrichmentPlotVisualizerNode,
    OncoPrintVisualizerNode,
    QqPlotVisualizerNode,
    MicrobiomeStackedBarVisualizerNode,
    PcoaScatterVisualizerNode,
    ChipAtacCoverageProfileVisualizerNode,
    MaPlotVisualizerNode,
    SpatialTissueOverlayVisualizerNode,
    KaplanMeierSurvivalVisualizerNode,
)

OUTPUT_BASE = root_dir / "results" / "test_visualizer_nodes"

@pytest.fixture(scope="module")
def setup_dirs():
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)


def test_01_volcano_and_ma_plots(setup_dirs):
    """Test Volcano Plot and MA Plot on real DEG outputs."""
    deg_csv = root_dir / "results" / "workflow_seqc_benchmark" / "deseq2_deg_results.csv"
    volcano_png = OUTPUT_BASE / "volcano_test.png"
    ma_png = OUTPUT_BASE / "ma_test.png"

    # Volcano
    volcano_node = VolcanoPlotVisualizerNode()
    v_path, v_tensor = volcano_node.run(
        deg_results_csv=str(deg_csv),
        log2fc_col="log2FoldChange",
        pvalue_col="pvalue",
        gene_col="gene_id",
        output_image_path=str(volcano_png),
        figure_title="SEQC DEG Volcano Plot",
        dpi=150,
    )
    assert Path(v_path).exists()
    print(f"\n[PASS] VolcanoPlotVisualizerNode rendered {v_path}.")

    # MA Plot
    ma_node = MaPlotVisualizerNode()
    ma_path, ma_tensor = ma_node.run(
        deg_results_csv=str(deg_csv),
        output_image_path=str(ma_png),
        figure_title="SEQC MA Plot",
        dpi=150,
    )
    assert Path(ma_path).exists()
    print(f"[PASS] MaPlotVisualizerNode rendered {ma_path}.")


def test_02_heatmap_and_gsea(setup_dirs):
    """Test Clustermap Heatmap and GSEA Enrichment curve."""
    heatmap_png = OUTPUT_BASE / "clustermap_test.png"
    gsea_png = OUTPUT_BASE / "gsea_test.png"

    # Heatmap
    heatmap_node = ClustermapHeatmapVisualizerNode()
    h_path, h_tensor = heatmap_node.run(
        matrix_csv="",
        output_image_path=str(heatmap_png),
        dpi=150,
    )
    assert Path(h_path).exists()
    print(f"\n[PASS] ClustermapHeatmapVisualizerNode rendered {h_path}.")

    # GSEA
    gsea_node = GseaEnrichmentPlotVisualizerNode()
    g_path, g_tensor = gsea_node.run(
        gsea_results_csv="",
        output_image_path=str(gsea_png),
        dpi=150,
    )
    assert Path(g_path).exists()
    print(f"[PASS] GseaEnrichmentPlotVisualizerNode rendered {g_path}.")


def test_03_gwas_and_survival_plots(setup_dirs):
    """Test Manhattan Plot, QQ Plot, and Kaplan-Meier survival curves."""
    manhattan_png = OUTPUT_BASE / "manhattan_test.png"
    qq_png = OUTPUT_BASE / "qq_test.png"
    km_png = OUTPUT_BASE / "km_test.png"

    # Manhattan
    manhattan_node = ManhattanPlotVisualizerNode()
    m_path, m_tensor = manhattan_node.run(
        gwas_summary_csv="",
        output_image_path=str(manhattan_png),
        dpi=150,
    )
    assert Path(m_path).exists()
    print(f"\n[PASS] ManhattanPlotVisualizerNode rendered {m_path}.")

    # QQ Plot
    qq_node = QqPlotVisualizerNode()
    q_path, q_tensor = qq_node.run(
        pvalues_csv="",
        output_image_path=str(qq_png),
        dpi=150,
    )
    assert Path(q_path).exists()
    print(f"[PASS] QqPlotVisualizerNode rendered {q_path}.")

    # Kaplan-Meier
    km_node = KaplanMeierSurvivalVisualizerNode()
    km_path, km_tensor = km_node.run(
        clinical_csv="",
        output_image_path=str(km_png),
        dpi=150,
    )
    assert Path(km_path).exists()
    print(f"[PASS] KaplanMeierSurvivalVisualizerNode rendered {km_path}.")


def test_04_single_cell_and_spatial_plots(setup_dirs):
    """Test UMAP, OncoPrint, and Spatial Tissue Overlay."""
    umap_png = OUTPUT_BASE / "umap_test.png"
    onco_png = OUTPUT_BASE / "oncoprint_test.png"
    spatial_png = OUTPUT_BASE / "spatial_test.png"

    # UMAP
    umap_node = UmapScatterVisualizerNode()
    u_path, u_tensor = umap_node.run(
        embedding_csv="",
        output_image_path=str(umap_png),
        dpi=150,
    )
    assert Path(u_path).exists()
    print(f"\n[PASS] UmapScatterVisualizerNode rendered {u_path}.")

    # OncoPrint
    onco_node = OncoPrintVisualizerNode()
    o_path, o_tensor = onco_node.run(
        mutation_maf_csv="",
        output_image_path=str(onco_png),
        dpi=150,
    )
    assert Path(o_path).exists()
    print(f"[PASS] OncoPrintVisualizerNode rendered {o_path}.")

    # Spatial Overlay
    spatial_node = SpatialTissueOverlayVisualizerNode()
    s_path, s_tensor = spatial_node.run(
        spatial_coords_csv="",
        output_image_path=str(spatial_png),
        dpi=150,
    )
    assert Path(s_path).exists()
    print(f"[PASS] SpatialTissueOverlayVisualizerNode rendered {s_path}.")


def test_05_microbiome_and_epigenomics_plots(setup_dirs):
    """Test Microbiome Stacked Bar, PCoA, and ATAC Coverage profiles."""
    microbiome_png = OUTPUT_BASE / "microbiome_test.png"
    pcoa_png = OUTPUT_BASE / "pcoa_test.png"
    atac_png = OUTPUT_BASE / "atac_profile_test.png"

    # Microbiome Stacked Bar
    mb_node = MicrobiomeStackedBarVisualizerNode()
    mb_path, mb_tensor = mb_node.run(
        abundance_csv="",
        output_image_path=str(microbiome_png),
        dpi=150,
    )
    assert Path(mb_path).exists()
    print(f"\n[PASS] MicrobiomeStackedBarVisualizerNode rendered {mb_path}.")

    # PCoA
    pcoa_node = PcoaScatterVisualizerNode()
    p_path, p_tensor = pcoa_node.run(
        distance_matrix_csv="",
        output_image_path=str(pcoa_png),
        dpi=150,
    )
    assert Path(p_path).exists()
    print(f"[PASS] PcoaScatterVisualizerNode rendered {p_path}.")

    # ATAC Coverage
    cov_node = ChipAtacCoverageProfileVisualizerNode()
    c_path, c_tensor = cov_node.run(
        coverage_matrix_gz="",
        output_image_path=str(atac_png),
        dpi=150,
    )
    assert Path(c_path).exists()
    print(f"[PASS] ChipAtacCoverageProfileVisualizerNode rendered {c_path}.")
