import importlib.util
import json
from pathlib import Path
import pandas as pd
import pytest

DESEQ2_MODULE = Path("nodes/class_2/deseq2.py").resolve()
PLOTS_MODULE = Path("nodes/class_1/plots.py").resolve()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_deseq2_module_loads_standalone_and_owns_mappings():
    mod = load_module(DESEQ2_MODULE, "standalone_deseq2")
    assert "DESeq2" in mod.NODE_CLASS_MAPPINGS
    assert "DESeq2" in mod.NODE_DISPLAY_NAME_MAPPINGS

    cls = mod.NODE_CLASS_MAPPINGS["DESeq2"]
    assert hasattr(cls, "INPUT_TYPES")
    assert hasattr(cls, "RETURN_TYPES")
    assert cls.CATEGORY == "ComfyBIO/RNA-Seq"
    assert cls.FUNCTION == "run"

    inputs = cls.INPUT_TYPES()
    assert "count_matrix_csv" in inputs["required"]
    assert "sample_metadata_csv" in inputs["required"]
    assert "extra_command" in inputs["optional"]

    # Backward compatibility aliases
    assert mod.DESeq2Node is cls
    assert mod.DESeq2Analysis is cls
    assert mod.DESeq2AnalysisNode is cls


def test_deseq2_differential_expression_synthetic_dataset(tmp_path):
    mod = load_module(DESEQ2_MODULE, "standalone_deseq2")
    node = mod.DESeq2()

    # Generate synthetic count matrix: 3 control samples, 3 treated samples
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]
    genes = [f"Gene_{i:02d}" for i in range(1, 21)]

    # Gene_01: strongly upregulated in trt
    # Gene_02: strongly downregulated in trt
    # Gene_03..20: house-keeping / invariant background
    data = []
    for g in genes:
        if g == "Gene_01":
            row = [g, 15, 18, 16, 520, 560, 540]
        elif g == "Gene_02":
            row = [g, 620, 590, 610, 12, 15, 14]
        else:
            row = [g, 200, 210, 195, 205, 198, 202]
        data.append(row)

    count_df = pd.DataFrame(data, columns=["gene_id"] + samples)
    counts_file = tmp_path / "raw_counts.csv"
    count_df.to_csv(counts_file, index=False)

    meta_df = pd.DataFrame({
        "sample_id": samples,
        "condition": ["control", "control", "control", "treated", "treated", "treated"],
    })
    meta_file = tmp_path / "sample_metadata.csv"
    meta_df.to_csv(meta_file, index=False)

    out_dir = tmp_path / "deseq2_out"
    res_csv, norm_csv, summary_json = node.run(
        count_matrix_csv=str(counts_file),
        sample_metadata_csv=str(meta_file),
        output_dir=str(out_dir),
        condition_col="condition",
        contrast_reference="control",
        contrast_target="treated",
        alpha=0.05,
        lfc_threshold=1.0,
    )

    assert Path(res_csv).is_file()
    assert Path(norm_csv).is_file()

    # Check results columns
    res_df = pd.read_csv(res_csv)
    for col in ["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"]:
        assert col in res_df.columns

    # Verify strong upregulation for Gene_01 and downregulation for Gene_02
    gene_1_res = res_df[res_df["gene_id"] == "Gene_01"].iloc[0]
    gene_2_res = res_df[res_df["gene_id"] == "Gene_02"].iloc[0]

    assert gene_1_res["log2FoldChange"] > 3.0
    assert gene_1_res["padj"] < 0.05

    assert gene_2_res["log2FoldChange"] < -3.0
    assert gene_2_res["padj"] < 0.05

    summary = json.loads(summary_json)
    assert summary["total_genes_tested"] == 20
    assert summary["significant_upregulated"] >= 1
    assert summary["significant_downregulated"] >= 1
    assert summary["total_significant_degs"] >= 2


def test_deseq2_interoperability_with_volcano_and_ma_plots(tmp_path):
    mod_deseq = load_module(DESEQ2_MODULE, "standalone_deseq2")
    mod_plots = load_module(PLOTS_MODULE, "standalone_plots")

    node_deseq = mod_deseq.DESeq2()
    volcano_node = mod_plots.VolcanoPlot()
    ma_node = mod_plots.MaPlot()

    samples = ["C1", "C2", "T1", "T2"]
    genes = [f"G_{i}" for i in range(10)]
    data = []
    for g in genes:
        if g == "G_0":
            data.append([g, 10, 10, 500, 500])
        elif g == "G_1":
            data.append([g, 500, 500, 10, 10])
        else:
            data.append([g, 100, 100, 100, 100])

    counts_file = tmp_path / "counts_small.csv"
    pd.DataFrame(data, columns=["gene_id"] + samples).to_csv(counts_file, index=False)

    meta_file = tmp_path / "meta_small.csv"
    pd.DataFrame({
        "sample": samples,
        "group": ["ctrl", "ctrl", "treat", "treat"],
    }).to_csv(meta_file, index=False)

    res_csv, norm_csv, _ = node_deseq.run(
        count_matrix_csv=str(counts_file),
        sample_metadata_csv=str(meta_file),
        output_dir=str(tmp_path / "deseq_res"),
        condition_col="group",
        contrast_reference="ctrl",
        contrast_target="treat",
    )

    # Feed deseq2_results.csv directly into VolcanoPlot
    volcano_png, volcano_json = volcano_node.run(
        deg_table_path=res_csv,
        output_dir=str(tmp_path / "plots"),
        title="DESeq2 Volcano",
    )
    assert Path(volcano_png).is_file()
    assert Path(volcano_png).stat().st_size > 1000

    # Feed deseq2_results.csv directly into MaPlot
    ma_png, ma_json = ma_node.run(
        deg_table_path=res_csv,
        output_dir=str(tmp_path / "plots"),
        title="DESeq2 MA Plot",
    )
    assert Path(ma_png).is_file()
    assert Path(ma_png).stat().st_size > 1000


def test_deseq2_rejects_missing_file():
    mod = load_module(DESEQ2_MODULE, "standalone_deseq2")
    with pytest.raises(FileNotFoundError, match="Gene Count Matrix is not a file"):
        mod.DESeq2().run(count_matrix_csv="nonexistent_counts.csv", sample_metadata_csv="nonexistent_meta.csv")
