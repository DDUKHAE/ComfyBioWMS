import importlib.util
from pathlib import Path
import pytest

MODULE = Path("nodes/class_1/plots.py").resolve()


def load_module():
    spec = importlib.util.spec_from_file_location("standalone_plots", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_plots_file_loads_standalone_and_owns_mappings():
    module = load_module()
    assert len(module.NODE_CLASS_MAPPINGS) == 20
    assert set(module.NODE_CLASS_MAPPINGS) == set(module.NODE_DISPLAY_NAME_MAPPINGS)
    for name, cls in module.NODE_CLASS_MAPPINGS.items():
        assert hasattr(cls, "INPUT_TYPES")
        assert hasattr(cls, "RETURN_TYPES")
        assert cls.RETURN_TYPES == ("STRING", "IMAGE")
        assert cls.FUNCTION == "run"
        assert cls.CATEGORY == "ComfyBIO/Visualization"


def test_volcano_plot_runs_on_real_deg_csv(tmp_path):
    module = load_module()
    csv_file = tmp_path / "deg.csv"
    csv_file.write_text(
        "gene,log2FoldChange,padj\n"
        "TP53,2.5,0.001\n"
        "EGFR,-3.1,0.0002\n"
        "ACTB,0.1,0.85\n"
        "GAPDH,-0.05,0.92\n"
    )
    node = module.VolcanoPlot()
    out_img, tensor = node.run(str(csv_file), output_dir=str(tmp_path / "out"))
    assert Path(out_img).is_file()
    assert tensor is not None


def test_ma_plot_runs_on_real_csv(tmp_path):
    module = load_module()
    csv_file = tmp_path / "deg.csv"
    csv_file.write_text(
        "gene,baseMean,log2FoldChange,padj\n"
        "TP53,1500,2.5,0.001\n"
        "EGFR,800,-3.1,0.0002\n"
        "ACTB,50000,0.1,0.85\n"
    )
    node = module.MaPlot()
    out_img, tensor = node.run(str(csv_file), output_dir=str(tmp_path / "out"))
    assert Path(out_img).is_file()
    assert tensor is not None


def test_clustermap_runs_on_real_matrix(tmp_path):
    module = load_module()
    csv_file = tmp_path / "counts.csv"
    csv_file.write_text(
        "gene,S1,S2,S3\n"
        "G1,10.5,12.3,15.1\n"
        "G2,5.2,4.8,6.1\n"
        "G3,100.2,95.4,110.0\n"
    )
    node = module.ClustermapHeatmap()
    out_img, tensor = node.run(str(csv_file), output_dir=str(tmp_path / "out"))
    assert Path(out_img).is_file()
    assert tensor is not None


def test_kaplan_meier_runs_on_real_clinical_data(tmp_path):
    module = load_module()
    csv_file = tmp_path / "clinical.csv"
    csv_file.write_text(
        "patient,time_months,vital_status,biomarker_strata\n"
        "P1,12,Dead,High\n"
        "P2,24,Alive,High\n"
        "P3,6,Dead,Low\n"
        "P4,18,Dead,Low\n"
    )
    node = module.KaplanMeierSurvival()
    out_img, tensor = node.run(str(csv_file), output_dir=str(tmp_path / "out"))
    assert Path(out_img).is_file()
    assert tensor is not None


def test_plots_reject_missing_files():
    module = load_module()
    with pytest.raises(FileNotFoundError, match="DEG Results Table is not a file"):
        module.VolcanoPlot().run("non_existent_file.csv")
