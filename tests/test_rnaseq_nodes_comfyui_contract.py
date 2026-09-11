import importlib
import inspect
from pathlib import Path
import pytest

RNASEQ_CLASS_1_FILES = [
    "nodes.class_1.tximport",
    "nodes.class_1.rnaseq_pca",
]

RNASEQ_CLASS_2_FILES = [
    "nodes.class_2.trim_galore",
    "nodes.class_2.trimmomatic",
    "nodes.class_2.star",
    "nodes.class_2.salmon",
    "nodes.class_2.kallisto",
    "nodes.class_2.multiqc",
    "nodes.class_2.cat_fastq",
    "nodes.class_2.strandedness",
    "nodes.class_2.umi_tools",
    "nodes.class_2.picard",
    "nodes.class_2.sortmerna",
    "nodes.class_2.bbtools",
    "nodes.class_2.stringtie",
    "nodes.class_2.bedtools",
    "nodes.class_2.ucsc_tools",
    "nodes.class_2.rseqc",
    "nodes.class_2.qualimap",
    "nodes.class_2.preseq",
    "nodes.class_2.dupradar",
    "nodes.class_2.kraken2",
    "nodes.class_2.bracken",
    "nodes.class_2.sylph",
    "nodes.class_2.hisat2",
    "nodes.class_2.rsem",
]

@pytest.mark.parametrize("mod_name", RNASEQ_CLASS_1_FILES + RNASEQ_CLASS_2_FILES)
def test_rnaseq_modules_export_mappings(mod_name):
    mod = importlib.import_module(mod_name)
    assert hasattr(mod, "NODE_CLASS_MAPPINGS"), f"{mod_name} missing NODE_CLASS_MAPPINGS"
    assert hasattr(mod, "NODE_DISPLAY_NAME_MAPPINGS"), f"{mod_name} missing NODE_DISPLAY_NAME_MAPPINGS"
    assert set(mod.NODE_CLASS_MAPPINGS) == set(mod.NODE_DISPLAY_NAME_MAPPINGS)

@pytest.mark.parametrize("mod_name", RNASEQ_CLASS_1_FILES + RNASEQ_CLASS_2_FILES)
def test_rnaseq_nodes_comply_with_comfyui_spec(mod_name):
    mod = importlib.import_module(mod_name)
    for name, node_cls in mod.NODE_CLASS_MAPPINGS.items():
        assert inspect.isclass(node_cls), f"{name} is not a class"
        assert hasattr(node_cls, "INPUT_TYPES"), f"{name} missing INPUT_TYPES"
        inputs = node_cls.INPUT_TYPES()
        assert "required" in inputs, f"{name} INPUT_TYPES missing required"
        assert hasattr(node_cls, "RETURN_TYPES"), f"{name} missing RETURN_TYPES"
        assert isinstance(node_cls.RETURN_TYPES, tuple), f"{name} RETURN_TYPES not tuple"
        assert hasattr(node_cls, "FUNCTION"), f"{name} missing FUNCTION"
        assert hasattr(node_cls, node_cls.FUNCTION), f"{name} missing function {node_cls.FUNCTION}"
        assert hasattr(node_cls, "CATEGORY"), f"{name} missing CATEGORY"
        assert node_cls.CATEGORY.startswith("ComfyBIO/"), f"{name} invalid CATEGORY {node_cls.CATEGORY}"

@pytest.mark.parametrize("mod_name", RNASEQ_CLASS_2_FILES)
def test_rnaseq_cli_nodes_have_extra_command_if_configurable(mod_name):
    mod = importlib.import_module(mod_name)
    for name, node_cls in mod.NODE_CLASS_MAPPINGS.items():
        # Nodes that take CLI options should offer extra_command
        inputs = node_cls.INPUT_TYPES()
        optional = inputs.get("optional", {})
        # Note: CatFastq and BedGraphToBigWig have fixed syntax with no extra CLI flags needed
        if name not in ("CatFastq", "BedGraphToBigWig", "RSeQCInferExperiment"):
            assert "extra_command" in optional, f"{name} missing extra_command in optional"
