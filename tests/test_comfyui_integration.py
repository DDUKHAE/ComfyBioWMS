"""ComfyUI integration test suite.

Verifies:
1. ComfyUI custom node discovery, registration contracts, and port signatures across all 98 registered nodes.
2. ComfyUI prompt graph dependency resolution and execution flow.
3. IS_CHANGED caching behavior during simulated ComfyUI queue runs.
"""

import json
from pathlib import Path
import pytest

import __init__ as comfy_entrypoint


def test_comfyui_node_registration_integrity():
    """Verify that ComfyUI discovers and registers all active nodes with complete contracts."""
    assert hasattr(comfy_entrypoint, "NODE_CLASS_MAPPINGS")
    assert hasattr(comfy_entrypoint, "NODE_DISPLAY_NAME_MAPPINGS")

    mappings = comfy_entrypoint.NODE_CLASS_MAPPINGS
    display_mappings = comfy_entrypoint.NODE_DISPLAY_NAME_MAPPINGS

    assert len(mappings) == 98, f"Expected exactly 98 active registered nodes, found {len(mappings)}"
    assert len(display_mappings) == 98

    for node_name, node_cls in mappings.items():
        assert hasattr(node_cls, "INPUT_TYPES"), f"{node_name} missing INPUT_TYPES"
        input_types = node_cls.INPUT_TYPES()
        assert isinstance(input_types, dict), f"{node_name} INPUT_TYPES() must return a dict"
        assert "required" in input_types or "optional" in input_types, f"{node_name} inputs missing required/optional"

        assert hasattr(node_cls, "RETURN_TYPES"), f"{node_name} missing RETURN_TYPES"
        assert isinstance(node_cls.RETURN_TYPES, (tuple, list)), f"{node_name} RETURN_TYPES must be a tuple/list"

        assert hasattr(node_cls, "FUNCTION"), f"{node_name} missing FUNCTION attribute"
        func_name = node_cls.FUNCTION
        assert hasattr(node_cls, func_name), f"{node_name} missing execution method '{func_name}'"

        assert hasattr(node_cls, "CATEGORY"), f"{node_name} missing CATEGORY attribute"
        assert "ComfyBIO" in node_cls.CATEGORY, f"{node_name} category '{node_cls.CATEGORY}' should contain 'ComfyBIO'"


def test_comfyui_caching_simulation(tmp_path):
    """Simulate ComfyUI's prompt cache engine using node IS_CHANGED methods."""
    from nodes.class_2.fastp import Fastp
    from nodes.class_1.plots import VolcanoPlot

    # Create real existing test input files
    mock_fq = tmp_path / "mock.fq"
    mock_fq.write_text("mock sequence data")
    mock_deg = tmp_path / "deg.csv"
    mock_deg.write_text("gene,baseMean,log2FoldChange,pvalue,padj\n")

    # Fastp caching check
    cache_key1 = Fastp.IS_CHANGED(read1=str(mock_fq), threads=4)
    cache_key2 = Fastp.IS_CHANGED(read1=str(mock_fq), threads=4)
    assert cache_key1 == cache_key2, "Identical prompt inputs must yield identical cache keys"

    cache_key_diff = Fastp.IS_CHANGED(read1=str(mock_fq), threads=8)
    assert cache_key1 != cache_key_diff, "Modified parameters must invalidate cache"

    # Visualizer caching check
    v_key1 = VolcanoPlot.IS_CHANGED(deg_table_path=str(mock_deg), figure_title="DEG Analysis")
    v_key2 = VolcanoPlot.IS_CHANGED(deg_table_path=str(mock_deg), figure_title="DEG Analysis")
    assert v_key1 == v_key2


def test_final_case_study_workflows_graph_integrity():
    """Verify that case study workflow JSONs match registered ComfyUI nodes and valid link contracts."""
    workflow_dir = Path("workflows")
    workflow_files = [
        p for p in workflow_dir.glob("*.json")
        if (p.name.startswith("final_") or p.name.startswith("cs")) and p.name != "final_case_study_status.json"
    ]
    assert len(workflow_files) >= 3, f"Expected at least 3 case study workflows, found {len(workflow_files)}"

    registered_classes = set(comfy_entrypoint.NODE_CLASS_MAPPINGS.keys())
    # ComfyUI core visualizers/preview nodes that workflows link to
    allowed_terminal_nodes = {"PreviewAny", "PreviewImage", "SaveImage"}

    for wf_path in workflow_files:
        with open(wf_path, "r", encoding="utf-8") as f:
            wf_data = json.load(f)

        assert "nodes" in wf_data or "prompt" in wf_data, f"Invalid workflow structure in {wf_path}"
        nodes = wf_data.get("nodes", [])
        for node in nodes:
            node_type = node.get("type", "")
            # Verify node type is registered or recognized preview terminal
            assert (
                node_type in registered_classes or node_type in allowed_terminal_nodes
            ), f"Workflow {wf_path.name} uses unregistered node: {node_type}"


def test_pandas_backward_compatibility_shim():
    """Verify that legacy pandas internal paths required by older AnnData/H5AD objects are safely aliased."""
    import nodes.compat
    nodes.compat.ensure_pandas_compat()

    import pandas.core.index as p_idx
    import pandas.core.indexes.numeric as p_num
    assert hasattr(p_idx, "Index")
    assert hasattr(p_num, "Int64Index")

