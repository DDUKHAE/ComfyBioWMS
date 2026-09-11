"""Registration checks for the verified standalone node catalog."""

import inspect
import pytest

from nodes.class_1 import CLASS_1_NODE_MAPPINGS
from nodes.class_2 import CLASS_2_NODE_MAPPINGS


def test_registered_nodes_inventory():
    assert len(CLASS_1_NODE_MAPPINGS) == 33
    assert len(CLASS_2_NODE_MAPPINGS) == 57


@pytest.mark.parametrize("mappings", [CLASS_1_NODE_MAPPINGS, CLASS_2_NODE_MAPPINGS])
def test_registered_nodes_follow_comfyui_contract(mappings):
    for name, node in mappings.items():
        assert inspect.isclass(node), name
        inputs = node.INPUT_TYPES()
        assert "required" in inputs, f"{name} missing required in INPUT_TYPES"
        assert isinstance(node.RETURN_TYPES, tuple), f"{name} RETURN_TYPES not tuple"
        assert isinstance(node.FUNCTION, str), f"{name} FUNCTION not str"
        assert callable(getattr(node, node.FUNCTION)), f"{name} missing callable function {node.FUNCTION}"
        assert isinstance(node.CATEGORY, str), f"{name} CATEGORY not str"
        assert node.CATEGORY.startswith("ComfyBIO/"), f"{name} invalid CATEGORY"


def test_every_registered_cli_node_offers_extra_command():
    for name, node in CLASS_2_NODE_MAPPINGS.items():
        if name != "ESMFold":  # Pure Python neural model
            assert "extra_command" in node.INPUT_TYPES()["optional"], name

