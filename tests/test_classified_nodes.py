"""Registration checks for the official-data verified node set."""

import inspect

import pytest

from nodes.class_1 import CLASS_1_NODE_MAPPINGS
from nodes.class_2 import CLASS_2_NODE_MAPPINGS


def test_only_verified_nodes_are_registered():
    assert set(CLASS_1_NODE_MAPPINGS) == {
        "BiopythonSeqIOStatsNode",
        "BiopythonAlignmentStatsNode",
    }
    assert set(CLASS_2_NODE_MAPPINGS) == {
        "FastpNode",
        "FastQCNode",
        "BwaMem2IndexNode",
        "BwaMem2AlignNode",
        "SamtoolsSortNode",
        "SamtoolsIndexNode",
        "SamtoolsMarkdupNode",
        "BcftoolsMpileupNode",
        "BcftoolsCallNode",
        "BcftoolsFilterNode",
        "SpadesNode",
        "QuastNode",
    }


@pytest.mark.parametrize("mappings", [CLASS_1_NODE_MAPPINGS, CLASS_2_NODE_MAPPINGS])
def test_registered_nodes_follow_comfyui_contract(mappings):
    for name, node in mappings.items():
        assert inspect.isclass(node), name
        inputs = node.INPUT_TYPES()
        assert "required" in inputs
        assert isinstance(node.RETURN_TYPES, tuple)
        assert isinstance(node.FUNCTION, str)
        assert callable(getattr(node, node.FUNCTION))
        assert isinstance(node.CATEGORY, str)


def test_every_registered_cli_node_offers_extra_command():
    for name, node in CLASS_2_NODE_MAPPINGS.items():
        assert "extra_command" in node.INPUT_TYPES()["optional"], name
