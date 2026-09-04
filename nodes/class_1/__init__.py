"""Verified standalone in-process nodes."""

from .biopython import (
    BiopythonAlignmentStatsNode,
    BiopythonSeqIOStatsNode,
    NODE_CLASS_MAPPINGS,
    NODE_DISPLAY_NAME_MAPPINGS,
)

CLASS_1_NODE_MAPPINGS = NODE_CLASS_MAPPINGS
CLASS_1_DISPLAY_NAME_MAPPINGS = NODE_DISPLAY_NAME_MAPPINGS

__all__ = [
    "BiopythonAlignmentStatsNode",
    "BiopythonSeqIOStatsNode",
    "CLASS_1_NODE_MAPPINGS",
    "CLASS_1_DISPLAY_NAME_MAPPINGS",
]
