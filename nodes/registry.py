"""Registry containing only nodes with official-data E2E verification."""

from .class_1 import CLASS_1_DISPLAY_NAME_MAPPINGS, CLASS_1_NODE_MAPPINGS
from .class_2 import CLASS_2_DISPLAY_NAME_MAPPINGS, CLASS_2_NODE_MAPPINGS

NODE_CLASS_MAPPINGS = {**CLASS_1_NODE_MAPPINGS, **CLASS_2_NODE_MAPPINGS}
NODE_DISPLAY_NAME_MAPPINGS = {
    **CLASS_1_DISPLAY_NAME_MAPPINGS,
    **CLASS_2_DISPLAY_NAME_MAPPINGS,
}
NODE_CLASS_ALIASES = {}


def resolve_node_class(node_type):
    try:
        return NODE_CLASS_MAPPINGS[node_type]
    except KeyError:
        raise KeyError(f"No verified ComfyBIO node registered for {node_type}") from None


__all__ = [
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "NODE_CLASS_ALIASES",
    "CLASS_1_NODE_MAPPINGS",
    "CLASS_2_NODE_MAPPINGS",
    "resolve_node_class",
]
