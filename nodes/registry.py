"""Registry containing only nodes with official-data E2E verification."""

from .class_1 import CLASS_1_DISPLAY_NAME_MAPPINGS, CLASS_1_NODE_MAPPINGS
from .class_2 import CLASS_2_DISPLAY_NAME_MAPPINGS, CLASS_2_NODE_MAPPINGS

NODE_CLASS_MAPPINGS = {**CLASS_1_NODE_MAPPINGS, **CLASS_2_NODE_MAPPINGS}
NODE_DISPLAY_NAME_MAPPINGS = {
    **CLASS_1_DISPLAY_NAME_MAPPINGS,
    **CLASS_2_DISPLAY_NAME_MAPPINGS,
}
NODE_CLASS_ALIASES = {f"{k}Node": v for k, v in NODE_CLASS_MAPPINGS.items()}


def resolve_node_class(node_type):
    if node_type in NODE_CLASS_MAPPINGS:
        return NODE_CLASS_MAPPINGS[node_type]
    if node_type in NODE_CLASS_ALIASES:
        return NODE_CLASS_ALIASES[node_type]
    raise KeyError(f"No verified ComfyBIO node registered for {node_type}")


__all__ = [
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "NODE_CLASS_ALIASES",
    "CLASS_1_NODE_MAPPINGS",
    "CLASS_2_NODE_MAPPINGS",
    "resolve_node_class",
]
