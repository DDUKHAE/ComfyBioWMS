"""Runnable ComfyBIO custom-node examples."""
import sys
from pathlib import Path

_ROOT_DIR = Path(__file__).resolve().parents[1]
_ENGINE_SRC = _ROOT_DIR / "engine" / "src"

if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
if str(_ENGINE_SRC) not in sys.path:
    sys.path.insert(0, str(_ENGINE_SRC))

from .registry import (
    NODE_CLASS_MAPPINGS,
    NODE_DISPLAY_NAME_MAPPINGS,
    CLASS_1_NODE_MAPPINGS,
    CLASS_2_NODE_MAPPINGS,
)

__all__ = [
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "CLASS_1_NODE_MAPPINGS",
    "CLASS_2_NODE_MAPPINGS",
]


