"""ComfyBIOWMS: High-Throughput Bioinformatics Workflow Management System for ComfyUI.

This module exposes the root entry point for ComfyUI custom nodes discovery,
including node class mappings, display names, and web frontend directory.
"""

import logging
import sys
from pathlib import Path

# Ensure repository root and engine/src are on sys.path for internal module resolution (e.g. bioflow runtime)
_ROOT_DIR = Path(__file__).resolve().parent
_ENGINE_SRC = _ROOT_DIR / "engine" / "src"

if str(_ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(_ROOT_DIR))
if str(_ENGINE_SRC) not in sys.path:
    sys.path.insert(0, str(_ENGINE_SRC))

# Apply backwards compatibility shims (e.g. pandas 2.x/3.x legacy unpickling support)
try:
    from .nodes import compat  # noqa: F401
except Exception:
    try:
        import nodes.compat  # noqa: F401
    except Exception:
        pass

logger = logging.getLogger("ComfyBIOWMS")

# Import all registered node classes and display names from the registry
try:
    from .nodes.registry import (
        NODE_CLASS_MAPPINGS,
        NODE_DISPLAY_NAME_MAPPINGS,
    )
except (ImportError, ValueError):
    from nodes.registry import (
        NODE_CLASS_MAPPINGS,
        NODE_DISPLAY_NAME_MAPPINGS,
    )

# Specify directory containing frontend extensions (LiteGraph custom widgets, styles, etc.)
WEB_DIRECTORY = "./web"

__all__ = [
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "WEB_DIRECTORY",
]

# Log initialization banner
try:
    total_nodes = len(NODE_CLASS_MAPPINGS)
    print(f"\033[32m[ComfyBIOWMS]\033[0m Successfully loaded {total_nodes} bioinformatics & CADD custom nodes.")
except Exception as e:
    logger.warning("ComfyBIOWMS node count logging error: %s", e)
