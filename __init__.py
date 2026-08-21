"""ComfyBIOWMS: High-Throughput Bioinformatics Workflow Management System for ComfyUI.

This module exposes the root entry point for ComfyUI custom nodes discovery,
including node class mappings, display names, and web frontend directory.
"""

import logging
from pathlib import Path

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
