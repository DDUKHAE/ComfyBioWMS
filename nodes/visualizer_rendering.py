"""Rendering and figure-to-tensor utilities for Publication Visualizers."""

import io
from pathlib import Path
from typing import Any, Optional, Tuple, Union


def configure_publication_style(
    style: str = "nature",
    font_size: int = 10,
    dpi: int = 300,
) -> None:
    """Apply publication-grade styling to Matplotlib RC parameters."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "Liberation Sans"],
            "font.size": font_size,
            "axes.labelsize": font_size + 1,
            "axes.titlesize": font_size + 2,
            "xtick.labelsize": font_size - 1,
            "ytick.labelsize": font_size - 1,
            "legend.fontsize": font_size - 1,
            "figure.titlesize": font_size + 3,
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "savefig.bbox": "tight",
            "axes.linewidth": 1.0,
            "grid.linewidth": 0.5,
            "lines.linewidth": 1.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        })
    except Exception:
        pass


def figure_to_image_tensor(fig: Any) -> Any:
    """Convert a Matplotlib figure to a ComfyUI standard IMAGE tensor [1, H, W, 3] (float32, 0.0-1.0)."""
    try:
        import numpy as np
        from PIL import Image
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=getattr(fig, "dpi", 300) or 300)
        buf.seek(0)
        pil_image = Image.open(buf).convert("RGB")
        np_arr = np.asarray(pil_image, dtype=np.float32) / 255.0
        buf.close()

        try:
            import torch
            return torch.from_numpy(np_arr)[None, :]
        except Exception:
            return np_arr[None, :]
    except Exception:
        return None


def save_and_tensorize_figure(
    fig: Any,
    output_path: Union[str, Path],
    close_fig: bool = True,
    dpi: int = 300,
) -> Tuple[str, Any]:
    """Save a Matplotlib figure to disk and return its path and ComfyUI IMAGE tensor."""
    try:
        import matplotlib.pyplot as plt
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(out_path), format="png", bbox_inches="tight", dpi=dpi)
        tensor = figure_to_image_tensor(fig)
        if close_fig:
            plt.close(fig)
        return str(out_path), tensor
    except Exception:
        return str(output_path), None
