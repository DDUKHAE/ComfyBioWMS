"""AnnData IO and Inspection node.

Python packages: anndata
External binaries: none
"""

import json
from pathlib import Path


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


class AnnDataInspect:
    CATEGORY = "ComfyBIO/Single-Cell"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("summary_json", "n_cells", "n_genes")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "h5ad_file": ("STRING", {"default": ""}),
            },
        }

    def run(self, h5ad_file: str):
        import anndata as ad
        path = _file(h5ad_file, "AnnData (.h5ad)")
        adata = ad.read_h5ad(str(path), backed="r")
        summary = {
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_keys": list(adata.obs.keys()),
            "var_keys": list(adata.var.keys()),
            "obsm_keys": list(adata.obsm.keys()),
            "layers_keys": list(adata.layers.keys()) if hasattr(adata, "layers") else [],
        }
        return (json.dumps(summary, indent=2), int(adata.n_obs), int(adata.n_vars))


NODE_CLASS_MAPPINGS = {"AnnDataInspect": AnnDataInspect}
NODE_DISPLAY_NAME_MAPPINGS = {"AnnDataInspect": "AnnData: Inspect H5AD Dataset"}


# Backward compatibility aliases
AnnDataInspectNode = AnnDataInspect

__all__ = ["AnnDataInspect",
    "AnnDataInspectNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
