"""Scanpy single-cell RNA-seq analysis nodes.

Python packages: scanpy, anndata, numpy, matplotlib
External binaries: none
"""

import json
from pathlib import Path

try:
    from ..compat import ensure_pandas_compat
except Exception:
    try:
        from nodes.compat import ensure_pandas_compat
    except Exception:
        def ensure_pandas_compat():
            pass

ensure_pandas_compat()


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path

def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and custom_dir.strip():
        return Path(custom_dir).expanduser().resolve()
    try:
        folder_paths = __import__("folder_paths")
        return Path(folder_paths.get_output_directory()) / node_name
    except Exception:
        return Path.cwd() / "ComfyUI" / "output" / node_name


class ScanpyQC:
    OUTPUT_NODE = True
    CATEGORY = "ComfyBIO/Single-Cell"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("filtered_h5ad", "qc_metrics_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_h5ad": ("STRING", {"default": ""}),
            },
            "optional": {
                "min_genes": ("INT", {"default": 200, "min": 0, "max": 10000}),
                "min_cells": ("INT", {"default": 3, "min": 0, "max": 1000}),
                "max_percent_mito": ("FLOAT", {"default": 20.0, "min": 0.0, "max": 100.0, "step": 0.5}),
            },
        }

    def run(self, input_h5ad: str, output_dir: str = "", min_genes: int = 200, min_cells: int = 3, max_percent_mito: float = 20.0):
        ensure_pandas_compat()
        import scanpy as sc
        h5ad_path = _file(input_h5ad, "Input AnnData (.h5ad)")
        out = _output_dir("ScanpyQC", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        adata = sc.read_h5ad(h5ad_path)
        n_cells_initial = adata.n_obs
        n_genes_initial = adata.n_vars

        if min_genes > 0:
            sc.pp.filter_cells(adata, min_genes=min_genes)
        if min_cells > 0:
            sc.pp.filter_genes(adata, min_cells=min_cells)

        # Mitochondrial genes
        adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

        if "pct_counts_mt" in adata.obs.columns:
            adata = adata[adata.obs["pct_counts_mt"] < max_percent_mito, :].copy()

        out_file = out / f"{h5ad_path.stem}_qc_filtered.h5ad"
        adata.write(out_file)

        summary = {
            "n_cells_initial": n_cells_initial,
            "n_genes_initial": n_genes_initial,
            "n_cells_filtered": adata.n_obs,
            "n_genes_filtered": adata.n_vars,
            "min_genes": min_genes,
            "min_cells": min_cells,
            "max_percent_mito": max_percent_mito,
        }

        return (str(out_file), json.dumps(summary))


class ScanpyNormalize:
    OUTPUT_NODE = True
    CATEGORY = "ComfyBIO/Single-Cell"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("normalized_h5ad",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_h5ad": ("STRING", {"default": ""}),
            },
            "optional": {
                "target_sum": ("FLOAT", {"default": 1e4, "min": 1e2, "max": 1e7}),
                "n_top_genes": ("INT", {"default": 2000, "min": 500, "max": 10000}),
            },
        }

    def run(self, input_h5ad: str, output_dir: str = "", target_sum: float = 1e4, n_top_genes: int = 2000):
        ensure_pandas_compat()
        import scanpy as sc
        h5ad_path = _file(input_h5ad, "Input AnnData (.h5ad)")
        out = _output_dir("ScanpyNormalize", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        adata = sc.read_h5ad(h5ad_path)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars))

        out_file = out / f"{h5ad_path.stem}_normalized.h5ad"
        adata.write(out_file)
        return (str(out_file),)


class ScanpyCluster:
    OUTPUT_NODE = True
    CATEGORY = "ComfyBIO/Single-Cell"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("clustered_h5ad", "umap_plot_png")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_h5ad": ("STRING", {"default": ""}),
            },
            "optional": {
                "n_pcs": ("INT", {"default": 30, "min": 5, "max": 100}),
                "resolution": ("FLOAT", {"default": 0.5, "min": 0.1, "max": 3.0, "step": 0.1}),
            },
        }

    def run(self, input_h5ad: str, output_dir: str = "", n_pcs: int = 30, resolution: float = 0.5):
        ensure_pandas_compat()
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import scanpy as sc

        h5ad_path = _file(input_h5ad, "Input AnnData (.h5ad)")
        out = _output_dir("ScanpyCluster", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        adata = sc.read_h5ad(h5ad_path)
        sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1, adata.n_obs - 1), svd_solver="arpack")
        sc.pp.neighbors(adata, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=resolution, flavor="igraph")

        umap_png = out / f"{h5ad_path.stem}_umap.png"
        fig = sc.pl.umap(adata, color=["leiden"], show=False, return_fig=True)
        fig.savefig(umap_png, dpi=200, bbox_inches="tight")
        plt.close(fig)

        out_file = out / f"{h5ad_path.stem}_clustered.h5ad"
        adata.write(out_file)
        return (str(out_file), str(umap_png))


NODE_CLASS_MAPPINGS = {
    "ScanpyQC": ScanpyQC,
    "ScanpyNormalize": ScanpyNormalize,
    "ScanpyCluster": ScanpyCluster,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "ScanpyQC": "Scanpy: Single-cell Quality Control",
    "ScanpyNormalize": "Scanpy: Normalize & Highly Variable Genes",
    "ScanpyCluster": "Scanpy: PCA, UMAP & Leiden Clustering",
}


# Backward compatibility aliases
ScanpyQCNode = ScanpyQC
ScanpyNormalizeNode = ScanpyNormalize
ScanpyClusterNode = ScanpyCluster
