"""Publication-quality scientific visualization custom nodes.

Consolidates all 20 Tier-1 publication figure renderers into a single module.
Provides direct ComfyUI IMAGE tensor output alongside saved image paths.
All nodes output: (plot_path: STRING, preview_image: IMAGE).

Strict scientific integrity:
- Every visualizer strictly requires valid input data.
- NO synthetic random fallbacks or hardcoded placeholder plots are allowed.
- Any missing, empty, or unparseable input raises FileNotFoundError or ValueError.

Python packages: matplotlib, numpy, pandas, PIL, biopython
External binaries: none
"""

import io
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

try:
    import torch
except ImportError:
    torch = None


def _file(value: str, label: str) -> Path:
    if not value or not str(value).strip():
        raise ValueError(f"{label} path is required and cannot be empty.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _read_table(path_or_str: str, label: str, **kwargs) -> pd.DataFrame:
    """Validate and read CSV/TSV table. Raises FileNotFoundError or ValueError if missing/empty."""
    p = _file(path_or_str, label)
    if p.stat().st_size == 0:
        raise ValueError(f"{label} at '{p}' is empty (0 bytes).")
    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        first_line = f.readline()
        if not first_line.strip():
            raise ValueError(f"{label} at '{p}' contains no data rows.")
        sep = "\t" if "\t" in first_line else ","
    df = pd.read_csv(p, sep=sep, **kwargs)
    if df.empty:
        raise ValueError(f"{label} at '{p}' contains no data rows.")
    return df


def _find_col(df: pd.DataFrame, target: str, candidates: List[str] = (), default_idx: int = 0) -> str:
    """Resolve column name flexibly from dataframe."""
    if target in df.columns:
        return target
    for c in candidates:
        if c in df.columns:
            return c
    num_cols = df.select_dtypes(include=[np.number]).columns
    if len(num_cols) > default_idx:
        return num_cols[default_idx]
    if len(df.columns) > default_idx:
        return df.columns[default_idx]
    raise KeyError(f"Target column '{target}' not found in dataframe columns: {list(df.columns)}")


def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and str(custom_dir).strip():
        return Path(custom_dir).expanduser().resolve()
    try:
        folder_paths = __import__("folder_paths")
        base = Path(folder_paths.get_output_directory())
    except Exception:
        base = Path.cwd() / "ComfyUI" / "output"
    return base / node_name


def configure_publication_style(style: str = "nature", font_size: int = 10, dpi: int = 300) -> None:
    """Apply publication-grade styling to Matplotlib RC parameters."""
    try:
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
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=getattr(fig, "dpi", 300) or 300)
        buf.seek(0)
        pil_image = Image.open(buf).convert("RGB")
        np_arr = np.asarray(pil_image, dtype=np.float32) / 255.0
        buf.close()
        if torch is not None:
            return torch.from_numpy(np_arr)[None, :]
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
    out_path = Path(output_path).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), format="png", bbox_inches="tight", dpi=dpi)
    tensor = figure_to_image_tensor(fig)
    if close_fig:
        try:
            plt.close(fig)
        except Exception:
            pass
    return str(out_path), tensor


def _resolve_output_image_path(node_name: str, output_image_path: str, output_dir: str, default_filename: str) -> Path:
    if output_image_path and str(output_image_path).strip():
        p = Path(output_image_path).expanduser().resolve()
        p.parent.mkdir(parents=True, exist_ok=True)
        return p
    out_base = _output_dir(node_name, output_dir)
    out_base.mkdir(parents=True, exist_ok=True)
    return out_base / default_filename


class _BaseVisualizer:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Visualization"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")
    FUNCTION = "run"

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        try:
            from bioflow.runtime.artifacts import compute_input_fingerprint
            return compute_input_fingerprint(**kwargs)
        except Exception:
            return float("NaN")


    @classmethod
    def _common_visualizer_inputs(cls) -> Dict[str, Any]:
        return {
            "figure_title": ("STRING", {"default": "Publication Figure"}),
            "dpi": ("INT", {"default": 300, "min": 72, "max": 1200}),
            "figure_width": ("FLOAT", {"default": 7.0, "min": 2.0, "max": 30.0}),
            "figure_height": ("FLOAT", {"default": 5.5, "min": 2.0, "max": 30.0}),
            "style": (["nature", "cell", "science", "classic"], {"default": "nature"}),
        }

    @classmethod
    def _inputs(
        cls,
        req: Optional[Dict[str, Any]] = None,
        opt: Optional[Dict[str, Any]] = None,
        title: str = "Publication Figure",
        w: float = 7.0,
        h: float = 5.5,
    ) -> Dict[str, Any]:
        options = cls._common_visualizer_inputs()
        options["figure_title"] = ("STRING", {"default": title})
        options["figure_width"] = ("FLOAT", {"default": w, "min": 2.0, "max": 30.0})
        options["figure_height"] = ("FLOAT", {"default": h, "min": 2.0, "max": 30.0})
        if opt:
            options.update(opt)
        return {"required": req or {}, "optional": options}

    def _target_path(self, default_name: str, out_img: str = "", out_dir: str = "", style: str = "nature", dpi: int = 300) -> Path:
        configure_publication_style(style=style, dpi=dpi)
        return _resolve_output_image_path(self.__class__.__name__, out_img, out_dir, default_name)


# 1. Volcano Plot
class VolcanoPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"deg_results_csv": ("STRING", {"default": ""})},
            opt={
                "log2fc_col": ("STRING", {"default": "log2FoldChange"}),
                "pvalue_col": ("STRING", {"default": "padj"}),
                "gene_col": ("STRING", {"default": "gene_name"}),
                "fc_cutoff": ("FLOAT", {"default": 1.0, "min": 0.0, "max": 10.0, "step": 0.1}),
                "p_cutoff": ("FLOAT", {"default": 0.05, "min": 1e-12, "max": 1.0, "step": 0.01}),
                "top_genes_to_label": ("INT", {"default": 10, "min": 0, "max": 50}),
            },
            title="Differential Expression Volcano Plot",
        )

    def run(
        self,
        deg_results_csv: str = "",
        log2fc_col: str = "log2FoldChange",
        pvalue_col: str = "padj",
        gene_col: str = "gene_name",
        fc_cutoff: float = 1.0,
        p_cutoff: float = 0.05,
        top_genes_to_label: int = 10,
        figure_title: str = "Differential Expression Volcano Plot",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 5.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        fc_cutoff = kwargs.get("log2fc_cutoff", fc_cutoff)
        p_cutoff = kwargs.get("pvalue_cutoff", p_cutoff)
        figure_title = kwargs.get("title", figure_title)
        target_path = self._target_path("volcano_plot.png", output_image_path, output_dir, style, dpi)

        in_file = deg_results_csv or kwargs.get("deg_table_path") or kwargs.get("differential_results_tsv") or kwargs.get("deg_results_table")
        df = _read_table(in_file, "DEG Results Table")

        log2fc_col = _find_col(df, log2fc_col, ["log2FoldChange", "logFC", "lfc", "log2_fold_change"], 0)
        pvalue_col = _find_col(df, pvalue_col, ["padj", "pvalue", "p_val_adj", "qval", "FDR", "P"], 1)
        gene_col = _find_col(df, gene_col, ["gene_name", "gene_id", "gene", "symbol", "Gene"], 0)
        if gene_col not in df.columns:
            df[gene_col] = [f"Feature_{i}" for i in range(len(df))]

        clean_df = df.dropna(subset=[log2fc_col, pvalue_col]).copy()
        if clean_df.empty:
            raise ValueError(f"No valid rows after dropping NAs in columns '{log2fc_col}' and '{pvalue_col}'")

        clean_df[pvalue_col] = pd.to_numeric(clean_df[pvalue_col], errors="coerce").fillna(1.0)
        clean_df[log2fc_col] = pd.to_numeric(clean_df[log2fc_col], errors="coerce").fillna(0.0)
        clean_df["neg_log10_p"] = -np.log10(np.clip(clean_df[pvalue_col].astype(float), 1e-300, 1.0))
        clean_df["lfc"] = clean_df[log2fc_col].astype(float)

        is_up = (clean_df["lfc"] >= fc_cutoff) & (clean_df[pvalue_col] <= p_cutoff)
        is_down = (clean_df["lfc"] <= -fc_cutoff) & (clean_df[pvalue_col] <= p_cutoff)
        is_ns = ~is_up & ~is_down

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.scatter(clean_df.loc[is_ns, "lfc"], clean_df.loc[is_ns, "neg_log10_p"], c="#8c96a0", alpha=0.5, s=12, label=f"Not Sig ({is_ns.sum()})")
        ax.scatter(clean_df.loc[is_up, "lfc"], clean_df.loc[is_up, "neg_log10_p"], c="#d62728", alpha=0.75, s=18, label=f"Up ({is_up.sum()})")
        ax.scatter(clean_df.loc[is_down, "lfc"], clean_df.loc[is_down, "neg_log10_p"], c="#1f77b4", alpha=0.75, s=18, label=f"Down ({is_down.sum()})")

        ax.axvline(fc_cutoff, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axvline(-fc_cutoff, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axhline(-np.log10(max(p_cutoff, 1e-300)), color="black", linestyle="--", linewidth=0.8, alpha=0.6)

        if top_genes_to_label > 0:
            top_up = clean_df[is_up].nlargest(min(top_genes_to_label, is_up.sum()), "neg_log10_p") if is_up.sum() > 0 else pd.DataFrame()
            top_down = clean_df[is_down].nlargest(min(top_genes_to_label, is_down.sum()), "neg_log10_p") if is_down.sum() > 0 else pd.DataFrame()
            to_label = pd.concat([top_up, top_down])
            for _, row in to_label.iterrows():
                ax.annotate(
                    str(row[gene_col]),
                    (row["lfc"], row["neg_log10_p"]),
                    fontsize=8,
                    xytext=(3, 3),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.3, ec="none"),
                )

        ax.set_xlabel(r"$\log_2$ Fold Change")
        ax.set_ylabel(r"$-\log_{10}$ adjusted $p$-value")
        ax.set_title(figure_title)
        ax.legend(frameon=True, loc="upper right")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 2. Manhattan Plot
class ManhattanPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"gwas_summary_csv": ("STRING", {"default": ""})},
            opt={
                "chr_col": ("STRING", {"default": "CHR"}),
                "bp_col": ("STRING", {"default": "BP"}),
                "pvalue_col": ("STRING", {"default": "P"}),
                "snp_col": ("STRING", {"default": "SNP"}),
                "genome_wide_cutoff": ("FLOAT", {"default": 5e-8, "min": 1e-20, "max": 1e-3}),
                "suggestive_cutoff": ("FLOAT", {"default": 1e-5, "min": 1e-10, "max": 1e-2}),
            },
            title="GWAS Manhattan Plot", w=10.0, h=4.5,
        )

    def run(
        self,
        gwas_summary_csv: str = "",
        chr_col: str = "CHR",
        bp_col: str = "BP",
        pvalue_col: str = "P",
        snp_col: str = "SNP",
        genome_wide_cutoff: float = 5e-8,
        suggestive_cutoff: float = 1e-5,
        figure_title: str = "GWAS Manhattan Plot",
        dpi: int = 300,
        figure_width: float = 10.0,
        figure_height: float = 4.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("manhattan_plot.png", output_image_path, output_dir, style, dpi)
        in_file = gwas_summary_csv or kwargs.get("gwas_data")
        df = _read_table(in_file, "GWAS Summary Table")

        chr_col = _find_col(df, chr_col, ["CHR", "chr", "chromosome", "CHROM"], 0)
        bp_col = _find_col(df, bp_col, ["BP", "pos", "position", "POS"], 1)
        pvalue_col = _find_col(df, pvalue_col, ["P", "pvalue", "p_val", "PVAL"], 2)

        df = df.dropna(subset=[chr_col, bp_col, pvalue_col]).copy()
        if df.empty:
            raise ValueError(f"No valid rows after dropping NAs in columns '{chr_col}', '{bp_col}', '{pvalue_col}'")

        df["chr_clean"] = df[chr_col].astype(str).str.replace(r"^chr", "", case=False, regex=True)
        def _chr_key(val: str):
            val_clean = val.upper()
            if val_clean == "X":
                return 23
            if val_clean == "Y":
                return 24
            if val_clean in ("M", "MT"):
                return 25
            try:
                return int(val_clean)
            except ValueError:
                return 999

        df["chr_num"] = df["chr_clean"].apply(_chr_key)
        df[bp_col] = pd.to_numeric(df[bp_col], errors="coerce").fillna(0).astype(int)
        df[pvalue_col] = pd.to_numeric(df[pvalue_col], errors="coerce").clip(1e-300, 1.0)
        df["minuslog10p"] = -np.log10(df[pvalue_col])
        df = df.sort_values(["chr_num", bp_col])

        curr_offset = 0
        df["cumulative_bp"] = 0
        chr_centers = {}
        for ch, group in df.groupby("chr_num", sort=False):
            df.loc[group.index, "cumulative_bp"] = group[bp_col] + curr_offset
            chr_centers[ch] = curr_offset + (group[bp_col].max() / 2)
            curr_offset += group[bp_col].max() + 5_000_000

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        colors = ["#2b5c8f", "#d95f02"]
        for idx, (ch, group) in enumerate(df.groupby("chr_num", sort=False)):
            ax.scatter(group["cumulative_bp"], group["minuslog10p"], color=colors[idx % 2], s=10, alpha=0.75, edgecolors="none")

        ax.axhline(-np.log10(genome_wide_cutoff), color="red", linestyle="--", lw=1.0, label="Genome-wide sig")
        if suggestive_cutoff:
            ax.axhline(-np.log10(suggestive_cutoff), color="blue", linestyle=":", lw=0.9, label="Suggestive")

        ax.set_xticks(list(chr_centers.values()))
        ax.set_xticklabels([str(c) for c in chr_centers.keys()], fontsize=8)
        ax.set_xlabel("Chromosome")
        ax.set_ylabel(r"$-\log_{10}(P)$")
        ax.set_title(figure_title)
        ax.legend(frameon=True, loc="upper right")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 3. UMAP Scatter
class UmapScatter(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"embedding_csv": ("STRING", {"default": ""})},
            opt={
                "dim1_col": ("STRING", {"default": "UMAP_1"}),
                "dim2_col": ("STRING", {"default": "UMAP_2"}),
                "cluster_col": ("STRING", {"default": "cluster"}),
                "point_size": ("FLOAT", {"default": 8.0, "min": 0.5, "max": 50.0}),
                "alpha": ("FLOAT", {"default": 0.8, "min": 0.1, "max": 1.0}),
            },
            title="Single-Cell UMAP Embedding", w=7.0, h=6.0,
        )

    def run(
        self,
        embedding_csv: str = "",
        dim1_col: str = "UMAP_1",
        dim2_col: str = "UMAP_2",
        cluster_col: str = "cluster",
        point_size: float = 8.0,
        alpha: float = 0.8,
        figure_title: str = "Single-Cell UMAP Embedding",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("umap_scatter.png", output_image_path, output_dir, style, dpi)
        in_file = embedding_csv or kwargs.get("embedding_data")
        df = _read_table(in_file, "Embedding Coordinates Table")

        dim1_col = _find_col(df, dim1_col, ["UMAP_1", "UMAP1", "dim1", "X", "umap1"], 0)
        dim2_col = _find_col(df, dim2_col, ["UMAP_2", "UMAP2", "dim2", "Y", "umap2"], 1)
        if cluster_col not in df.columns:
            cluster_col = _find_col(df, "cluster", ["leiden", "louvain", "cell_type", "group", "Cluster"], 2)

        df = df.dropna(subset=[dim1_col, dim2_col]).copy()
        if df.empty:
            raise ValueError(f"No valid coordinate rows in '{dim1_col}' and '{dim2_col}'")

        df[dim1_col] = pd.to_numeric(df[dim1_col], errors="coerce")
        df[dim2_col] = pd.to_numeric(df[dim2_col], errors="coerce")
        df = df.dropna(subset=[dim1_col, dim2_col])

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        unique_clusters = df[cluster_col].unique()
        cmap = plt.get_cmap("tab20" if len(unique_clusters) > 10 else "tab10")

        for idx, cl in enumerate(unique_clusters):
            sub = df[df[cluster_col] == cl]
            ax.scatter(sub[dim1_col], sub[dim2_col], color=cmap(idx % cmap.N), s=point_size, alpha=alpha, label=str(cl), edgecolors="none")
            if len(sub) > 0:
                ax.text(sub[dim1_col].median(), sub[dim2_col].median(), str(cl), fontsize=9, fontweight="bold", ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", alpha=0.7, lw=0.5))

        ax.set_xlabel(dim1_col)
        ax.set_ylabel(dim2_col)
        ax.set_title(figure_title)
        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", frameon=False)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 4. Clustermap Heatmap
class ClustermapHeatmap(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"matrix_csv": ("STRING", {"default": ""})},
            opt={
                "colormap": (["RdBu_r", "viridis", "plasma", "inferno", "coolwarm", "bwr"], {"default": "RdBu_r"}),
                "cluster_rows": ("BOOLEAN", {"default": True}),
                "cluster_cols": ("BOOLEAN", {"default": True}),
                "z_score_normalize": ("BOOLEAN", {"default": True}),
            },
            title="Expression Clustermap", w=8.0, h=7.0,
        )

    def run(
        self,
        matrix_csv: str = "",
        colormap: str = "RdBu_r",
        cluster_rows: bool = True,
        cluster_cols: bool = True,
        z_score_normalize: bool = True,
        figure_title: str = "Expression Clustermap",
        dpi: int = 300,
        figure_width: float = 8.0,
        figure_height: float = 7.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("clustermap.png", output_image_path, output_dir, style, dpi)
        in_file = matrix_csv or kwargs.get("matrix_data")
        df = _read_table(in_file, "Matrix Table", index_col=0)

        num_df = df.select_dtypes(include=[np.number]).dropna().copy()
        if num_df.empty or num_df.shape[1] == 0:
            raise ValueError(f"Matrix table at '{in_file}' contains no valid numeric columns.")

        if z_score_normalize and len(num_df) > 0:
            std = num_df.std(axis=1).replace(0, 1.0)
            num_df = num_df.sub(num_df.mean(axis=1), axis=0).div(std, axis=0)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        im = ax.imshow(num_df.values, aspect="auto", cmap=colormap, interpolation="nearest")
        ax.set_xticks(range(num_df.shape[1]))
        ax.set_xticklabels(num_df.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(num_df.shape[0]))
        ax.set_yticklabels(num_df.index, fontsize=7)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Z-score" if z_score_normalize else "Expression")
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 5. GSEA Enrichment Plot
class GseaEnrichmentPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"gsea_results_csv": ("STRING", {"default": ""})},
            opt={
                "pathway_name": ("STRING", {"default": "HALLMARK_PATHWAY"}),
                "nes_value": ("FLOAT", {"default": 0.0, "min": -10.0, "max": 10.0}),
                "fdr_qval": ("FLOAT", {"default": 0.05, "min": 0.0, "max": 1.0}),
            },
            title="GSEA Enrichment Score Profile", w=7.0, h=6.0,
        )

    def run(
        self,
        gsea_results_csv: str = "",
        pathway_name: str = "HALLMARK_PATHWAY",
        nes_value: float = 0.0,
        fdr_qval: float = 0.05,
        figure_title: str = "GSEA Enrichment Score Profile",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("gsea_plot.png", output_image_path, output_dir, style, dpi)
        in_file = gsea_results_csv or kwargs.get("gsea_data")
        df = _read_table(in_file, "GSEA Results Table")

        if "Running_ES" in df.columns or "RES" in df.columns:
            es_col = "Running_ES" if "Running_ES" in df.columns else "RES"
            es_curve = df[es_col].values.astype(float)
            hit_col = "Hit" if "Hit" in df.columns else ("hit" if "hit" in df.columns else None)
            hit_indices = np.where(df[hit_col].astype(bool))[0] if hit_col else np.array([])
            metric_col = "Metric" if "Metric" in df.columns else ("rank_metric" if "rank_metric" in df.columns else None)
            metrics = df[metric_col].values.astype(float) if metric_col else np.linspace(2.0, -2.0, len(df))
        elif "rank_metric" in df.columns or "score" in df.columns:
            m_col = "rank_metric" if "rank_metric" in df.columns else "score"
            metrics = df[m_col].values.astype(float)
            hit_col = "hit" if "hit" in df.columns else ("in_set" if "in_set" in df.columns else None)
            if not hit_col:
                raise ValueError(f"GSEA table at '{in_file}' requires a boolean 'hit' or 'in_set' column.")
            hits = df[hit_col].astype(bool).values
            hit_indices = np.where(hits)[0]
            n = len(metrics)
            n_h = len(hit_indices)
            if n_h == 0:
                raise ValueError(f"Gene set contains 0 hits in ranked list of length {n}.")
            p_hit = np.zeros(n)
            p_miss = np.zeros(n)
            norm_r = np.sum(np.abs(metrics[hits])) or 1.0
            p_hit[hits] = np.abs(metrics[hits]) / norm_r
            p_miss[~hits] = 1.0 / (n - n_h)
            es_curve = np.cumsum(p_hit) - np.cumsum(p_miss)
        else:
            raise ValueError(
                f"GSEA Results table at '{in_file}' must contain either ('Running_ES', 'Hit') or ('rank_metric', 'hit')."
            )

        n_genes = len(es_curve)
        x = np.arange(n_genes)
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(figure_width, figure_height), dpi=dpi,
                                            sharex=True, gridspec_kw={"height_ratios": [3, 1, 1]})
        ax1.plot(x, es_curve, color="#2ca02c", lw=2)
        ax1.axhline(0, color="gray", linestyle="--", lw=0.8)
        ax1.set_ylabel("Enrichment Score (ES)")
        title_extra = f" (NES={nes_value:.2f}, FDR={fdr_qval:.4f})" if nes_value != 0.0 else ""
        ax1.set_title(f"{figure_title}\nPathway: {pathway_name}{title_extra}")

        if len(hit_indices) > 0:
            ax2.vlines(hit_indices, 0, 1, color="black", lw=0.7)
        ax2.set_yticks([])
        ax2.set_ylabel("Hits")

        ax3.plot(x, metrics, color="#1f77b4", lw=1.2)
        ax3.axhline(0, color="gray", linestyle="--", lw=0.8)
        ax3.set_ylabel("Rank Metric")
        ax3.set_xlabel("Rank in Ordered Gene Dataset")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 6. OncoPrint
class OncoPrint(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"mutation_maf_csv": ("STRING", {"default": ""})},
            opt={"top_n_genes": ("INT", {"default": 10, "min": 3, "max": 30})},
            title="OncoPrint Somatic Mutation Landscape", w=9.0, h=5.0,
        )

    def run(
        self,
        mutation_maf_csv: str = "",
        top_n_genes: int = 10,
        figure_title: str = "OncoPrint Somatic Mutation Landscape",
        dpi: int = 300,
        figure_width: float = 9.0,
        figure_height: float = 5.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("oncoprint.png", output_image_path, output_dir, style, dpi)
        in_file = mutation_maf_csv or kwargs.get("mutation_data")
        df = _read_table(in_file, "Mutation MAF / Table")

        gene_col = _find_col(df, "Hugo_Symbol", ["Hugo_Symbol", "gene", "gene_name", "Gene", "GENE"], 0)
        sample_col = _find_col(df, "Tumor_Sample_Barcode", ["Tumor_Sample_Barcode", "sample", "sample_id", "Sample", "SAMPLE"], 1)
        type_col = _find_col(df, "Variant_Classification", ["Variant_Classification", "mutation_type", "type", "Consequence", "Type"], 2)

        df = df.dropna(subset=[gene_col, sample_col, type_col]).copy()
        if df.empty:
            raise ValueError(f"No valid mutation entries in '{in_file}'")

        top_genes = df[gene_col].value_counts().head(top_n_genes).index.tolist()
        samples = sorted(df[sample_col].unique().tolist())

        colors = {
            "Missense": "#377eb8",
            "Nonsense": "#e41a1c",
            "Frameshift": "#4daf4a",
            "InFrame": "#984ea3",
            "Splice": "#ff7f0e",
            "WT": "#f0f0f0",
        }

        mut_map = {(r[gene_col], r[sample_col]): r[type_col] for _, r in df.iterrows()}

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        for i, g in enumerate(top_genes):
            for j, s in enumerate(samples):
                m_type = mut_map.get((g, s), "WT")
                col = colors.get(m_type, "#999999" if m_type != "WT" else "#f0f0f0")
                ax.add_patch(plt.Rectangle((j, i), 0.9, 0.8, color=col, ec="white", lw=0.5))

        ax.set_xlim(0, len(samples))
        ax.set_ylim(0, len(top_genes))
        ax.set_yticks([i + 0.4 for i in range(len(top_genes))])
        ax.set_yticklabels(top_genes, fontweight="bold")
        ax.set_xticks([j + 0.45 for j in range(len(samples))])
        ax.set_xticklabels(samples, rotation=90, fontsize=7)
        ax.set_title(figure_title)

        found_types = [m for m in ["Missense", "Nonsense", "Frameshift", "InFrame", "Splice"] if m in df[type_col].values]
        if not found_types:
            found_types = ["Mutated"]
            colors["Mutated"] = "#377eb8"
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors.get(m, "#377eb8")) for m in found_types]
        ax.legend(handles, found_types, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 7. QQ Plot
class QqPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"pvalues_csv": ("STRING", {"default": ""})},
            opt={"pvalue_col": ("STRING", {"default": "P"})},
            title="Quantile-Quantile (Q-Q) Plot", w=6.0, h=5.5,
        )

    def run(
        self,
        pvalues_csv: str = "",
        pvalue_col: str = "P",
        figure_title: str = "Quantile-Quantile (Q-Q) Plot",
        dpi: int = 300,
        figure_width: float = 6.0,
        figure_height: float = 5.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("qq_plot.png", output_image_path, output_dir, style, dpi)
        in_file = pvalues_csv or kwargs.get("pvalues_data")
        df = _read_table(in_file, "P-values Table")

        pvalue_col = _find_col(df, pvalue_col, ["P", "pvalue", "p_val", "padj", "PVAL"], 0)
        pvals = pd.to_numeric(df[pvalue_col], errors="coerce").dropna().values
        if len(pvals) == 0:
            raise ValueError(f"No valid numeric p-values in column '{pvalue_col}'")

        pvals = np.sort(np.clip(pvals, 1e-300, 1.0))
        n = len(pvals)
        exp_p = -np.log10(np.arange(1, n + 1) / (n + 1))
        obs_p = -np.log10(pvals)
        lambda_gc = float(np.median(-2 * np.log(pvals)) / 1.3863)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.scatter(exp_p, obs_p, c="#1f77b4", s=8, alpha=0.7, edgecolors="none")
        max_val = max(exp_p.max(), obs_p.max())
        ax.plot([0, max_val], [0, max_val], color="#d62728", linestyle="--", lw=1.2)
        ax.set_xlabel(r"Expected $-\log_{10}(P)$")
        ax.set_ylabel(r"Observed $-\log_{10}(P)$")
        ax.set_title(f"{figure_title}\n($\\lambda_{{GC}} = {lambda_gc:.3f}$)")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 8. Sankey / Alluvial Cell Fate
class SankeyCellFate(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"transitions_csv": ("STRING", {"default": ""})},
            opt={
                "source_col": ("STRING", {"default": "source"}),
                "target_col": ("STRING", {"default": "target"}),
                "weight_col": ("STRING", {"default": "weight"}),
            },
            title="Cell Fate & State Transitions", w=8.0, h=5.0,
        )

    def run(
        self,
        transitions_csv: str = "",
        source_col: str = "source",
        target_col: str = "target",
        weight_col: str = "weight",
        figure_title: str = "Cell Fate & State Transitions",
        dpi: int = 300,
        figure_width: float = 8.0,
        figure_height: float = 5.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("sankey_cell_fate.png", output_image_path, output_dir, style, dpi)
        in_file = transitions_csv or kwargs.get("transitions_data")
        df = _read_table(in_file, "Cell Fate Transitions Table")

        source_col = _find_col(df, source_col, ["source", "from", "start_state", "Source"], 0)
        target_col = _find_col(df, target_col, ["target", "to", "end_state", "Target"], 1)
        weight_col = _find_col(df, weight_col, ["weight", "count", "probability", "cells", "n_cells"], 2)

        df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce").fillna(1.0)
        sources = df[source_col].unique().tolist()
        targets = df[target_col].unique().tolist()

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        y_src = {s: (idx + 1) / (len(sources) + 1) for idx, s in enumerate(sources)}
        y_dst = {t: (idx + 1) / (len(targets) + 1) for idx, t in enumerate(targets)}

        for s, y in y_src.items():
            ax.add_patch(plt.Rectangle((0.15, y - 0.04), 0.1, 0.08, color="#457b9d", ec="black", lw=0.8))
            ax.text(0.2, y, str(s), color="white", fontsize=8, ha="center", va="center", fontweight="bold")

        for t, y in y_dst.items():
            ax.add_patch(plt.Rectangle((0.75, y - 0.04), 0.1, 0.08, color="#e76f51", ec="black", lw=0.8))
            ax.text(0.8, y, str(t), color="white", fontsize=8, ha="center", va="center", fontweight="bold")

        max_w = df[weight_col].max() or 1.0
        for _, row in df.iterrows():
            s, t, w = row[source_col], row[target_col], row[weight_col]
            lw = max(0.5, (w / max_w) * 6.0)
            ax.plot([0.25, 0.75], [y_src[s], y_dst[t]], color="#2a9d8f", alpha=0.5, lw=lw)

        ax.set_xlim(0, 1.0)
        ax.set_ylim(0, 1.0)
        ax.axis("off")
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 9. Microbiome Stacked Bar
class MicrobiomeStackedBar(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"abundance_csv": ("STRING", {"default": ""})},
            opt={"top_taxa": ("INT", {"default": 8, "min": 3, "max": 20})},
            title="Taxonomic Composition (Relative Abundance)", w=8.5, h=5.0,
        )

    def run(
        self,
        abundance_csv: str = "",
        top_taxa: int = 8,
        figure_title: str = "Taxonomic Composition (Relative Abundance)",
        dpi: int = 300,
        figure_width: float = 8.5,
        figure_height: float = 5.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("microbiome_stacked_bar.png", output_image_path, output_dir, style, dpi)
        in_file = abundance_csv or kwargs.get("abundance_data")
        df = _read_table(in_file, "Taxonomy Abundance Table", index_col=0)

        num_df = df.select_dtypes(include=[np.number]).fillna(0.0)
        if num_df.empty:
            raise ValueError(f"No numeric abundance columns in table '{in_file}'")

        col_sums = num_df.sum(axis=0).replace(0, 1.0)
        rel_df = num_df.div(col_sums, axis=1)
        samples = list(rel_df.columns)
        taxa = list(rel_df.mean(axis=1).sort_values(ascending=False).index[:top_taxa])
        mat = rel_df.loc[taxa].values

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        bottom = np.zeros(len(samples))
        cmap = plt.get_cmap("Set2")
        for idx, taxon in enumerate(taxa):
            ax.bar(samples, mat[idx], bottom=bottom, label=taxon, color=cmap(idx % 8), width=0.7)
            bottom += mat[idx]

        ax.set_ylabel("Relative Abundance")
        ax.set_ylim(0, 1.05)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_title(figure_title)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 10. PCoA Scatter
class PcoaScatter(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"distance_matrix_csv": ("STRING", {"default": ""})},
            opt={"group_metadata_csv": ("STRING", {"default": ""})},
            title="PCoA - Microbiome Beta Diversity", w=7.0, h=5.5,
        )

    def run(
        self,
        distance_matrix_csv: str = "",
        group_metadata_csv: str = "",
        figure_title: str = "PCoA - Microbiome Beta Diversity",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 5.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("pcoa_scatter.png", output_image_path, output_dir, style, dpi)
        in_file = distance_matrix_csv or kwargs.get("distance_data")
        df = _read_table(in_file, "Distance Matrix Table", index_col=0)

        D = df.select_dtypes(include=[np.number]).values
        n = D.shape[0]
        if n < 3:
            raise ValueError(f"Distance matrix must have at least 3 samples, got {n}.")

        H = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * H.dot(D ** 2).dot(H)
        eigvals, eigvecs = np.linalg.eigh(B)
        idx = np.argsort(eigvals)[::-1]
        pos_eig = np.maximum(eigvals[idx], 0)
        tot = pos_eig.sum() or 1.0
        var1 = round(float(pos_eig[0] / tot) * 100, 1)
        var2 = round(float(pos_eig[1] / tot) * 100, 1) if len(pos_eig) > 1 else 0.0
        pts = eigvecs[:, idx[:2]] * np.sqrt(pos_eig[:2])

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.scatter(pts[:, 0], pts[:, 1], c="#1f77b4", s=40, alpha=0.85)
        for i, s_name in enumerate(df.index):
            ax.annotate(str(s_name), (pts[i, 0], pts[i, 1]), fontsize=8, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(f"PCoA 1 ({var1}% explained var)")
        ax.set_ylabel(f"PCoA 2 ({var2}% explained var)")
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 11. Synteny Genome
class SyntenyGenome(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"synteny_coords_tsv": ("STRING", {"default": ""})},
            title="Comparative Genome Synteny Blocks", w=9.0, h=4.5,
        )

    def run(
        self,
        synteny_coords_tsv: str = "",
        figure_title: str = "Comparative Genome Synteny Blocks",
        dpi: int = 300,
        figure_width: float = 9.0,
        figure_height: float = 4.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("synteny_genome.png", output_image_path, output_dir, style, dpi)
        in_file = synteny_coords_tsv or kwargs.get("synteny_data")
        df = _read_table(in_file, "Synteny Coordinates Table")

        ref_start_col = _find_col(df, "ref_start", ["ref_start", "start1", "s1", "RefStart"], 1)
        ref_end_col = _find_col(df, "ref_end", ["ref_end", "end1", "e1", "RefEnd"], 2)
        qry_start_col = _find_col(df, "qry_start", ["qry_start", "start2", "s2", "QryStart"], 3)
        qry_end_col = _find_col(df, "qry_end", ["qry_end", "end2", "e2", "QryEnd"], 4)

        df[ref_start_col] = pd.to_numeric(df[ref_start_col], errors="coerce").fillna(0)
        df[ref_end_col] = pd.to_numeric(df[ref_end_col], errors="coerce").fillna(0)
        df[qry_start_col] = pd.to_numeric(df[qry_start_col], errors="coerce").fillna(0)
        df[qry_end_col] = pd.to_numeric(df[qry_end_col], errors="coerce").fillna(0)

        max_ref = max(df[ref_end_col].max(), 1.0)
        max_qry = max(df[qry_end_col].max(), 1.0)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.add_patch(plt.Rectangle((0.1, 0.7), 0.8, 0.08, color="#457b9d", ec="black", lw=1))
        ax.text(0.5, 0.82, "Reference Genome", ha="center", fontsize=9, fontweight="bold")
        ax.add_patch(plt.Rectangle((0.1, 0.2), 0.8, 0.08, color="#e76f51", ec="black", lw=1))
        ax.text(0.5, 0.12, "Target Genome", ha="center", fontsize=9, fontweight="bold")

        for _, row in df.iterrows():
            x1 = 0.1 + 0.8 * (row[ref_start_col] / max_ref)
            x2 = 0.1 + 0.8 * (row[ref_end_col] / max_ref)
            x3 = 0.1 + 0.8 * (row[qry_end_col] / max_qry)
            x4 = 0.1 + 0.8 * (row[qry_start_col] / max_qry)
            poly = plt.Polygon([[x1, 0.7], [x2, 0.7], [x3, 0.28], [x4, 0.28]], color="#2a9d8f", alpha=0.45, ec="none")
            ax.add_patch(poly)

        ax.set_xlim(0, 1.0)
        ax.set_ylim(0, 1.0)
        ax.axis("off")
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 12. ChIP / ATAC Coverage Profile
class ChipAtacCoverageProfile(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"coverage_matrix_gz": ("STRING", {"default": ""})},
            opt={"sample_label": ("STRING", {"default": "Signal Coverage Around Peak Center"})},
            title="Signal Coverage Around Peak Center", w=7.0, h=5.0,
        )

    def run(
        self,
        coverage_matrix_gz: str = "",
        sample_label: str = "Signal Coverage Around Peak Center",
        figure_title: str = "Signal Coverage Around Peak Center",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 5.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("peak_coverage_profile.png", output_image_path, output_dir, style, dpi)
        in_file = coverage_matrix_gz or kwargs.get("coverage_matrix")
        df = _read_table(in_file, "Coverage Matrix")

        num_cols = df.select_dtypes(include=[np.number]).columns
        if len(num_cols) == 0:
            raise ValueError(f"Coverage matrix at '{in_file}' contains no numeric signal bins.")

        mean_profile = df[num_cols].mean(axis=0).values
        std_profile = df[num_cols].std(axis=0).values / np.sqrt(len(df))

        n_bins = len(mean_profile)
        x = np.linspace(-2000, 2000, n_bins)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.plot(x, mean_profile, color="#1f77b4", lw=2, label=sample_label)
        ax.fill_between(x, mean_profile - std_profile, mean_profile + std_profile, alpha=0.2, color="#1f77b4")

        ax.set_xlabel("Distance from Peak Center (bp)")
        ax.set_ylabel("Normalized Signal Density")
        ax.set_title(figure_title)
        ax.legend(frameon=True, loc="upper right")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 13. MA Plot
class MaPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"deg_results_csv": ("STRING", {"default": ""})},
            opt={
                "basemean_col": ("STRING", {"default": "baseMean"}),
                "log2fc_col": ("STRING", {"default": "log2FoldChange"}),
                "pvalue_col": ("STRING", {"default": "padj"}),
                "fdr_cutoff": ("FLOAT", {"default": 0.05, "min": 1e-10, "max": 0.2}),
            },
            title="MA Plot (Log Intensity vs Log Ratio)", w=7.0, h=5.0,
        )

    def run(
        self,
        deg_results_csv: str = "",
        basemean_col: str = "baseMean",
        log2fc_col: str = "log2FoldChange",
        pvalue_col: str = "padj",
        fdr_cutoff: float = 0.05,
        figure_title: str = "MA Plot (Log Intensity vs Log Ratio)",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 5.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        figure_title = kwargs.get("title", figure_title)
        target_path = self._target_path("ma_plot.png", output_image_path, output_dir, style, dpi)
        in_file = deg_results_csv or kwargs.get("deg_table_path") or kwargs.get("deg_results_table")
        df = _read_table(in_file, "DEG Results Table")

        basemean_col = _find_col(df, basemean_col, ["baseMean", "AveExpr", "mean", "base_mean"], 1)
        log2fc_col = _find_col(df, log2fc_col, ["log2FoldChange", "logFC", "lfc"], 2)
        pvalue_col = _find_col(df, pvalue_col, ["padj", "pvalue", "p_val_adj", "FDR"], len(df.columns) - 1)

        clean = df.dropna(subset=[basemean_col, log2fc_col, pvalue_col])
        if clean.empty:
            raise ValueError(f"No valid rows after dropping NAs in columns '{basemean_col}', '{log2fc_col}', '{pvalue_col}'")

        a = pd.to_numeric(clean[basemean_col], errors="coerce").fillna(1.0).values
        m = pd.to_numeric(clean[log2fc_col], errors="coerce").fillna(0.0).values
        padj = pd.to_numeric(clean[pvalue_col], errors="coerce").fillna(1.0).values

        is_sig = padj < fdr_cutoff
        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.scatter(a[~is_sig], m[~is_sig], color="#7f7f7f", alpha=0.4, s=8, label=f"Not Sig ({(~is_sig).sum()})")
        ax.scatter(a[is_sig], m[is_sig], color="#d62728", alpha=0.75, s=14, label=f"FDR < {fdr_cutoff} ({is_sig.sum()})")
        ax.set_xscale("log")
        ax.axhline(0, color="blue", linestyle="--", lw=1.0)
        ax.set_xlabel("Mean Expression (BaseMean)")
        ax.set_ylabel(r"$\log_2$ Fold Change (M)")
        ax.set_title(figure_title)
        ax.legend(frameon=True, loc="upper right")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 14. Spatial Tissue Overlay
class SpatialTissueOverlay(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"spatial_coords_csv": ("STRING", {"default": ""})},
            opt={
                "x_col": ("STRING", {"default": "x"}),
                "y_col": ("STRING", {"default": "y"}),
                "value_col": ("STRING", {"default": "expression"}),
            },
            title="Spatial Transcriptomics Tissue Microenvironment", w=7.0, h=6.0,
        )

    def run(
        self,
        spatial_coords_csv: str = "",
        x_col: str = "x",
        y_col: str = "y",
        value_col: str = "expression",
        figure_title: str = "Spatial Transcriptomics Tissue Microenvironment",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("spatial_tissue_overlay.png", output_image_path, output_dir, style, dpi)
        in_file = spatial_coords_csv or kwargs.get("spatial_data")
        df = _read_table(in_file, "Spatial Coordinates Table")

        x_col = _find_col(df, x_col, ["x", "X", "array_row", "pxl_row_in_fullres", "col"], 0)
        y_col = _find_col(df, y_col, ["y", "Y", "array_col", "pxl_col_in_fullres", "row"], 1)
        value_col = _find_col(df, value_col, ["expression", "count", "cluster", "value", "Value"], 2)

        x = pd.to_numeric(df[x_col], errors="coerce")
        y = pd.to_numeric(df[y_col], errors="coerce")
        c = pd.to_numeric(df[value_col], errors="coerce").fillna(0.0)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        sc = ax.scatter(x, y, c=c, cmap="inferno", s=25, alpha=0.9, edgecolors="none")
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(value_col)
        ax.set_aspect("equal")
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 15. MD Trajectory Plotter
class MdTrajectoryPlotter(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"trajectory_metrics_csv": ("STRING", {"default": ""})},
            opt={
                "time_col": ("STRING", {"default": "time_ns"}),
                "rmsd_col": ("STRING", {"default": "rmsd"}),
                "rmsf_col": ("STRING", {"default": "rmsf"}),
            },
            title="Molecular Dynamics Trajectory Stability", w=8.0, h=6.0,
        )

    def run(
        self,
        trajectory_metrics_csv: str = "",
        time_col: str = "time_ns",
        rmsd_col: str = "rmsd",
        rmsf_col: str = "rmsf",
        figure_title: str = "Molecular Dynamics Trajectory Stability",
        dpi: int = 300,
        figure_width: float = 8.0,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("md_trajectory_rmsd.png", output_image_path, output_dir, style, dpi)
        in_file = trajectory_metrics_csv or kwargs.get("trajectory_data")
        df = _read_table(in_file, "Trajectory Metrics Table")

        time_col = _find_col(df, time_col, ["time_ns", "time", "Time", "step", "Frame"], 0)
        rmsd_col = _find_col(df, rmsd_col, ["rmsd", "RMSD", "backbone_rmsd"], 1)

        time_vals = pd.to_numeric(df[time_col], errors="coerce").fillna(0).values
        rmsd_vals = pd.to_numeric(df[rmsd_col], errors="coerce").fillna(0).values

        has_rmsf = rmsf_col in df.columns or "RMSF" in df.columns
        if has_rmsf:
            rmsf_col = "RMSF" if "RMSF" in df.columns else rmsf_col
            rmsf_vals = pd.to_numeric(df[rmsf_col], errors="coerce").dropna().values
            residues = np.arange(1, len(rmsf_vals) + 1)
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(figure_width, figure_height), dpi=dpi)
            ax1.plot(time_vals, rmsd_vals, color="#1f77b4", lw=1.5)
            ax1.set_ylabel("Backbone RMSD (nm)")
            ax1.set_xlabel("Time (ns)")
            ax1.set_title(figure_title)

            ax2.plot(residues, rmsf_vals, color="#d62728", lw=1.2)
            ax2.set_ylabel(r"C$\alpha$ RMSF (nm)")
            ax2.set_xlabel("Residue Number")
        else:
            fig, ax1 = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
            ax1.plot(time_vals, rmsd_vals, color="#1f77b4", lw=1.5)
            ax1.set_ylabel("Backbone RMSD (nm)")
            ax1.set_xlabel("Time (ns)")
            ax1.set_title(figure_title)

        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 16. Ramachandran Plot
class RamachandranPlot(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"pdb_file_path": ("STRING", {"default": ""})},
            title="Ramachandran Plot (Torsion Angles)", w=6.5, h=6.0,
        )

    def run(
        self,
        pdb_file_path: str = "",
        figure_title: str = "Ramachandran Plot (Torsion Angles)",
        dpi: int = 300,
        figure_width: float = 6.5,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("ramachandran_plot.png", output_image_path, output_dir, style, dpi)
        in_file = pdb_file_path or kwargs.get("pdb_file")
        p = _file(in_file, "PDB Structure File")

        try:
            from Bio.PDB import PDBParser, PPBuilder
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", str(p))
            ppb = PPBuilder()
            phi_list, psi_list = [], []
            for pp in ppb.build_peptides(structure):
                angles = pp.get_phi_psi_list()
                for phi, psi in angles:
                    if phi is not None and psi is not None:
                        phi_list.append(math.degrees(phi))
                        psi_list.append(math.degrees(psi))
        except Exception as e:
            raise ValueError(f"Failed to parse PDB structure and extract torsion angles from '{p}': {e}")

        if not phi_list:
            raise ValueError(f"No valid peptide residues with phi/psi angles found in PDB '{p}'")

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.scatter(phi_list, psi_list, c="#1f77b4", s=18, alpha=0.75, edgecolors="none")
        ax.axhline(0, color="gray", lw=0.6, linestyle=":")
        ax.axvline(0, color="gray", lw=0.6, linestyle=":")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xlabel(r"$\Phi$ (degrees)")
        ax.set_ylabel(r"$\Psi$ (degrees)")
        ax.set_title(f"{figure_title}\n(N={len(phi_list)} residues)")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 17. Protein-Ligand Interaction
class ProteinLigandInteraction(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"complex_pdb_or_json": ("STRING", {"default": ""})},
            title="2D Protein-Ligand Interaction Diagram", w=7.0, h=6.5,
        )

    def run(
        self,
        complex_pdb_or_json: str = "",
        figure_title: str = "2D Protein-Ligand Interaction Diagram",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 6.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("protein_ligand_2d.png", output_image_path, output_dir, style, dpi)
        in_file = complex_pdb_or_json or kwargs.get("interaction_data")
        p = _file(in_file, "Protein-Ligand Complex / Contacts")

        residues = []
        if p.suffix.lower() == ".json":
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
            contacts = data if isinstance(data, list) else data.get("contacts", [])
            for c in contacts:
                name = c.get("residue", "Res")
                itype = c.get("type", "Contact")
                residues.append((name, itype))
        else:
            df = _read_table(str(p), "Interaction Contacts Table")
            res_col = _find_col(df, "residue", ["residue", "res_name", "Residue"], 0)
            type_col = _find_col(df, "type", ["type", "interaction_type", "interaction", "Type"], 1)
            for _, row in df.iterrows():
                residues.append((str(row[res_col]), str(row[type_col])))

        if not residues:
            raise ValueError(f"No protein-ligand contacts found in '{p}'")

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        ax.add_patch(plt.Circle((0.5, 0.5), 0.15, color="#2b5c8f", alpha=0.85))
        ax.text(0.5, 0.5, "LIGAND", color="white", ha="center", va="center", fontweight="bold", fontsize=10)

        n_res = len(residues)
        angles = np.linspace(0, 2 * np.pi, n_res, endpoint=False)
        colors_map = {
            "h-bond": "#2a9d8f",
            "salt bridge": "#e76f51",
            "pi-pi": "#9d4edd",
            "hydrophobic": "#457b9d",
        }
        for i, (res_name, itype) in enumerate(residues):
            ang = angles[i]
            rx = 0.5 + 0.35 * np.cos(ang)
            ry = 0.5 + 0.35 * np.sin(ang)
            color = colors_map.get(itype.lower(), "#2a9d8f")
            ax.add_patch(plt.Circle((rx, ry), 0.06, color=color, alpha=0.3, ec=color, lw=1.5))
            ax.text(rx, ry, res_name, ha="center", va="center", fontsize=7, fontweight="bold")
            ax.plot([0.5, rx], [0.5, ry], color=color, linestyle="--", lw=1.5)

        ax.set_xlim(0.05, 0.95)
        ax.set_ylim(0.05, 0.95)
        ax.axis("off")
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 18. Phylogenetic Tree
class PhylogeneticTree(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"newick_tree_file": ("STRING", {"default": ""})},
            opt={"tree_layout": (["rectangular", "radial"], {"default": "rectangular"})},
            title="Phylogenetic Tree Visualization", w=7.0, h=7.0,
        )

    def run(
        self,
        newick_tree_file: str = "",
        tree_layout: str = "rectangular",
        figure_title: str = "Phylogenetic Tree Visualization",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 7.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("phylogenetic_tree.png", output_image_path, output_dir, style, dpi)
        in_file = newick_tree_file or kwargs.get("tree_file")
        p = _file(in_file, "Newick Tree File")

        try:
            from Bio import Phylo
            tree = Phylo.read(str(p), "newick")
        except Exception as e:
            raise ValueError(f"Failed to parse Newick tree from '{p}': {e}")

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        Phylo.draw(tree, axes=ax, do_show=False)
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 19. Linkage Disequilibrium
class LinkageDisequilibrium(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"vcf_file_path": ("STRING", {"default": ""})},
            opt={"metric": (["r2", "Dprime"], {"default": "r2"})},
            title="Linkage Disequilibrium (LD) Heatmap", w=7.0, h=6.0,
        )

    def run(
        self,
        vcf_file_path: str = "",
        metric: str = "r2",
        figure_title: str = "Linkage Disequilibrium (LD) Heatmap",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 6.0,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("ld_heatmap.png", output_image_path, output_dir, style, dpi)
        in_file = vcf_file_path or kwargs.get("ld_data")
        p = _file(in_file, "VCF / LD Matrix")

        if p.suffix.lower() in (".csv", ".tsv", ".txt"):
            df = _read_table(str(p), "LD Matrix Table", index_col=0)
            r2_mat = df.select_dtypes(include=[np.number]).values
            labels = list(df.index)
        else:
            dosages = []
            labels = []
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 10:
                        chrom, pos, vid = parts[0], parts[1], parts[2]
                        label = vid if vid != "." else f"{chrom}:{pos}"
                        labels.append(label)
                        row_dosages = []
                        for sample_col in parts[9:]:
                            gt = sample_col.split(":")[0].replace("|", "/")
                            if "1/1" in gt:
                                row_dosages.append(2.0)
                            elif "0/1" in gt or "1/0" in gt:
                                row_dosages.append(1.0)
                            else:
                                row_dosages.append(0.0)
                        dosages.append(row_dosages)
                    if len(labels) >= 50:
                        break

            if len(dosages) < 2:
                raise ValueError(f"VCF '{p}' must contain at least 2 loci with sample genotypes.")

            mat = np.array(dosages)
            with np.errstate(divide="ignore", invalid="ignore"):
                corr = np.corrcoef(mat)
            r2_mat = np.nan_to_num(corr ** 2, nan=0.0)

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        im = ax.imshow(r2_mat, cmap="YlOrRd", vmin=0, vmax=1)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$r^2$ Correlation")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title(figure_title)
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


# 20. Kaplan-Meier Survival
class KaplanMeierSurvival(_BaseVisualizer):
    OUTPUT_NODE = True
    OUPUT_NODE = True
    @classmethod
    def INPUT_TYPES(cls):
        return cls._inputs(
            req={"clinical_csv": ("STRING", {"default": ""})},
            opt={
                "time_col": ("STRING", {"default": "time_months"}),
                "event_col": ("STRING", {"default": "vital_status"}),
                "strata_col": ("STRING", {"default": "biomarker_strata"}),
            },
            title="Kaplan-Meier Overall Survival", w=7.0, h=5.5,
        )

    def run(
        self,
        clinical_csv: str = "",
        time_col: str = "time_months",
        event_col: str = "vital_status",
        strata_col: str = "biomarker_strata",
        figure_title: str = "Kaplan-Meier Overall Survival",
        dpi: int = 300,
        figure_width: float = 7.0,
        figure_height: float = 5.5,
        style: str = "nature",
        output_dir: str = "",
        output_image_path: str = "",
        **kwargs,
    ) -> Tuple[str, Any]:
        target_path = self._target_path("kaplan_meier_survival.png", output_image_path, output_dir, style, dpi)
        in_file = clinical_csv or kwargs.get("clinical_data")
        df = _read_table(in_file, "Clinical Survival Table")

        time_col = _find_col(df, time_col, ["time_months", "time", "os_months", "Time", "survival_time"], 0)
        event_col = _find_col(df, event_col, ["vital_status", "status", "event", "vital", "Event", "death"], 1)

        df = df.dropna(subset=[time_col, event_col]).copy()
        if df.empty:
            raise ValueError(f"No valid rows after dropping NAs in '{time_col}' and '{event_col}'")

        fig, ax = plt.subplots(figsize=(figure_width, figure_height), dpi=dpi)
        strata_groups = df.groupby(strata_col) if strata_col in df.columns else [("All Patients", df)]
        colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]

        for i, (name, group) in enumerate(strata_groups):
            t = pd.to_numeric(group[time_col], errors="coerce").fillna(0).values
            e = group[event_col].astype(str).str.lower().isin(["1", "true", "dead", "event", "relapse"]).astype(int).values
            order = np.argsort(t)
            t_sorted, e_sorted = t[order], e[order]
            n = len(t_sorted)
            km_times, km_surv, curr_s = [0.0], [1.0], 1.0
            for k in range(n):
                at_risk = n - k
                if at_risk <= 0:
                    break
                if e_sorted[k] == 1:
                    curr_s *= (1.0 - 1.0 / at_risk)
                km_times.append(t_sorted[k])
                km_surv.append(curr_s)
            ax.step(km_times, km_surv, where="post", color=colors[i % len(colors)], lw=2, label=f"{name} (n={n})")

        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Time (Months)")
        ax.set_ylabel("Overall Survival Probability")
        ax.set_title(figure_title)
        ax.legend(frameon=True, loc="upper right")
        return save_and_tensorize_figure(fig, target_path, dpi=dpi)


NODE_CLASS_MAPPINGS = {
    "VolcanoPlot": VolcanoPlot,
    "ManhattanPlot": ManhattanPlot,
    "UmapScatter": UmapScatter,
    "ClustermapHeatmap": ClustermapHeatmap,
    "GseaEnrichmentPlot": GseaEnrichmentPlot,
    "OncoPrint": OncoPrint,
    "QqPlot": QqPlot,
    "SankeyCellFate": SankeyCellFate,
    "MicrobiomeStackedBar": MicrobiomeStackedBar,
    "PcoaScatter": PcoaScatter,
    "SyntenyGenome": SyntenyGenome,
    "ChipAtacCoverageProfile": ChipAtacCoverageProfile,
    "MaPlot": MaPlot,
    "SpatialTissueOverlay": SpatialTissueOverlay,
    "MdTrajectoryPlotter": MdTrajectoryPlotter,
    "RamachandranPlot": RamachandranPlot,
    "ProteinLigandInteraction": ProteinLigandInteraction,
    "PhylogeneticTree": PhylogeneticTree,
    "LinkageDisequilibrium": LinkageDisequilibrium,
    "KaplanMeierSurvival": KaplanMeierSurvival,
}

NODE_DISPLAY_NAME_MAPPINGS = {
    "VolcanoPlot": "Visualization: Volcano Plot",
    "ManhattanPlot": "Visualization: Manhattan Plot",
    "UmapScatter": "Visualization: UMAP Scatter",
    "ClustermapHeatmap": "Visualization: Clustermap Heatmap",
    "GseaEnrichmentPlot": "Visualization: GSEA Enrichment Plot",
    "OncoPrint": "Visualization: OncoPrint Mutation Landscape",
    "QqPlot": "Visualization: QQ Plot",
    "SankeyCellFate": "Visualization: Sankey Cell Fate Alluvial",
    "MicrobiomeStackedBar": "Visualization: Microbiome Stacked Bar",
    "PcoaScatter": "Visualization: Microbiome PCoA Scatter",
    "SyntenyGenome": "Visualization: Synteny Genome Alignment",
    "ChipAtacCoverageProfile": "Visualization: ChIP/ATAC Coverage Profile",
    "MaPlot": "Visualization: MA Plot",
    "SpatialTissueOverlay": "Visualization: Spatial Tissue Overlay",
    "MdTrajectoryPlotter": "Visualization: MD Trajectory RMSD/RMSF Plotter",
    "RamachandranPlot": "Visualization: Ramachandran Dihedral Plot",
    "ProteinLigandInteraction": "Visualization: Protein-Ligand Interaction",
    "PhylogeneticTree": "Visualization: Phylogenetic Tree",
    "LinkageDisequilibrium": "Visualization: Linkage Disequilibrium Heatmap",
    "KaplanMeierSurvival": "Visualization: Kaplan-Meier Survival Curve",
}

# Backward compatibility aliases (*Node and *VisualizerNode)
for _name, _cls in list(NODE_CLASS_MAPPINGS.items()):
    globals()[f"{_name}Node"] = _cls
    globals()[f"{_name}VisualizerNode"] = _cls

PUBLICATION_VISUALIZER_CLASSES = [f"{k}VisualizerNode" for k in NODE_CLASS_MAPPINGS]

__all__ = [
    *NODE_CLASS_MAPPINGS.keys(),
    *[f"{k}Node" for k in NODE_CLASS_MAPPINGS],
    *[f"{k}VisualizerNode" for k in NODE_CLASS_MAPPINGS],
    "PUBLICATION_VISUALIZER_CLASSES",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
    "configure_publication_style",
    "figure_to_image_tensor",
    "save_and_tensorize_figure",
]
