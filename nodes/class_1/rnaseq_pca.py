"""RNA-seq Exploratory Data Analysis & Sample QC node (PCA & Correlation Heatmap).

Python packages: numpy, pandas, matplotlib
External binaries: none
"""

import csv
import json
from pathlib import Path


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path

def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and str(custom_dir).strip():
        out = Path(custom_dir).expanduser().resolve()
    else:
        try:
            folder_paths = __import__("folder_paths")
            base = Path(folder_paths.get_output_directory())
        except Exception:
            base = Path.cwd() / "ComfyUI" / "output"
        out = base / node_name
    return out



class DESeq2SampleQC:
    CATEGORY = "ComfyBIO/Quality Control"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("pca_plot_png", "correlation_heatmap_png", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "counts_tsv": ("STRING", {"default": ""}),
            },
            "optional": {
                "metadata_tsv": ("STRING", {"default": ""}),
            },
        }

    def run(self, counts_tsv: str, output_dir: str = "", metadata_tsv: str = ""):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        counts_path = _file(counts_tsv, "Gene Counts TSV")
        out = _output_dir("DESeq2SampleQC", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        # Parse counts table
        with open(counts_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            sample_names = header[1:]
            genes = []
            matrix_rows = []
            for row in reader:
                if not row or not row[0].strip():
                    continue
                genes.append(row[0])
                matrix_rows.append([float(x) for x in row[1:]])

        data = np.array(matrix_rows)  # shape: (genes, samples)
        n_genes, n_samples = data.shape

        if n_samples < 2:
            raise ValueError(f"PCA and correlation heatmap require at least 2 samples, got {n_samples}")

        # Library size normalization (CPM) + log2(x + 1)
        col_sums = data.sum(axis=0)
        col_sums[col_sums == 0] = 1.0
        cpm = (data / col_sums) * 1e6
        log_cpm = np.log2(cpm + 1.0)

        # Filter out low variance genes (top 500 variable genes)
        gene_vars = np.var(log_cpm, axis=1)
        top_idx = np.argsort(gene_vars)[::-1][:min(500, n_genes)]
        subset = log_cpm[top_idx, :]  # shape: (features, samples)

        # Center data per gene
        centered = subset - subset.mean(axis=1, keepdims=True)

        # SVD / PCA on samples: transpose so samples are rows (n_samples, n_features)
        X = centered.T
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        pcs = U * S  # shape: (n_samples, n_components)

        explained_var = (S ** 2) / (S ** 2).sum()
        pc1_var = round(float(explained_var[0]) * 100, 1) if len(explained_var) > 0 else 0.0
        pc2_var = round(float(explained_var[1]) * 100, 1) if len(explained_var) > 1 else 0.0

        # Plot PCA
        pca_png = out / "sample_pca.png"
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(pcs[:, 0], pcs[:, 1], c="steelblue", s=80, edgecolors="black", zorder=3)
        for i, s_name in enumerate(sample_names):
            ax.annotate(s_name, (pcs[i, 0], pcs[i, 1]), textcoords="offset points", xytext=(5, 5))
        ax.set_xlabel(f"PC1 ({pc1_var}% variance)")
        ax.set_ylabel(f"PC2 ({pc2_var}% variance)")
        ax.set_title("RNA-seq Sample PCA (log2-CPM)")
        ax.grid(True, linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.savefig(pca_png, dpi=200)
        plt.close(fig)

        # Sample correlation heatmap
        # Compute Pearson correlation matrix across samples
        corr = np.corrcoef(subset, rowvar=False)  # shape: (n_samples, n_samples)
        heatmap_png = out / "sample_correlation_heatmap.png"

        fig, ax = plt.subplots(figsize=(max(6, n_samples * 0.8), max(5, n_samples * 0.8)))
        im = ax.imshow(corr, cmap="YlOrRd", vmin=max(0.5, float(corr.min())), vmax=1.0)
        ax.set_xticks(range(n_samples))
        ax.set_yticks(range(n_samples))
        ax.set_xticklabels(sample_names, rotation=45, ha="right")
        ax.set_yticklabels(sample_names)
        fig.colorbar(im, ax=ax, label="Pearson Correlation")
        ax.set_title("Sample Correlation Heatmap")
        plt.tight_layout()
        plt.savefig(heatmap_png, dpi=200)
        plt.close(fig)

        summary = {
            "num_samples": n_samples,
            "samples": sample_names,
            "pc1_explained_variance_pct": pc1_var,
            "pc2_explained_variance_pct": pc2_var,
        }

        return str(pca_png), str(heatmap_png), json.dumps(summary)


NODE_CLASS_MAPPINGS = {"DESeq2SampleQC": DESeq2SampleQC}
NODE_DISPLAY_NAME_MAPPINGS = {"DESeq2SampleQC": "RNA-Seq: Sample PCA & Correlation QC"}


# Backward compatibility aliases
DESeq2SampleQCNode = DESeq2SampleQC

__all__ = ["DESeq2SampleQC",
    "DESeq2SampleQCNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
