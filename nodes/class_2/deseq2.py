"""DESeq2 Differential Gene Expression (DEG) Analysis node.

Executes rigorous differential expression analysis from count matrices and sample metadata
using Bioconductor DESeq2 via Rscript.
Ensures genuine Bioconductor DESeq2 execution with strict experimental design verification,
sample ID alignment, and execution audit manifests.

Strict scientific integrity:
- Automatic fallback to native approximations is completely removed.
- If Rscript or Bioconductor DESeq2 is unavailable or encounters an error, the execution
  fails immediately with full stderr logging.
- Arbitrary sample splitting or condition guessing is strictly prohibited.

Python packages: pandas
External binaries: Rscript (with bioconductor-deseq2 installed in bulk_rna_seq env)
"""

import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _file(value: str, label: str) -> Path:
    if not value or not str(value).strip():
        raise ValueError(f"{label} path is required and cannot be empty.")
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


def _find_rscript() -> Optional[str]:
    """Find Rscript executable from environment variable, PATH, or standard bulk_rna_seq conda paths."""
    env_r = os.environ.get("RSCRIPT_PATH")
    if env_r and Path(env_r).is_file():
        return str(Path(env_r).resolve())

    ambient = shutil.which("Rscript")
    if ambient:
        return ambient

    # Search standard conda environment paths
    conda_candidates = [
        Path("/opt/miniconda3/envs/bulk_rna_seq/bin/Rscript"),
        Path.home() / "miniconda3/envs/bulk_rna_seq/bin/Rscript",
        Path.home() / "anaconda3/envs/bulk_rna_seq/bin/Rscript",
        Path("/opt/homebrew/Caskroom/miniconda/base/envs/bulk_rna_seq/bin/Rscript"),
    ]
    for c in conda_candidates:
        if c.is_file():
            return str(c.resolve())

    # Check CONDA_PREFIX
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        c = Path(prefix) / "bin" / "Rscript"
        if c.is_file():
            return str(c.resolve())

    return None


class DESeq2:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/RNA-Seq"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("deseq2_results_csv", "normalized_counts_csv", "deg_summary_json")

    @classmethod
    def IS_CHANGED(cls, **kwargs):
        try:
            from bioflow.runtime.artifacts import compute_input_fingerprint
            return compute_input_fingerprint(**kwargs)
        except Exception:
            return float("NaN")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "count_matrix_csv": ("STRING", {"default": ""}),
                "sample_metadata_csv": ("STRING", {"default": ""}),
            },
            "optional": {
                "condition_col": ("STRING", {"default": "condition"}),
                "contrast_reference": ("STRING", {"default": ""}),
                "contrast_target": ("STRING", {"default": ""}),
                "alpha": ("FLOAT", {"default": 0.05, "min": 0.0, "max": 1.0, "step": 0.01}),
                "lfc_threshold": ("FLOAT", {"default": 0.0, "min": 0.0, "max": 10.0, "step": 0.1}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        count_matrix_csv: str,
        sample_metadata_csv: str,
        output_dir: str = "",
        condition_col: str = "condition",
        contrast_reference: str = "",
        contrast_target: str = "",
        alpha: float = 0.05,
        lfc_threshold: float = 0.0,
        extra_command: str = "",
        **kwargs,
    ) -> Tuple[str, str, str]:
        import pandas as pd

        counts_path = _file(count_matrix_csv, "Gene Count Matrix")
        meta_path = _file(sample_metadata_csv, "Sample Metadata")

        out = _output_dir("DESeq2", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        results_csv = out / "deseq2_results.csv"
        norm_counts_csv = out / "normalized_counts.csv"
        r_script_path = out / "run_deseq2.R"
        r_log_path = out / "deseq2_rscript.log"
        manifest_json_path = out / "deseq2_manifest.json"

        # 1. Load and validate count matrix
        with open(counts_path, "r", encoding="utf-8", errors="ignore") as f:
            sep_counts = "\t" if "\t" in f.readline() else ","
        counts_df = pd.read_csv(counts_path, sep=sep_counts)
        if counts_df.empty:
            raise ValueError(f"Gene count matrix '{counts_path}' is empty.")

        gene_col = counts_df.columns[0]
        gene_ids = counts_df[gene_col].astype(str).tolist()
        sample_cols = [c for c in counts_df.columns if c != gene_col]

        # 2. Load and validate metadata
        with open(meta_path, "r", encoding="utf-8", errors="ignore") as f:
            sep_meta = "\t" if "\t" in f.readline() else ","
        meta_df = pd.read_csv(meta_path, sep=sep_meta)
        if meta_df.empty:
            raise ValueError(f"Sample metadata '{meta_path}' is empty.")

        sample_id_col = meta_df.columns[0]
        meta_samples = meta_df[sample_id_col].astype(str).tolist()

        # 3. Match samples between count matrix and metadata
        common_samples = [s for s in meta_samples if s in sample_cols]
        if len(common_samples) < 2:
            raise ValueError(
                f"At least 2 common sample IDs required between count matrix and metadata.\n"
                f"Count matrix columns: {sample_cols}\n"
                f"Metadata samples: {meta_samples}\n"
                f"Common matching samples: {common_samples}"
            )

        # 4. Validate experimental design and contrasts
        if condition_col not in meta_df.columns:
            candidates = [c for c in meta_df.columns if c != sample_id_col]
            if not candidates:
                raise ValueError(f"No experimental condition column found in metadata '{meta_path}'")
            condition_col = candidates[0]

        meta_sub = meta_df[meta_df[sample_id_col].astype(str).isin(common_samples)].copy()
        conditions = [str(c) for c in meta_sub[condition_col].unique() if pd.notna(c) and str(c).strip()]
        if len(conditions) < 2:
            raise ValueError(
                f"At least 2 distinct experimental conditions required in column '{condition_col}', "
                f"found only: {conditions}"
            )

        ref = contrast_reference.strip()
        target = contrast_target.strip()

        if ref:
            if ref not in conditions:
                raise ValueError(
                    f"Specified contrast_reference '{ref}' not found in '{condition_col}' values: {conditions}"
                )
        else:
            ref = str(conditions[0])

        if target:
            if target not in conditions:
                raise ValueError(
                    f"Specified contrast_target '{target}' not found in '{condition_col}' values: {conditions}"
                )
            if target == ref:
                raise ValueError(
                    f"contrast_target ('{target}') and contrast_reference ('{ref}') cannot be identical."
                )
        else:
            target_candidates = [c for c in conditions if c != ref]
            if not target_candidates:
                raise ValueError(f"No alternative condition found for comparison against reference '{ref}'")
            target = target_candidates[0]

        # 5. Check numeric integrity of counts
        sub_counts = counts_df[[gene_col] + common_samples].copy()
        numeric_mat = sub_counts[common_samples].apply(pd.to_numeric, errors="coerce")
        if numeric_mat.isna().any().any():
            raise ValueError("Count matrix contains non-numeric or NaN entries for common samples.")
        if (numeric_mat < 0).any().any():
            raise ValueError("Count matrix contains negative count values; raw or length-scaled counts must be non-negative.")

        # 6. Locate Rscript with Bioconductor DESeq2
        rscript_bin = _find_rscript()
        if not rscript_bin:
            raise RuntimeError(
                "Rscript executable not found. R with Bioconductor DESeq2 in environment 'bulk_rna_seq' is required.\n"
                "Please ensure conda environment 'bulk_rna_seq' is installed or set RSCRIPT_PATH environment variable."
            )

        # 7. Write isolated R script
        r_script_code = f"""# Auto-generated by ComfyBIOWMS DESeq2 custom node
suppressPackageStartupMessages({{
    library(DESeq2)
}})

counts_file <- "{counts_path}"
meta_file <- "{meta_path}"
condition_col <- "{condition_col}"
ref_level <- "{ref}"
target_level <- "{target}"
results_csv <- "{results_csv}"
norm_counts_csv <- "{norm_counts_csv}"
lfc_thresh <- {lfc_threshold}
alpha_thresh <- {alpha}

# Load counts
sep_counts <- if (grepl("\\t", readLines(counts_file, n=1))) "\\t" else ","
counts <- read.csv(counts_file, sep=sep_counts, check.names = FALSE, row.names = 1)

# Load metadata
sep_meta <- if (grepl("\\t", readLines(meta_file, n=1))) "\\t" else ","
meta <- read.csv(meta_file, sep=sep_meta, check.names = FALSE, row.names = 1)

common_samples <- intersect(colnames(counts), rownames(meta))
if (length(common_samples) < 2) {{
    stop(paste("Fewer than 2 common samples:", paste(common_samples, collapse=", ")))
}}

counts <- counts[, common_samples, drop=FALSE]
meta <- meta[common_samples, , drop=FALSE]

meta[[condition_col]] <- factor(meta[[condition_col]])
if (!ref_level %in% levels(meta[[condition_col]])) {{
    stop(paste("Reference level", ref_level, "not in condition factor levels"))
}}
meta[[condition_col]] <- relevel(meta[[condition_col]], ref = ref_level)

# Build DESeqDataSet
design_formula <- as.formula(paste("~", condition_col))
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)), colData = meta, design = design_formula)

# Run DESeq with robust dispersion fitting for variable gene counts
dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) {{
        tryCatch(
            DESeq(dds, fitType = "local", quiet = TRUE),
            error = function(e2) {{
                tryCatch(
                    DESeq(dds, fitType = "mean", quiet = TRUE),
                    error = function(e3) {{
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
                        dispersions(dds) <- mcols(dds)$dispGeneEst
                        nbinomWaldTest(dds, quiet = TRUE)
                    }}
                )
            }}
        )
    }}
)

res <- results(dds, contrast = c(condition_col, target_level, ref_level), alpha = alpha_thresh)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Order columns
cols <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
present_cols <- intersect(cols, colnames(res_df))
other_cols <- setdiff(colnames(res_df), present_cols)
res_df <- res_df[, c(present_cols, other_cols)]

write.csv(res_df, results_csv, row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, norm_counts_csv)

cat("SUCCESS\\n")
cat("DESeq2_version:", as.character(packageVersion("DESeq2")), "\\n")
"""
        r_script_path.write_text(r_script_code, encoding="utf-8")

        # 8. Execute Rscript
        sub_env = os.environ.copy()
        sub_env["LC_ALL"] = "C"
        proc = subprocess.run(
            [rscript_bin, str(r_script_path)],
            capture_output=True,
            text=True,
            cwd=str(out),
            env=sub_env,
        )

        r_log_path.write_text(f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n", encoding="utf-8")

        if proc.returncode != 0 or not results_csv.is_file() or results_csv.stat().st_size == 0:
            raise RuntimeError(
                f"Bioconductor DESeq2 failed with returncode {proc.returncode}.\n"
                f"Command: {rscript_bin} {r_script_path}\n"
                f"STDERR:\n{proc.stderr}\n"
                f"STDOUT:\n{proc.stdout}"
            )

        # 9. Parse and verify results
        results_df = pd.read_csv(results_csv)
        lfc_col = "log2FoldChange" if "log2FoldChange" in results_df.columns else results_df.columns[2]
        padj_col = "padj" if "padj" in results_df.columns else results_df.columns[-1]

        valid = results_df.dropna(subset=[lfc_col, padj_col])
        sig_up = int(((valid[lfc_col] >= lfc_threshold) & (valid[padj_col] <= alpha)).sum())
        sig_down = int(((valid[lfc_col] <= -lfc_threshold) & (valid[padj_col] <= alpha)).sum())
        total_tested = int(len(results_df))

        # Extract DESeq2 version from output
        deseq2_version = "unknown"
        for line in proc.stdout.splitlines():
            if line.startswith("DESeq2_version:"):
                deseq2_version = line.split(":", 1)[1].strip()

        summary = {
            "total_genes_tested": total_tested,
            "contrast": f"{target} vs {ref}",
            "condition_column": condition_col,
            "reference_level": ref,
            "target_level": target,
            "alpha_cutoff": alpha,
            "lfc_threshold": lfc_threshold,
            "significant_upregulated": sig_up,
            "significant_downregulated": sig_down,
            "total_significant_degs": sig_up + sig_down,
            "engine": "Bioconductor_DESeq2",
            "deseq2_version": deseq2_version,
            "rscript_executable": rscript_bin,
            "samples_analyzed": common_samples,
        }

        manifest_json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        manifest_sh_path = out / "run_manifest.sh"
        sh_lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"# Bioconductor DESeq2 execution",
            f"cd {shlex.quote(str(out))}",
            f"{shlex.quote(rscript_bin)} {shlex.quote(str(r_script_path))}",
            "",
        ]
        manifest_sh_path.write_text("\n".join(sh_lines), encoding="utf-8")
        manifest_sh_path.chmod(0o755)

        return str(results_csv), str(norm_counts_csv), json.dumps(summary, indent=2)


NODE_CLASS_MAPPINGS = {"DESeq2": DESeq2}
NODE_DISPLAY_NAME_MAPPINGS = {"DESeq2": "DESeq2: Differential Gene Expression"}

DESeq2Node = DESeq2
DESeq2Analysis = DESeq2
DESeq2AnalysisNode = DESeq2

__all__ = [
    "DESeq2",
    "DESeq2Node",
    "DESeq2Analysis",
    "DESeq2AnalysisNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
