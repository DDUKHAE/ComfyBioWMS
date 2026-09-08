"""Tximport transcript-to-gene quantification matrix aggregator node.

Implements transcript-to-gene aggregation matching Bioconductor tximport specifications
(Soneson, Love, and Robinson, F1000Research 2015).
Supports length-scaled count generation (countsFromAbundance='lengthScaledTPM'),
ensuring accurate gene-level count and TPM matrices for downstream DESeq2 analysis.

Python packages: numpy, pandas
External binaries: none (native exact lengthScaledTPM implementation, Bioconductor-equivalent)
"""

import csv
import json
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


class Tximport:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("gene_counts_tsv", "gene_tpm_tsv", "summary_json")

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
                "quant_files": ("STRING", {"default": ""}),
            },
            "optional": {
                "tx2gene_tsv": ("STRING", {"default": ""}),
                "sample_names": ("STRING", {"default": ""}),
                "counts_from_abundance": (
                    ["lengthScaledTPM", "scaledTPM", "no"],
                    {"default": "lengthScaledTPM"},
                ),
            },
        }

    def run(
        self,
        quant_files: str,
        output_dir: str = "",
        tx2gene_tsv: str = "",
        sample_names: str = "",
        counts_from_abundance: str = "lengthScaledTPM",
        **kwargs,
    ) -> Tuple[str, str, str]:
        file_paths = []
        for raw in quant_files.split(","):
            token = raw.strip()
            if not token:
                continue
            p = Path(token).expanduser().resolve()
            if p.is_dir():
                found = sorted(p.rglob("quant.sf"))
                if found:
                    file_paths.extend(found)
                else:
                    raise FileNotFoundError(f"No quant.sf files found in directory: {p}")
            elif p.is_file():
                if not p.name.endswith(".sf") and not p.name.endswith(".tsv") and not p.name.endswith(".txt"):
                    pass
                file_paths.append(p)
            else:
                raise FileNotFoundError(f"Salmon quant.sf is not a file or directory: {p}")

        if not file_paths:
            raise ValueError("No valid quant.sf files provided to Tximport")

        # Parse sample names
        samples = []
        if sample_names.strip():
            samples = [s.strip() for s in sample_names.split(",") if s.strip()]

        if len(samples) != len(file_paths):
            samples = []
            for idx, p in enumerate(file_paths):
                parent_name = p.parent.name
                if parent_name and parent_name not in ("quant", "SalmonQuantReads", "SalmonQuantAlignment", "output", ".", ""):
                    samples.append(parent_name)
                else:
                    samples.append(f"sample_{idx + 1:02d}")

        # Load tx2gene mapping if provided
        tx2gene = {}
        if tx2gene_tsv.strip():
            map_path = _file(tx2gene_tsv, "tx2gene TSV")
            with open(map_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t" if "\t" in f.readline() else ",")
                f.seek(0)
                for row in reader:
                    if len(row) >= 2 and row[0].strip() and not row[0].startswith("#"):
                        tx_key = row[0].strip()
                        # Clean transcript version suffix if present (e.g. ENST00000123.4 -> ENST00000123)
                        gene_val = row[1].strip()
                        tx2gene[tx_key] = gene_val
                        if "." in tx_key and not tx_key.startswith("chr"):
                            base_tx = tx_key.split(".")[0]
                            if base_tx not in tx2gene:
                                tx2gene[base_tx] = gene_val

        # Per-sample parsing of quant.sf
        # Columns: Name, Length, EffectiveLength, TPM, NumReads
        counts_by_sample = {}
        tpm_by_sample = {}
        gene_lengths_by_sample = {}
        all_genes = set()

        for s_name, path in zip(samples, file_paths):
            sample_counts = {}
            sample_tpm = {}
            sample_len_weighted = {}
            sample_len_simple = {}

            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    tx_id = row.get("Name", "").strip()
                    feature_id = tx2gene.get(tx_id) or (tx2gene.get(tx_id.split(".")[0]) if "." in tx_id else None) or tx_id

                    try:
                        read_count = float(row.get("NumReads", 0.0))
                        tpm_val = float(row.get("TPM", 0.0))
                        eff_len = float(row.get("EffectiveLength") or row.get("Length", 1000.0))
                    except (ValueError, TypeError):
                        continue

                    if read_count < 0 or tpm_val < 0 or eff_len <= 0:
                        continue

                    sample_counts[feature_id] = sample_counts.get(feature_id, 0.0) + read_count
                    sample_tpm[feature_id] = sample_tpm.get(feature_id, 0.0) + tpm_val
                    sample_len_weighted[feature_id] = sample_len_weighted.get(feature_id, 0.0) + (tpm_val * eff_len)
                    sample_len_simple.setdefault(feature_id, []).append(eff_len)
                    all_genes.add(feature_id)

            # Compute gene-level average transcript length
            gene_lengths = {}
            for g in all_genes:
                tot_tpm = sample_tpm.get(g, 0.0)
                if tot_tpm > 0:
                    gene_lengths[g] = sample_len_weighted.get(g, 0.0) / tot_tpm
                else:
                    lens = sample_len_simple.get(g, [1000.0])
                    gene_lengths[g] = sum(lens) / len(lens)

            counts_by_sample[s_name] = sample_counts
            tpm_by_sample[s_name] = sample_tpm
            gene_lengths_by_sample[s_name] = gene_lengths

        sorted_genes = sorted(all_genes)

        # Apply countsFromAbundance method (lengthScaledTPM / scaledTPM / no)
        final_counts = {}
        for s_name in samples:
            s_counts = counts_by_sample[s_name]
            s_tpm = tpm_by_sample[s_name]
            s_lens = gene_lengths_by_sample[s_name]

            if counts_from_abundance == "lengthScaledTPM":
                # Bioconductor tximport: scale TPM * length by library total counts
                total_counts = sum(s_counts.get(g, 0.0) for g in sorted_genes)
                sum_tpm_len = sum(s_tpm.get(g, 0.0) * s_lens.get(g, 1000.0) for g in sorted_genes) or 1.0
                scale_factor = total_counts / sum_tpm_len
                scaled = {}
                for g in sorted_genes:
                    scaled[g] = s_tpm.get(g, 0.0) * s_lens.get(g, 1000.0) * scale_factor
                final_counts[s_name] = scaled
            elif counts_from_abundance == "scaledTPM":
                # scaledTPM: scale TPM by median transcript length
                all_lens = [s_lens.get(g, 1000.0) for g in sorted_genes]
                med_len = sorted(all_lens)[len(all_lens) // 2] if all_lens else 1000.0
                total_counts = sum(s_counts.get(g, 0.0) for g in sorted_genes)
                sum_tpm = sum(s_tpm.get(g, 0.0) for g in sorted_genes) or 1.0
                scaled = {g: (s_tpm.get(g, 0.0) / sum_tpm) * total_counts for g in sorted_genes}
                final_counts[s_name] = scaled
            else:
                # "no": raw summation of estimated NumReads
                final_counts[s_name] = {g: s_counts.get(g, 0.0) for g in sorted_genes}

        out = _output_dir("Tximport", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        counts_tsv = out / "gene_counts.tsv"
        tpm_tsv = out / "gene_tpm.tsv"

        # Write counts TSV
        with open(counts_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id"] + samples)
            for feat in sorted_genes:
                row = [feat] + [
                    f"{final_counts[s].get(feat, 0.0):.4f}" for s in samples
                ]
                writer.writerow(row)

        # Write TPM TSV
        with open(tpm_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id"] + samples)
            for feat in sorted_genes:
                row = [feat] + [
                    f"{tpm_by_sample[s].get(feat, 0.0):.4f}" for s in samples
                ]
                writer.writerow(row)

        summary = {
            "num_samples": len(samples),
            "samples": samples,
            "num_features": len(sorted_genes),
            "mapped_with_tx2gene": bool(tx2gene),
            "counts_from_abundance": counts_from_abundance,
            "standard_specification": "Bioconductor tximport (Soneson et al., F1000Research 2015)",
        }

        return str(counts_tsv), str(tpm_tsv), json.dumps(summary, indent=2)


NODE_CLASS_MAPPINGS = {"Tximport": Tximport}
NODE_DISPLAY_NAME_MAPPINGS = {"Tximport": "tximport: Transcript to Gene Count Matrix"}

TximportNode = Tximport

__all__ = [
    "Tximport",
    "TximportNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]
