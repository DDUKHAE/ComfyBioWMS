"""Tximport transcript-to-gene quantification matrix aggregator node.

Python packages: pandas
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



class Tximport:
    CATEGORY = "ComfyBIO/Quantification"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("gene_counts_tsv", "gene_tpm_tsv", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "quant_files": ("STRING", {"default": ""}),
            },
            "optional": {
                "tx2gene_tsv": ("STRING", {"default": ""}),
                "sample_names": ("STRING", {"default": ""}),
            },
        }

    def run(
        self,
        quant_files: str,
        output_dir: str = "",
        tx2gene_tsv: str = "",
        sample_names: str = "",
    ):
        file_paths = []
        for raw in quant_files.split(","):
            token = raw.strip()
            if token:
                file_paths.append(_file(token, "Salmon quant.sf"))

        if not file_paths:
            raise ValueError("No valid quant.sf files provided to Tximport")

        # Parse sample names
        samples = []
        if sample_names.strip():
            samples = [s.strip() for s in sample_names.split(",") if s.strip()]

        if len(samples) != len(file_paths):
            # Derive sample name from parent folder name or file stem
            samples = []
            for p in file_paths:
                parent_name = p.parent.name
                if parent_name and parent_name not in ("quant", ".", ""):
                    samples.append(parent_name)
                else:
                    samples.append(p.stem)

        # Load tx2gene mapping if provided
        tx2gene = {}
        if tx2gene_tsv.strip():
            map_path = _file(tx2gene_tsv, "tx2gene TSV")
            with open(map_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t" if "\t" in f.readline() else ",")
                f.seek(0)
                for row in reader:
                    if len(row) >= 2 and row[0].strip() and not row[0].startswith("#"):
                        tx2gene[row[0].strip()] = row[1].strip()

        # Parse each quant.sf
        # Columns: Name, Length, EffectiveLength, TPM, NumReads
        counts_by_sample = {}
        tpm_by_sample = {}
        all_features = set()

        for s_name, path in zip(samples, file_paths):
            counts = {}
            tpms = {}
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    tx_id = row.get("Name", "").strip()
                    feature_id = tx2gene.get(tx_id, tx_id) if tx2gene else tx_id
                    try:
                        read_count = float(row.get("NumReads", 0.0))
                        tpm_val = float(row.get("TPM", 0.0))
                    except ValueError:
                        continue

                    counts[feature_id] = counts.get(feature_id, 0.0) + read_count
                    tpms[feature_id] = tpms.get(feature_id, 0.0) + tpm_val
                    all_features.add(feature_id)

            counts_by_sample[s_name] = counts
            tpm_by_sample[s_name] = tpms

        sorted_features = sorted(all_features)
        out = _output_dir("Tximport", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        counts_tsv = out / "gene_counts.tsv"
        tpm_tsv = out / "gene_tpm.tsv"

        # Write counts TSV
        with open(counts_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id"] + samples)
            for feat in sorted_features:
                row = [feat] + [
                    f"{counts_by_sample[s].get(feat, 0.0):.4f}" for s in samples
                ]
                writer.writerow(row)

        # Write TPM TSV
        with open(tpm_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id"] + samples)
            for feat in sorted_features:
                row = [feat] + [
                    f"{tpm_by_sample[s].get(feat, 0.0):.4f}" for s in samples
                ]
                writer.writerow(row)

        summary = {
            "num_samples": len(samples),
            "samples": samples,
            "num_features": len(sorted_features),
            "mapped_with_tx2gene": bool(tx2gene),
        }

        return str(counts_tsv), str(tpm_tsv), json.dumps(summary)


NODE_CLASS_MAPPINGS = {"Tximport": Tximport}
NODE_DISPLAY_NAME_MAPPINGS = {"Tximport": "tximport: Transcript to Gene Count Matrix"}


# Backward compatibility aliases
TximportNode = Tximport

__all__ = ["Tximport",
    "TximportNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
