"""GSEApy Gene Set Enrichment Analysis & Enrichr node.

Python packages: gseapy
External binaries: none
"""

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



class GSEAPathway:
    CATEGORY = "ComfyBIO/Pathways"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("enrichment_results_tsv", "summary_json")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "gene_list_file": ("STRING", {"default": ""}),
                "gene_set": (["KEGG_2021_Human", "GO_Biological_Process_2023", "MSigDB_Hallmark_2020", "WikiPathways_2024_Human"], {"default": "KEGG_2021_Human"}),
            },
        }

    def run(self, gene_list_file: str, gene_set: str = "KEGG_2021_Human", output_dir: str = ""):
        import gseapy as gp
        path = _file(gene_list_file, "Gene list file")
        out = _output_dir("GSEAPathway", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        genes = [l.strip().upper() for l in open(path).readlines() if l.strip()]

        enr = gp.enrichr(gene_list=genes, gene_sets=gene_set, outdir=str(out), no_plot=True)
        res_tsv = out / f"{gene_set}.enrichr.reports.txt"
        if not res_tsv.is_file():
            candidates = list(out.glob("*.txt")) + list(out.glob("*.csv"))
            res_tsv = candidates[0] if candidates else out / "enrichment.tsv"
            if not res_tsv.is_file():
                enr.res2d.to_csv(str(res_tsv), sep="\t", index=False)

        summary = {
            "n_genes_input": len(genes),
            "gene_set": gene_set,
            "top_terms": enr.res2d["Term"].head(5).tolist() if hasattr(enr, "res2d") and "Term" in enr.res2d.columns else [],
        }

        return (str(res_tsv), json.dumps(summary, indent=2))


NODE_CLASS_MAPPINGS = {"GSEAPathway": GSEAPathway}
NODE_DISPLAY_NAME_MAPPINGS = {"GSEAPathway": "GSEApy: Gene Set Enrichment Analysis"}


# Backward compatibility aliases
GSEAPathwayNode = GSEAPathway

__all__ = ["GSEAPathway",
    "GSEAPathwayNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
