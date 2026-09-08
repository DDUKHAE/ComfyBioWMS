"""cyvcf2 high-performance VCF/BCF parsing and statistics node.

Python packages: cyvcf2
External binaries: none
"""

import json
from pathlib import Path


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


class Cyvcf2Stats:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Genomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT", "INT")
    RETURN_NAMES = ("summary_json", "total_variants", "snps_count", "indels_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "vcf_file": ("STRING", {"default": ""}),
            },
        }

    def run(self, vcf_file: str):
        import cyvcf2
        vcf_path = _file(vcf_file, "Input VCF/BCF")
        vcf = cyvcf2.VCF(str(vcf_path))
        n_snps = 0
        n_indels = 0
        n_other = 0
        samples = list(vcf.samples)
        n_transitions = 0
        n_transversions = 0
        transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        transversions = {
            ("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
            ("C", "G"), ("G", "C"), ("G", "T"), ("T", "G"),
        }

        for var in vcf:
            if var.is_snp:
                n_snps += 1
                ref = str(var.REF).upper() if var.REF else ""
                alt_list = var.ALT or []
                for alt_allele in alt_list:
                    alt = str(alt_allele).upper()
                    if len(ref) == 1 and len(alt) == 1:
                        pair = (ref, alt)
                        if pair in transitions:
                            n_transitions += 1
                        elif pair in transversions:
                            n_transversions += 1
            elif var.is_indel:
                n_indels += 1
            else:
                n_other += 1

        total = n_snps + n_indels + n_other
        ti_tv = (n_transitions / n_transversions) if n_transversions > 0 else (float(n_transitions) if n_transitions > 0 else None)
        summary = {
            "total_variants": total,
            "snps_count": n_snps,
            "indels_count": n_indels,
            "other_count": n_other,
            "transitions_count": n_transitions,
            "transversions_count": n_transversions,
            "samples": samples,
            "ti_tv_ratio": round(ti_tv, 4) if ti_tv is not None else None,
        }

        return (json.dumps(summary, indent=2), total, n_snps, n_indels)


NODE_CLASS_MAPPINGS = {"Cyvcf2Stats": Cyvcf2Stats}
NODE_DISPLAY_NAME_MAPPINGS = {"Cyvcf2Stats": "cyvcf2: Variant File Statistics"}


# Backward compatibility aliases
Cyvcf2StatsNode = Cyvcf2Stats

__all__ = ["Cyvcf2Stats",
    "Cyvcf2StatsNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
