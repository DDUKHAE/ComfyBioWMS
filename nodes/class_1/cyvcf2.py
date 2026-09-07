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

        for var in vcf:
            if var.is_snp:
                n_snps += 1
            elif var.is_indel:
                n_indels += 1
            else:
                n_other += 1

        total = n_snps + n_indels + n_other
        summary = {
            "total_variants": total,
            "snps_count": n_snps,
            "indels_count": n_indels,
            "other_count": n_other,
            "samples": samples,
            "ti_tv_ratio": None,
        }

        return (json.dumps(summary, indent=2), total, n_snps, n_indels)


NODE_CLASS_MAPPINGS = {"Cyvcf2Stats": Cyvcf2Stats}
NODE_DISPLAY_NAME_MAPPINGS = {"Cyvcf2Stats": "cyvcf2: Variant File Statistics"}


# Backward compatibility aliases
Cyvcf2StatsNode = Cyvcf2Stats

__all__ = ["Cyvcf2Stats",
    "Cyvcf2StatsNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
