"""Pysam BAM/SAM/CRAM statistics and inspection node.

Python packages: pysam
External binaries: none
"""

import json
from pathlib import Path


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


class PysamStats:
    CATEGORY = "ComfyBIO/Genomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT", "INT", "FLOAT")
    RETURN_NAMES = ("summary_json", "total_reads", "mapped_reads", "mapping_rate_pct")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bam_file": ("STRING", {"default": ""}),
            },
        }

    def run(self, bam_file: str):
        import pysam
        bam_path = _file(bam_file, "Input BAM/CRAM")

        with pysam.AlignmentFile(str(bam_path), "rb") as af:
            total = af.mapped + af.unmapped
            mapped = af.mapped
            rate = round(100.0 * mapped / total, 2) if total > 0 else 0.0
            references = list(af.references)
            lengths = list(af.lengths)

        summary = {
            "total_reads": total,
            "mapped_reads": mapped,
            "unmapped_reads": total - mapped,
            "mapping_rate_pct": rate,
            "n_references": len(references),
            "references": references[:20],
        }

        return (json.dumps(summary, indent=2), total, mapped, rate)


NODE_CLASS_MAPPINGS = {"PysamStats": PysamStats}
NODE_DISPLAY_NAME_MAPPINGS = {"PysamStats": "Pysam: Alignment Statistics & Summary"}


# Backward compatibility aliases
PysamStatsNode = PysamStats

__all__ = ["PysamStats",
    "PysamStatsNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
