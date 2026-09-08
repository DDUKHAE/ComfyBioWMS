"""pybedtools genomic interval operations node.

Python packages: pybedtools
External binaries: bedtools
"""

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



class PybedtoolsIntersect:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Genomics"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("intersected_bed", "feature_count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bed_a": ("STRING", {"default": ""}),
                "bed_b": ("STRING", {"default": ""}),
            },
            "optional": {
                "write_original_a": ("BOOLEAN", {"default": True}),
            },
        }

    def run(self, bed_a: str, bed_b: str, output_dir: str = "", write_original_a: bool = True):
        import pybedtools
        a_path = _file(bed_a, "BED file A")
        b_path = _file(bed_b, "BED file B")

        out = _output_dir("PybedtoolsIntersect", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        dest = out / f"{a_path.stem}_intersect_{b_path.stem}.bed"

        a = pybedtools.BedTool(str(a_path))
        b = pybedtools.BedTool(str(b_path))
        c = a.intersect(b, u=write_original_a)
        c.saveas(str(dest))

        count = len(c)
        return (str(dest), count)


NODE_CLASS_MAPPINGS = {"PybedtoolsIntersect": PybedtoolsIntersect}
NODE_DISPLAY_NAME_MAPPINGS = {"PybedtoolsIntersect": "pybedtools: Interval Intersect"}


# Backward compatibility aliases
PybedtoolsIntersectNode = PybedtoolsIntersect

__all__ = ["PybedtoolsIntersect",
    "PybedtoolsIntersectNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
