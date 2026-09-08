"""Path collector node.

ComfyUI ports carry one value, so a per-sample tool (Salmon, RSEM, ...) cannot
feed a multi-sample aggregator (Tximport) directly. This joins up to 8 upstream
path outputs into the comma-separated string those aggregators already parse,
and — unlike typing the paths into a widget — keeps the DAG dependency, so the
aggregator waits for every branch.

Python packages: none
External binaries: none
"""

MAX_INPUTS = 8


class JoinPaths:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Utility"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("joined_paths", "count")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "path_1": ("STRING", {"default": ""}),
            },
            "optional": {
                **{f"path_{i}": ("STRING", {"default": ""}) for i in range(2, MAX_INPUTS + 1)},
                "separator": ("STRING", {"default": ","}),
            },
        }

    def run(self, path_1: str = "", separator: str = ",", **kwargs):
        ordered = [path_1] + [kwargs.get(f"path_{i}", "") for i in range(2, MAX_INPUTS + 1)]
        paths = [p.strip() for p in ordered if isinstance(p, str) and p.strip()]
        if not paths:
            raise ValueError("JoinPaths received no non-empty path")
        return (separator.join(paths), len(paths))


JoinPathsNode = JoinPaths

NODE_CLASS_MAPPINGS = {"JoinPaths": JoinPaths}
NODE_DISPLAY_NAME_MAPPINGS = {"JoinPaths": "Join Paths (multi-sample collector)"}

__all__ = [
    "JoinPaths",
    "JoinPathsNode",
    "NODE_CLASS_MAPPINGS",
    "NODE_DISPLAY_NAME_MAPPINGS",
]


if __name__ == "__main__":
    n = JoinPaths()
    assert n.run("a/quant.sf", path_2="b/quant.sf") == ("a/quant.sf,b/quant.sf", 2)
    # holes in the middle collapse; a disconnected optional port is just ""
    assert n.run("a", path_3=" c ", path_5="") == ("a,c", 2)
    assert n.run("a", separator=";", path_2="b") == ("a;b", 2)
    try:
        n.run("   ")
    except ValueError:
        pass
    else:
        raise AssertionError("empty input must raise")
    print("JoinPaths self-check OK")
