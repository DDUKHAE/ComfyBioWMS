"""UCSC toolkit node (bedGraphToBigWig).

Python packages: none
External binaries: bedGraphToBigWig
"""

import shlex
import shutil
import subprocess
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



def _sort_bedgraph_if_needed(in_bg: Path, sorted_bg: Path) -> None:
    with open(in_bg, "r", encoding="utf-8") as rf:
        lines = rf.readlines()
    
    rows = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("track"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            chrom, start, end, val = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            rows.append((chrom, start, end, val))

    rows.sort(key=lambda x: (x[0], x[1]))

    with open(sorted_bg, "w", encoding="utf-8") as wf:
        for r in rows:
            wf.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")


class BedGraphToBigWig:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/Visualization"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("bigwig_file",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "bedgraph_file": ("STRING", {"default": ""}),
                "chrom_sizes": ("STRING", {"default": ""}),
            },
            "optional": {
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, bedgraph_file: str, chrom_sizes: str, output_dir: str = "", extra_command: str = ""):
        bg_path = _file(bedgraph_file, "bedGraph file")
        sizes_path = _file(chrom_sizes, "Chromosome sizes file")

        executable = shutil.which("bedGraphToBigWig")
        if not executable:
            raise RuntimeError(
                "bedGraphToBigWig executable not found on PATH; install bioconda package ucsc-bedgraphtobigwig"
            )

        out = _output_dir("BedGraphToBigWig", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        sorted_bg = out / f"sorted_{bg_path.name}"
        _sort_bedgraph_if_needed(bg_path, sorted_bg)

        out_bw = out / f"{bg_path.stem}.bigWig"

        argv = [executable]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())
        argv += [str(sorted_bg), str(sizes_path), str(out_bw)]

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"bedGraphToBigWig failed ({proc.returncode}): {proc.stderr}")

        if not out_bw.is_file() or out_bw.stat().st_size == 0:
            raise RuntimeError(f"bedGraphToBigWig failed to produce valid bigWig at {out_bw}")

        return (str(out_bw),)


NODE_CLASS_MAPPINGS = {"BedGraphToBigWig": BedGraphToBigWig}
NODE_DISPLAY_NAME_MAPPINGS = {"BedGraphToBigWig": "UCSC: bedGraph to BigWig"}


# Backward compatibility aliases
BedGraphToBigWigNode = BedGraphToBigWig

__all__ = ["BedGraphToBigWig",
    "BedGraphToBigWigNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
