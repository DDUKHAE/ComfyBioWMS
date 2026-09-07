"""FastQ Concatenation node for merging lane-split or multi-run fastq files.

Python packages: none
External binaries: cat (or gzip / python gzip stream)
"""

import shutil
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



class CatFastq:
    CATEGORY = "ComfyBIO/Preprocessing"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("merged_fastq",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "fastq_files": ("STRING", {"default": ""}),
            },
            "optional": {
                "output_filename": ("STRING", {"default": "merged.fastq.gz"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, fastq_files: str, output_dir: str = "", output_filename: str = "merged.fastq.gz", extra_command: str = ""):
        file_list = []
        for raw in fastq_files.split(","):
            token = raw.strip()
            if token:
                file_list.append(_file(token, "Input FastQ"))

        if not file_list:
            raise ValueError("No valid fastq files provided to CatFastq")

        out = _output_dir("CatFastq", output_dir)
        out.mkdir(parents=True, exist_ok=True)
        dest = out / output_filename

        with open(dest, "wb") as wfd:
            for fpath in file_list:
                with open(fpath, "rb") as rfd:
                    shutil.copyfileobj(rfd, wfd, length=1024 * 1024)

        if not dest.is_file() or dest.stat().st_size == 0:
            raise RuntimeError(f"CatFastq failed to create nonempty merged file at {dest}")

        return (str(dest),)


NODE_CLASS_MAPPINGS = {"CatFastq": CatFastq}
NODE_DISPLAY_NAME_MAPPINGS = {"CatFastq": "FastQ: Concatenate Split Files"}


# Backward compatibility aliases
CatFastqNode = CatFastq

__all__ = ["CatFastq",
    "CatFastqNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
