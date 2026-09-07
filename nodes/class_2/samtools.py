"""samtools sorting, indexing, and duplicate-marking nodes.

Python packages: none
External binary: samtools==1.24 (Conda: bioconda::samtools=1.24;
Apt: samtools)
Galaxy wrappers: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tool_collections/samtools/samtools_sort, samtools_index, and samtools_markdup
"""

import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path



def _file(value, label):
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



def _executable():
    executable = shutil.which("samtools")
    if not executable:
        raise RuntimeError(
            "samtools executable not found on PATH; install conda package samtools=1.24"
        )
    return executable


def _run(argv, cwd, partial=None):
    try:
        subprocess.run(argv, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        if partial:
            Path(partial).unlink(missing_ok=True)
        raise RuntimeError(f"samtools exited {e.returncode}: {shlex.join(argv)}") from e
    return ""


class SamtoolsSort:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("sorted_bam",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_alignment": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "sort_order": (["coordinate", "name", "lexicographical"], {"default": "coordinate"}),
                "memory_per_thread": ("STRING", {"default": "768M"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_alignment, output_bam: str = "", threads=1, sort_order="coordinate", memory_per_thread="768M", extra_command=""):
        source = _file(input_alignment, "Input alignment")
        executable = _executable()
        output = Path(output_bam).expanduser().resolve() if (output_bam and str(output_bam).strip()) else (_output_dir("SamtoolsSort") / f"{_stem(source)}.sorted.bam")
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "sort", "-@", str(threads), "-m", memory_per_thread, "-O", "BAM", "-o", str(output)]
        if sort_order == "name":
            argv.append("-n")
        elif sort_order == "lexicographical":
            argv.append("-N")
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        argv.append(str(source))
        _run(argv, output.parent, output)
        _run([executable, "quickcheck", str(output)], output.parent)
        return (str(output),)


class SamtoolsIndex:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("alignment_index",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {"input_bam": ("STRING", {"default": ""})},
            "optional": {
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "index_format": (["bai", "csi"], {"default": "bai"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bam, output_index="", threads=1, index_format="bai", extra_command=""):
        source = _file(input_bam, "Input BAM")
        executable = _executable()
        output = Path(output_index).expanduser().resolve() if output_index else Path(str(source) + (".csi" if index_format == "csi" else ".bai")).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "index", "-@", str(threads), "-c" if index_format == "csi" else "-b", "-o", str(output)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        argv.append(str(source))
        _run(argv, output.parent, output)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"samtools did not create a nonempty index: {output}")
        _run([executable, "idxstats", str(source)], output.parent)
        return (str(output),)


class SamtoolsMarkdup:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("marked_bam",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_bam": ("STRING", {"default": ""}),
            },
            "optional": {
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "remove_duplicates": ("BOOLEAN", {"default": False}),
                "mode": (["template", "sequence"], {"default": "template"}),
                "optical_distance": ("INT", {"default": 100, "min": 0}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bam, output_bam: str = "", threads=1, remove_duplicates=False, mode="template", optical_distance=100, extra_command=""):
        source = _file(input_bam, "Input BAM")
        executable = _executable()
        output = Path(output_bam).expanduser().resolve() if (output_bam and str(output_bam).strip()) else (_output_dir("SamtoolsMarkdup") / f"{_stem(source)}.markdup.bam")
        output.parent.mkdir(parents=True, exist_ok=True)
        extra = shlex.split(extra_command) if extra_command.strip() else []
        with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
            temporary = Path(temporary)
            namesort = temporary / "namesort.bam"
            fixmate = temporary / "fixmate.bam"
            coordsort = temporary / "coordsort.bam"
            commands = [
                [executable, "collate", "-@", str(threads), "-o", str(namesort), str(source)],
                [executable, "fixmate", "-@", str(threads), "-m", str(namesort), str(fixmate)],
                [executable, "sort", "-@", str(threads), "-O", "BAM", "-o", str(coordsort), str(fixmate)],
                [executable, "markdup", "-@", str(threads), "--mode", {"template": "t", "sequence": "s"}[mode], "-d", str(optical_distance), *(["-r"] if remove_duplicates else []), *extra, str(coordsort), str(output)],
            ]
            for argv in commands:
                _run(argv, output.parent, output if argv is commands[-1] else None)
        _run([executable, "quickcheck", str(output)], output.parent)
        return (str(output),)


NODE_CLASS_MAPPINGS = {
    "SamtoolsSort": SamtoolsSort,
    "SamtoolsIndex": SamtoolsIndex,
    "SamtoolsMarkdup": SamtoolsMarkdup,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "SamtoolsSort": "samtools: Sort Alignment",
    "SamtoolsIndex": "samtools: Index Alignment",
    "SamtoolsMarkdup": "samtools: Mark Duplicates",
}


# Backward compatibility aliases
SamtoolsSortNode = SamtoolsSort
SamtoolsIndexNode = SamtoolsIndex
SamtoolsMarkdupNode = SamtoolsMarkdup

__all__ = [*NODE_CLASS_MAPPINGS, "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
