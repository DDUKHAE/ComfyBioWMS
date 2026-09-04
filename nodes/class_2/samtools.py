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
import threading
from pathlib import Path


_SORT_MANAGED = {"-@": 1, "-m": 1, "-n": 0, "-N": 0, "-o": 1, "-O": 1}
_INDEX_MANAGED = {"-@": 1, "-b": 0, "-c": 0, "-o": 1}
_MARKDUP_MANAGED = {"-@": 1, "-r": 0, "--mode": 1, "-d": 1, "-f": 1}


def _file(value, label):
    path = Path(value).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _executable():
    executable = shutil.which("samtools")
    if not executable:
        raise RuntimeError(
            "samtools executable not found on PATH; install conda package samtools=1.24"
        )
    return executable


def _filter_extra(text, managed):
    tokens, kept, ignored = shlex.split(text), [], []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        match = next(
            (
                flag
                for flag in managed
                if token == flag
                or token.startswith(flag + "=")
                or (
                    len(flag) == 2
                    and managed[flag] == 1
                    and token.startswith(flag)
                    and token != flag
                )
            ),
            None,
        )
        if match:
            ignored.append(token)
            i += 1 + (managed[match] if token == match else 0)
        else:
            kept.append(token)
            i += 1
    return kept, ignored


def _run(argv, cwd, partial=None):
    process = subprocess.Popen(
        argv,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
    )
    stdout, stderr = [], []

    def drain(pipe, target, collected):
        for line in iter(pipe.readline, ""):
            collected.append(line)
            print(line, end="", file=target, flush=True)

    threads = [
        threading.Thread(target=drain, args=(process.stdout, sys.stdout, stdout)),
        threading.Thread(target=drain, args=(process.stderr, sys.stderr, stderr)),
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    code = process.wait()
    if code:
        if partial:
            Path(partial).unlink(missing_ok=True)
        raise RuntimeError(f"samtools exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")
    return "".join(stdout)


def _extras(text, managed, command):
    kept, ignored = _filter_extra(text, managed)
    if ignored:
        print(
            f"[samtools {command}] ignored node-managed extra options: {' '.join(ignored)}",
            file=sys.stderr,
        )
    return kept


class SamtoolsSortNode:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("sorted_bam",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_alignment": ("STRING", {"default": ""}),
                "output_bam": ("STRING", {"default": "alignment/sorted.bam"}),
            },
            "optional": {
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "sort_order": (["coordinate", "name", "lexicographical"], {"default": "coordinate"}),
                "memory_per_thread": ("STRING", {"default": "768M"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_alignment, output_bam, threads=1, sort_order="coordinate", memory_per_thread="768M", extra_command=""):
        source = _file(input_alignment, "Input alignment")
        executable = _executable()
        output = Path(output_bam).expanduser()
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "sort", "-@", str(threads), "-m", memory_per_thread, "-O", "BAM", "-o", str(output)]
        if sort_order == "name":
            argv.append("-n")
        elif sort_order == "lexicographical":
            argv.append("-N")
        argv += _extras(extra_command, _SORT_MANAGED, "sort") + [str(source)]
        _run(argv, output.parent, output)
        _run([executable, "quickcheck", str(output)], output.parent)
        return (str(output),)


class SamtoolsIndexNode:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("alignment_index",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {"input_bam": ("STRING", {"default": ""})},
            "optional": {
                "output_index": ("STRING", {"default": ""}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "index_format": (["bai", "csi"], {"default": "bai"}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bam, output_index="", threads=1, index_format="bai", extra_command=""):
        source = _file(input_bam, "Input BAM")
        executable = _executable()
        output = Path(output_index).expanduser() if output_index else Path(str(source) + (".csi" if index_format == "csi" else ".bai"))
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "index", "-@", str(threads), "-c" if index_format == "csi" else "-b", "-o", str(output)]
        argv += _extras(extra_command, _INDEX_MANAGED, "index") + [str(source)]
        _run(argv, output.parent, output)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"samtools did not create a nonempty index: {output}")
        _run([executable, "idxstats", str(source)], output.parent)
        return (str(output),)


class SamtoolsMarkdupNode:
    CATEGORY = "ComfyBIO/Alignment"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("marked_bam",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_bam": ("STRING", {"default": ""}),
                "output_bam": ("STRING", {"default": "alignment/marked.bam"}),
            },
            "optional": {
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "remove_duplicates": ("BOOLEAN", {"default": False}),
                "mode": (["template", "sequence"], {"default": "template"}),
                "optical_distance": ("INT", {"default": 100, "min": 0}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bam, output_bam, threads=1, remove_duplicates=False, mode="template", optical_distance=100, extra_command=""):
        source = _file(input_bam, "Input BAM")
        executable = _executable()
        output = Path(output_bam).expanduser()
        output.parent.mkdir(parents=True, exist_ok=True)
        extra = _extras(extra_command, _MARKDUP_MANAGED, "markdup")
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
    "SamtoolsSortNode": SamtoolsSortNode,
    "SamtoolsIndexNode": SamtoolsIndexNode,
    "SamtoolsMarkdupNode": SamtoolsMarkdupNode,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "SamtoolsSortNode": "samtools: Sort Alignment",
    "SamtoolsIndexNode": "samtools: Index Alignment",
    "SamtoolsMarkdupNode": "samtools: Mark Duplicates",
}

__all__ = [*NODE_CLASS_MAPPINGS, "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
