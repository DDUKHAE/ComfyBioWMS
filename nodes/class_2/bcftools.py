"""bcftools mpileup, variant calling, and filtering nodes.

Python packages: none
External binary: bcftools==1.24 (Conda: bioconda::bcftools=1.24;
Apt: bcftools)
Galaxy wrappers: galaxyproject/tools-iuc@6a1769b029357f74e43c73b3da515b5e6f02a608,
tools/bcftools/bcftools_mpileup.xml, bcftools_call.xml, and bcftools_filter.xml
"""

import shlex
import shutil
import subprocess
import sys
import threading
from pathlib import Path


_MPILEUP_MANAGED = {
    "-f": 1, "--fasta-ref": 1, "-d": 1, "--max-depth": 1,
    "-Q": 1, "--min-BQ": 1, "-q": 1, "--min-MQ": 1,
    "--threads": 1, "-O": 1, "--output-type": 1, "-o": 1, "--output": 1,
}
_CALL_MANAGED = {
    "-m": 0, "--multiallelic-caller": 0, "-c": 0, "--consensus-caller": 0,
    "-v": 0, "--variants-only": 0, "--ploidy": 1, "--threads": 1,
    "-O": 1, "--output-type": 1, "-o": 1, "--output": 1,
}
_FILTER_MANAGED = {
    "--threads": 1, "-O": 1, "--output-type": 1, "-o": 1, "--output": 1,
    "-e": 1, "--exclude": 1, "-i": 1, "--include": 1,
    "-s": 1, "--soft-filter": 1, "-g": 1, "--SnpGap": 1,
    "-G": 1, "--IndelGap": 1,
}


def _file(value, label):
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path


def _executable():
    executable = shutil.which("bcftools")
    if not executable:
        raise RuntimeError(
            "bcftools executable not found on PATH; install conda package bcftools=1.24"
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


def _extras(text, managed, command):
    kept, ignored = _filter_extra(text, managed)
    if ignored:
        print(
            f"[bcftools {command}] ignored node-managed extra options: {' '.join(ignored)}",
            file=sys.stderr,
        )
    return kept


def _run(argv, cwd, partial=None):
    process = subprocess.Popen(
        argv, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1
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
        raise RuntimeError(f"bcftools exited {code}: {shlex.join(argv)}\n{''.join(stderr)}")
    return "".join(stdout)


def _validate(executable, output):
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"bcftools did not create a nonempty output: {output}")
    _run([executable, "view", "-h", str(output)], output.parent)


class BcftoolsMpileupNode:
    CATEGORY = "ComfyBIO/Variants"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("pileup_bcf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "reference_fasta": ("STRING", {"default": ""}),
                "input_bam": ("STRING", {"default": ""}),
                "output_bcf": ("STRING", {"default": "variants/pileup.bcf"}),
            },
            "optional": {
                "max_depth": ("INT", {"default": 250, "min": 0}),
                "min_base_quality": ("INT", {"default": 13, "min": 0}),
                "min_mapping_quality": ("INT", {"default": 0, "min": 0}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, reference_fasta, input_bam, output_bcf, max_depth=250, min_base_quality=13, min_mapping_quality=0, threads=1, extra_command=""):
        reference = _file(reference_fasta, "Reference FASTA")
        bam = _file(input_bam, "Input BAM")
        _file(str(reference) + ".fai", "Reference FASTA index")
        executable = _executable()
        output = Path(output_bcf).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "mpileup", "-f", str(reference), "-d", str(max_depth), "-Q", str(min_base_quality), "-q", str(min_mapping_quality), "--threads", str(threads), "-Ob", "-o", str(output)]
        argv += _extras(extra_command, _MPILEUP_MANAGED, "mpileup") + [str(bam)]
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


class BcftoolsCallNode:
    CATEGORY = "ComfyBIO/Variants"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("called_vcf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_bcf": ("STRING", {"default": ""}),
                "output_vcf": ("STRING", {"default": "variants/called.vcf"}),
            },
            "optional": {
                "calling_method": (["multiallelic", "consensus"], {"default": "multiallelic"}),
                "variants_only": ("BOOLEAN", {"default": True}),
                "ploidy": ("STRING", {"default": "default"}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bcf, output_vcf, calling_method="multiallelic", variants_only=True, ploidy="default", threads=1, extra_command=""):
        source = _file(input_bcf, "Input BCF")
        executable = _executable()
        output = Path(output_vcf).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "call", "-m" if calling_method == "multiallelic" else "-c"]
        if variants_only:
            argv.append("-v")
        if ploidy != "default":
            argv += ["--ploidy", ploidy]
        argv += ["--threads", str(threads), "-Ov", "-o", str(output)]
        argv += _extras(extra_command, _CALL_MANAGED, "call") + [str(source)]
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


class BcftoolsFilterNode:
    CATEGORY = "ComfyBIO/Variants"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("filtered_vcf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_vcf": ("STRING", {"default": ""}),
                "output_vcf": ("STRING", {"default": "variants/filtered.vcf"}),
            },
            "optional": {
                "exclude": ("STRING", {"default": "QUAL<10"}),
                "include": ("STRING", {"default": ""}),
                "soft_filter": ("STRING", {"default": ""}),
                "snp_gap": ("INT", {"default": 0, "min": 0}),
                "indel_gap": ("INT", {"default": 0, "min": 0}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_vcf, output_vcf, exclude="QUAL<10", include="", soft_filter="", snp_gap=0, indel_gap=0, threads=1, extra_command=""):
        if include and exclude:
            raise ValueError("include and exclude expressions are mutually exclusive")
        source = _file(input_vcf, "Input VCF/BCF")
        executable = _executable()
        output = Path(output_vcf).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "filter", "--threads", str(threads), "-Ov", "-o", str(output)]
        if exclude:
            argv += ["-e", exclude]
        if include:
            argv += ["-i", include]
        if soft_filter:
            argv += ["-s", soft_filter]
        if snp_gap:
            argv += ["-g", str(snp_gap)]
        if indel_gap:
            argv += ["-G", str(indel_gap)]
        argv += _extras(extra_command, _FILTER_MANAGED, "filter") + [str(source)]
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


NODE_CLASS_MAPPINGS = {
    "BcftoolsMpileupNode": BcftoolsMpileupNode,
    "BcftoolsCallNode": BcftoolsCallNode,
    "BcftoolsFilterNode": BcftoolsFilterNode,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "BcftoolsMpileupNode": "bcftools: Generate Pileup",
    "BcftoolsCallNode": "bcftools: Call Variants",
    "BcftoolsFilterNode": "bcftools: Filter Variants",
}

__all__ = [*NODE_CLASS_MAPPINGS, "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
