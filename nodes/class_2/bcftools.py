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
    executable = shutil.which("bcftools")
    if not executable:
        raise RuntimeError(
            "bcftools executable not found on PATH; install conda package bcftools=1.24"
        )
    return executable


def _run(argv, cwd, partial=None):
    try:
        subprocess.run(argv, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        if partial:
            Path(partial).unlink(missing_ok=True)
        raise RuntimeError(f"bcftools exited {e.returncode}: {shlex.join(argv)}") from e
    return ""


def _validate(executable, output):
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"bcftools did not create a nonempty output: {output}")
    _run([executable, "view", "-h", str(output)], output.parent)


class BcftoolsMpileup:
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
            },
            "optional": {
                "max_depth": ("INT", {"default": 250, "min": 0}),
                "min_base_quality": ("INT", {"default": 13, "min": 0}),
                "min_mapping_quality": ("INT", {"default": 0, "min": 0}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, reference_fasta, input_bam, output_bcf: str = "", max_depth=250, min_base_quality=13, min_mapping_quality=0, threads=1, extra_command=""):
        reference = _file(reference_fasta, "Reference FASTA")
        bam = _file(input_bam, "Input BAM")
        _file(str(reference) + ".fai", "Reference FASTA index")
        executable = _executable()
        output = Path(output_bcf).expanduser().resolve() if (output_bcf and str(output_bcf).strip()) else (_output_dir("BcftoolsMpileup") / f"{Path(input_bam).stem}.bcf")
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "mpileup", "-f", str(reference), "-d", str(max_depth), "-Q", str(min_base_quality), "-q", str(min_mapping_quality), "--threads", str(threads), "-Ob", "-o", str(output)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        argv.append(str(bam))
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


class BcftoolsCall:
    CATEGORY = "ComfyBIO/Variants"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("called_vcf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_bcf": ("STRING", {"default": ""}),
            },
            "optional": {
                "calling_method": (["multiallelic", "consensus"], {"default": "multiallelic"}),
                "variants_only": ("BOOLEAN", {"default": True}),
                "ploidy": ("STRING", {"default": "default"}),
                "threads": ("INT", {"default": 1, "min": 1, "max": 256}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, input_bcf, output_vcf: str = "", calling_method="multiallelic", variants_only=True, ploidy="default", threads=1, extra_command=""):
        source = _file(input_bcf, "Input BCF")
        executable = _executable()
        output = Path(output_vcf).expanduser().resolve() if (output_vcf and str(output_vcf).strip()) else (_output_dir("BcftoolsCall") / f"{Path(input_bcf).stem}.calls.vcf.gz")
        output.parent.mkdir(parents=True, exist_ok=True)
        argv = [executable, "call", "-m" if calling_method == "multiallelic" else "-c"]
        if variants_only:
            argv.append("-v")
        if ploidy != "default":
            argv += ["--ploidy", ploidy]
        argv += ["--threads", str(threads), "-Ov", "-o", str(output)]
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        argv.append(str(source))
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


class BcftoolsFilter:
    CATEGORY = "ComfyBIO/Variants"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING",)
    RETURN_NAMES = ("filtered_vcf",)

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "input_vcf": ("STRING", {"default": ""}),
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

    def run(self, input_vcf, output_vcf: str = "", exclude="QUAL<10", include="", soft_filter="", snp_gap=0, indel_gap=0, threads=1, extra_command=""):
        if include and exclude:
            raise ValueError("include and exclude expressions are mutually exclusive")
        source = _file(input_vcf, "Input VCF/BCF")
        executable = _executable()
        output = Path(output_vcf).expanduser().resolve() if (output_vcf and str(output_vcf).strip()) else (_output_dir("BcftoolsFilter") / f"{Path(input_vcf).stem}.filtered.vcf.gz")
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
        if extra_command.strip():
            argv.extend(shlex.split(extra_command))
        argv.append(str(source))
        _run(argv, output.parent, output)
        _validate(executable, output)
        return (str(output),)


NODE_CLASS_MAPPINGS = {
    "BcftoolsMpileup": BcftoolsMpileup,
    "BcftoolsCall": BcftoolsCall,
    "BcftoolsFilter": BcftoolsFilter,
}
NODE_DISPLAY_NAME_MAPPINGS = {
    "BcftoolsMpileup": "bcftools: Generate Pileup",
    "BcftoolsCall": "bcftools: Call Variants",
    "BcftoolsFilter": "bcftools: Filter Variants",
}


# Backward compatibility aliases
BcftoolsMpileupNode = BcftoolsMpileup
BcftoolsCallNode = BcftoolsCall
BcftoolsFilterNode = BcftoolsFilter

__all__ = [*NODE_CLASS_MAPPINGS, "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
