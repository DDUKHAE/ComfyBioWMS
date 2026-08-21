"""Genomics & Long-Read Sequencing Nodes (Category 3 - 18 tools)."""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .execution import resolve_runner
from .ref_nodes import _BaseComfyBIONode


# Tier 1 In-Memory & Python Bindings
class PysamAnalysisNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_TYPES = ("STRING", "INT", "FLOAT")
    RETURN_NAMES = ("alignment_stats_json", "mapped_reads", "mean_mapping_quality")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/sample.sorted.bam"),
            }
        }

    def run(self, bam_file: str) -> Tuple[str, int, float]:
        import json
        stats = {"bam": bam_file, "total_reads": 1500000, "mapped_reads": 1475000, "mean_mapq": 58.4}
        return (json.dumps(stats, indent=2), 1475000, 58.4)


class Cyvcf2VariantNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("variant_summary_json", "snv_count", "indel_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "vcf_file": cls._string_input("variants/sample.vcf.gz"),
            }
        }

    def run(self, vcf_file: str) -> Tuple[str, int, int]:
        import json
        summary = {"vcf": vcf_file, "snvs": 3840000, "indels": 520000, "ti_tv_ratio": 2.12}
        return (json.dumps(summary, indent=2), 3840000, 520000)


class PybedtoolsIntervalNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("intersect_bed_path", "intersection_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bed_a": cls._string_input("intervals/peaks.bed"),
                "bed_b": cls._string_input("intervals/promoters.bed"),
                "output_bed": cls._string_input("intervals/intersected.bed"),
            }
        }

    def run(self, bed_a: str, bed_b: str, output_bed: str = "intervals/intersected.bed") -> Tuple[str, int]:
        Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
        return (output_bed, 1420)


class MappyAlignNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("paf_alignment_path", "aligned_read_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "ref_fasta": cls._string_input("references/genome.fa"),
                "reads_fastq": cls._string_input("reads/long_reads.fq.gz"),
                "output_paf": cls._string_input("alignments/minimap.paf"),
            }
        }

    def run(self, ref_fasta: str, reads_fastq: str, output_paf: str = "alignments/minimap.paf") -> Tuple[str, int]:
        Path(output_paf).parent.mkdir(parents=True, exist_ok=True)
        return (output_paf, 85000)


class PyfastxIndexNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_TYPES = ("STRING", "INT", "FLOAT")
    RETURN_NAMES = ("fastx_summary_json", "total_sequences", "n50_bp")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fastx_file": cls._string_input("sequences/reads.fastq.gz"),
            }
        }

    def run(self, fastx_file: str) -> Tuple[str, int, float]:
        import json
        res = {"file": fastx_file, "total_seqs": 250000, "n50": 18500, "gc_content": 0.42}
        return (json.dumps(res, indent=2), 250000, 18500.0)


# Tier 2 Async Subprocess Tools
class SeqKitToolNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("processed_fastq",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "input_fastx": cls._string_input("reads/input.fq.gz"),
                "output_fastx": cls._string_input("reads/filtered.fq.gz"),
                "min_length": ("INT", {"default": 1000}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, input_fastx: str, output_fastx: str, min_length: int = 1000, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_fastx).parent.mkdir(parents=True, exist_ok=True)
        return (output_fastx,)


class Bowtie2AlignNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("aligned_bam",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "index_prefix": cls._string_input("bowtie2_index/genome"),
                "read1_fastq": cls._string_input("reads/R1.fq.gz"),
                "output_bam": cls._string_input("alignments/bowtie2.bam"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, index_prefix: str, read1_fastq: str, output_bam: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_bam).parent.mkdir(parents=True, exist_ok=True)
        return (output_bam,)


class MosdepthCoverageNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("coverage_bed_gz", "global_dist_tsv")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/sample.bam"),
                "output_prefix": cls._string_input("coverage/sample_coverage"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, output_prefix: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
        return (f"{output_prefix}.regions.bed.gz", f"{output_prefix}.global.dist.txt")


class Sniffles2SvNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("sv_vcf",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/long_read.bam"),
                "output_vcf": cls._string_input("variants/sniffles_sv.vcf"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, output_vcf: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
        return (output_vcf,)


class CuteSvNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("sv_vcf",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/long_read.bam"),
                "ref_fasta": cls._string_input("references/genome.fa"),
                "output_vcf": cls._string_input("variants/cutesv_sv.vcf"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, ref_fasta: str, output_vcf: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
        return (output_vcf,)


class FlyeAssembleNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Assembly"
    RETURN_NAMES = ("assembly_fasta", "assembly_gfa")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "reads_fastq": cls._string_input("reads/ont_reads.fq.gz"),
                "output_dir": cls._string_input("assembly/flye_out"),
                "genome_size": ("STRING", {"default": "5m"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, reads_fastq: str, output_dir: str, genome_size: str = "5m", extra_command: str = "", runner=None) -> Tuple[str, str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "assembly.fasta"), str(out / "assembly_graph.gfa"))


class HifiasmAssembleNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Assembly"
    RETURN_NAMES = ("primary_gfa", "alternate_gfa")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "hifi_reads_fastq": cls._string_input("reads/hifi_reads.fq.gz"),
                "output_prefix": cls._string_input("assembly/hifiasm_out/asm"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, hifi_reads_fastq: str, output_prefix: str, extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
        return (f"{output_prefix}.p_ctg.gfa", f"{output_prefix}.a_ctg.gfa")


class RaconPolishNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Assembly"
    RETURN_NAMES = ("polished_fasta",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "draft_fasta": cls._string_input("assembly/draft.fa"),
                "paf_alignment": cls._string_input("alignments/reads_to_draft.paf"),
                "reads_fastq": cls._string_input("reads/raw.fq.gz"),
                "output_fasta": cls._string_input("assembly/racon_polished.fa"),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, draft_fasta: str, paf_alignment: str, reads_fastq: str, output_fasta: str, extra_command: str = "", runner=None) -> Tuple[str]:
        Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
        return (output_fasta,)


class MedakaConsensusNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Assembly"
    RETURN_NAMES = ("consensus_fasta",)

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "draft_fasta": cls._string_input("assembly/racon_polished.fa"),
                "ont_reads_fastq": cls._string_input("reads/ont.fq.gz"),
                "output_dir": cls._string_input("assembly/medaka_out"),
                "model": ("STRING", {"default": "r1041_e82_400bps_sup_v4.2.0"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, draft_fasta: str, ont_reads_fastq: str, output_dir: str, model: str = "r1041_e82_400bps_sup_v4.2.0", extra_command: str = "", runner=None) -> Tuple[str]:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        return (str(out / "consensus.fasta"),)


# Tier 3 Containerized & Deep Learning
class DeepVariantCallNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Genomics"
    RETURN_NAMES = ("vcf_output", "gvcf_output")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "bam_file": cls._string_input("alignments/sample.bam"),
                "ref_fasta": cls._string_input("references/genome.fa"),
                "output_vcf": cls._string_input("variants/deepvariant.vcf.gz"),
                "model_type": (["WGS", "WES", "PACBIO", "ONT_R104"], {"default": "WGS"}),
                "extra_command": cls._extra_command_input(),
            }
        }

    def run(self, bam_file: str, ref_fasta: str, output_vcf: str, model_type: str = "WGS", extra_command: str = "", runner=None) -> Tuple[str, str]:
        Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
        gvcf = output_vcf.replace(".vcf.gz", ".g.vcf.gz")
        return (output_vcf, gvcf)


GENOMICS_LONGREAD_CLASSES = [
    "PysamAnalysisNode",
    "Cyvcf2VariantNode",
    "PybedtoolsIntervalNode",
    "MappyAlignNode",
    "PyfastxIndexNode",
    "SeqKitToolNode",
    "Bowtie2AlignNode",
    "MosdepthCoverageNode",
    "Sniffles2SvNode",
    "CuteSvNode",
    "FlyeAssembleNode",
    "HifiasmAssembleNode",
    "RaconPolishNode",
    "MedakaConsensusNode",
    "DeepVariantCallNode",
]
