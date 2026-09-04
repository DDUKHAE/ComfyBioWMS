import importlib.util
import subprocess
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/bcftools.py").resolve()
BWA_MODULE = Path("nodes/class_2/bwa_mem2.py").resolve()
SAMTOOLS_MODULE = Path("nodes/class_2/samtools.py").resolve()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def bcftools_module():
    return load_module("standalone_bcftools", MODULE)


def test_bcftools_file_owns_three_nodes_and_optional_extra_command(bcftools_module):
    assert set(bcftools_module.NODE_CLASS_MAPPINGS) == {
        "BcftoolsMpileupNode",
        "BcftoolsCallNode",
        "BcftoolsFilterNode",
    }
    for node in bcftools_module.NODE_CLASS_MAPPINGS.values():
        assert "extra_command" in node.INPUT_TYPES()["optional"]
    assert bcftools_module.BcftoolsMpileupNode.INPUT_TYPES()["optional"]["max_depth"][1]["default"] == 250
    assert bcftools_module.BcftoolsFilterNode.INPUT_TYPES()["optional"]["exclude"][1]["default"] == "QUAL<10"


def test_bcftools_filters_managed_short_long_and_equals_options(bcftools_module):
    cases = [
        ("--threads=8 -d500 -a FORMAT/DP", bcftools_module._MPILEUP_MANAGED, ["-a", "FORMAT/DP"]),
        ("-m -v --ploidy 2 --keep-alts", bcftools_module._CALL_MANAGED, ["--keep-alts"]),
        ("--exclude=QUAL<20 -s LowQual --IndelGap 5 --SnpGap=3 -m +x", bcftools_module._FILTER_MANAGED, ["-m", "+x"]),
    ]
    for command, managed, expected in cases:
        kept, ignored = bcftools_module._filter_extra(command, managed)
        assert kept == expected
        assert ignored


def test_bcftools_filter_rejects_include_and_exclude_together(bcftools_module, tmp_path):
    source = tmp_path / "input.vcf"
    source.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(ValueError, match="include and exclude"):
        bcftools_module.BcftoolsFilterNode().run(
            str(source), str(tmp_path / "output.vcf"), include="QUAL>20", exclude="QUAL<10"
        )


def test_bcftools_missing_binary_does_not_create_output(bcftools_module, tmp_path, monkeypatch):
    source = tmp_path / "input.bcf"
    source.write_bytes(b"BCF")
    output = tmp_path / "missing" / "output.vcf"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="bcftools executable not found on PATH"):
        bcftools_module.BcftoolsCallNode().run(str(source), str(output))
    assert not output.parent.exists()


@pytest.mark.e2e
def test_bcftools_mpileup_call_filter_official_sarek_alignment(bcftools_module, tmp_path):
    bwa = load_module("standalone_bwa_for_bcftools", BWA_MODULE)
    samtools = load_module("standalone_samtools_for_bcftools", SAMTOOLS_MODULE)
    reference = fetch_official_data("sarek_ref.fasta", tmp_path)
    read1 = fetch_official_data("sarek_R1.fastq.gz", tmp_path)
    read2 = fetch_official_data("sarek_R2.fastq.gz", tmp_path)
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    (indexed_reference,) = bwa.BwaMem2IndexNode().run(str(reference), str(tmp_path / "bwa"))
    (sam_path,) = bwa.BwaMem2AlignNode().run(
        indexed_reference, str(read1), str(tmp_path / "reads.sam"), read2=str(read2)
    )
    (bam_path,) = samtools.SamtoolsSortNode().run(sam_path, str(tmp_path / "reads.bam"))
    samtools.SamtoolsIndexNode().run(bam_path)

    (bcf_path,) = bcftools_module.BcftoolsMpileupNode().run(
        str(reference), bam_path, str(tmp_path / "calls.bcf"), extra_command="--threads=99 -a FORMAT/DP"
    )
    (called_path,) = bcftools_module.BcftoolsCallNode().run(
        bcf_path, str(tmp_path / "called.vcf"), extra_command="-m --keep-alts"
    )
    (filtered_path,) = bcftools_module.BcftoolsFilterNode().run(
        called_path, str(tmp_path / "filtered.vcf"), exclude="QUAL<10"
    )

    for artifact in (bcf_path, called_path, filtered_path):
        subprocess.run(["bcftools", "view", "-h", artifact], check=True, capture_output=True)
    count = lambda path: int(subprocess.run(
        ["bcftools", "view", "-H", path], check=True, text=True, capture_output=True
    ).stdout.count("\n"))
    assert count(filtered_path) <= count(called_path)
