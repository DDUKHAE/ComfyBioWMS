import importlib.util
import subprocess
from pathlib import Path

import pytest

from tests.official_data import fetch_official_data


MODULE = Path("nodes/class_2/samtools.py").resolve()
BWA_MODULE = Path("nodes/class_2/bwa_mem2.py").resolve()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def samtools_module():
    return load_module("standalone_samtools", MODULE)


def test_samtools_file_owns_three_nodes_and_optional_extra_command(samtools_module):
    assert set(samtools_module.NODE_CLASS_MAPPINGS) == {
        "SamtoolsSortNode",
        "SamtoolsIndexNode",
        "SamtoolsMarkdupNode",
    }
    for node in samtools_module.NODE_CLASS_MAPPINGS.values():
        assert "extra_command" in node.INPUT_TYPES()["optional"]


def test_samtools_filters_node_managed_options(samtools_module):
    cases = [
        ("-@8 -m=2G -n --write-index", samtools_module._SORT_MANAGED, ["--write-index"]),
        ("-@ 8 -c -o other.csi --verbose", samtools_module._INDEX_MANAGED, ["--verbose"]),
        ("-@8 -r --mode sequence -d 10 -f stats.txt --no-multi-dup",
         samtools_module._MARKDUP_MANAGED, ["--no-multi-dup"]),
    ]
    for command, managed, expected in cases:
        kept, ignored = samtools_module._filter_extra(command, managed)
        assert kept == expected
        assert ignored


def test_samtools_missing_binary_does_not_create_output(samtools_module, tmp_path, monkeypatch):
    source = tmp_path / "input.sam"
    source.write_text("@HD\tVN:1.6\n")
    output = tmp_path / "missing" / "output.bam"
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="samtools executable not found on PATH"):
        samtools_module.SamtoolsSortNode().run(str(source), str(output))
    assert not output.parent.exists()


@pytest.mark.e2e
def test_samtools_sort_index_markdup_official_sarek_alignment(samtools_module, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    bwa = load_module("standalone_bwa_for_samtools", BWA_MODULE)
    reference = fetch_official_data("sarek_ref.fasta", tmp_path)
    read1 = fetch_official_data("sarek_R1.fastq.gz", tmp_path)
    read2 = fetch_official_data("sarek_R2.fastq.gz", tmp_path)
    (indexed_reference,) = bwa.BwaMem2IndexNode().run(str(reference), "bwa")
    (sam_path,) = bwa.BwaMem2AlignNode().run(
        indexed_reference, str(read1), "reads.sam", read2=str(read2)
    )

    (sorted_bam,) = samtools_module.SamtoolsSortNode().run(
        sam_path, "sorted.bam", threads=2, extra_command="-@ 99"
    )
    (index_path,) = samtools_module.SamtoolsIndexNode().run(sorted_bam, threads=2)
    (marked_bam,) = samtools_module.SamtoolsMarkdupNode().run(
        sorted_bam, "marked.bam", threads=2
    )

    assert Path(sorted_bam).is_absolute()
    assert Path(index_path).is_absolute()
    assert Path(marked_bam).is_absolute()
    assert Path(index_path).is_file()
    subprocess.run(["samtools", "quickcheck", sorted_bam, marked_bam], check=True)
    idxstats = subprocess.run(
        ["samtools", "idxstats", sorted_bam], check=True, text=True, capture_output=True
    ).stdout
    assert idxstats.strip()
