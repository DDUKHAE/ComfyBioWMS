import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "engine" / "src"))

from bioflow.benchmark.graph_checks import check_dag

_LINEAR = {
    "1": {"class_type": "fastp", "inputs": {"reads": "data/x.fastq"}},
    "2": {"class_type": "salmon", "inputs": {"reads": ["1", 0]}},
    "3": {"class_type": "deseq2", "inputs": {"quant": ["2", 0]}},
}

_CYCLIC = {
    "1": {"class_type": "fastp", "inputs": {"reads": ["3", 0]}},
    "2": {"class_type": "spades", "inputs": {"reads": ["1", 0]}},
    "3": {"class_type": "quast", "inputs": {"contigs": ["2", 0]}},
}

_DANGLING = {
    "1": {"class_type": "fastp", "inputs": {"reads": "data/x.fastq"}},
    "2": {"class_type": "salmon", "inputs": {"reads": ["99", 0]}},
}


def test_linear_pipeline_is_a_dag():
    result = check_dag(_LINEAR)
    assert result["is_dag"] is True
    assert result["has_cycle"] is False
    assert result["dangling_links"] == []
    assert result["node_count"] == 3


def test_cycle_is_detected():
    result = check_dag(_CYCLIC)
    assert result["has_cycle"] is True
    assert result["is_dag"] is False


def test_dangling_link_is_reported():
    result = check_dag(_DANGLING)
    assert result["dangling_links"] == ["2.reads -> 99"]
    assert result["is_dag"] is False


def test_empty_workflow_is_not_a_dag():
    for empty in (None, {}, "not a dict"):
        result = check_dag(empty)
        assert result["is_dag"] is False
        assert result["node_count"] == 0


from bioflow.benchmark.graph_checks import check_port_types

_TOOL_TYPES = {
    "fastp": {"outputs": ["FASTQ"], "inputs": {"reads": "FASTQ"}},
    "bwa_mem2": {"outputs": ["BAM"], "inputs": {"reads": "FASTQ"}},
    "bcftools": {"outputs": ["VCF"], "inputs": {"alignment": "BAM"}},
}


def test_matching_port_types_score_full_fidelity():
    workflow = {
        "1": {"class_type": "fastp", "inputs": {"reads": "data/x.fastq"}},
        "2": {"class_type": "bwa_mem2", "inputs": {"reads": ["1", 0]}},
        "3": {"class_type": "bcftools", "inputs": {"alignment": ["2", 0]}},
    }
    result = check_port_types(workflow, _TOOL_TYPES)
    assert result["mismatches"] == []
    assert result["checked_links"] == 2
    assert result["fidelity_pct"] == 100.0


def test_fastq_into_a_bam_slot_is_a_mismatch():
    workflow = {
        "1": {"class_type": "fastp", "inputs": {"reads": "data/x.fastq"}},
        "2": {"class_type": "bcftools", "inputs": {"alignment": ["1", 0]}},
    }
    result = check_port_types(workflow, _TOOL_TYPES)
    assert result["mismatches"] == ["2.alignment expects BAM but 1 (fastp) emits FASTQ"]
    assert result["fidelity_pct"] == 0.0


def test_unknown_tools_are_not_counted_as_passing():
    workflow = {
        "1": {"class_type": "mystery_tool", "inputs": {}},
        "2": {"class_type": "another_mystery", "inputs": {"x": ["1", 0]}},
    }
    result = check_port_types(workflow, _TOOL_TYPES)
    assert result["checked_links"] == 0
    assert result["fidelity_pct"] == 0.0


from bioflow.benchmark.graph_checks import check_custom_node_source

_VALID_NODE = '''
class SeuratClusteringNode:
    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {"matrix_path": ("STRING", {})}}
    RETURN_TYPES = ("UMAP_PLOT",)
    FUNCTION = "run"
    def run(self, matrix_path):
        return (matrix_path,)
'''

_MISSING_RETURN_TYPES = '''
class HalfBakedNode:
    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {}}
    FUNCTION = "run"
    def run(self):
        return ()
'''


def test_valid_custom_node_is_api_compliant():
    result = check_custom_node_source(_VALID_NODE)
    assert result["parses"] is True
    assert result["class_name"] == "SeuratClusteringNode"
    assert result["api_compliant"] is True


def test_missing_return_types_is_not_compliant():
    result = check_custom_node_source(_MISSING_RETURN_TYPES)
    assert result["parses"] is True
    assert result["has_return_types"] is False
    assert result["api_compliant"] is False


def test_syntax_error_is_not_compliant():
    result = check_custom_node_source("class Broken(:\n  pass")
    assert result["parses"] is False
    assert result["api_compliant"] is False


def test_absent_source_is_not_compliant():
    result = check_custom_node_source(None)
    assert result["api_compliant"] is False


