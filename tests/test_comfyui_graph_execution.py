"""
ComfyUI Runtime Graph Execution, Caching, and Fault-Propagation Test Suite.
Verifies the execution matrix and cache invalidation boundaries mandated by Section 6 (T3) of docs/review-remediation-guideline.md.
"""

import copy
import os
import shutil
import tempfile
import time
from pathlib import Path
import pytest

import __init__ as comfy_entrypoint
from bioflow.benchmark.graph_checks import check_dag


class SimpleComfyExecutor:
    """
    Simulates the ComfyUI execution engine lifecycle:
    1. Dependency graph resolution from output targets.
    2. IS_CHANGED-based input/parameter caching.
    3. Error handling: node failure immediately halts execution and marks dependent downstream nodes as 'blocked'.
    4. Records execution status for each node: 'executed', 'cached', 'failed', 'blocked', 'not_evaluated'.
    """

    def __init__(self, node_classes: dict):
        self.node_classes = node_classes
        self.cache_store = {}  # node_id -> (cache_key, outputs)

    def execute_prompt(self, prompt: dict, output_nodes: list[str]) -> dict:
        """
        Executes a ComfyUI prompt dictionary.
        Returns execution results, node statuses, and logs.
        """
        # 1. Structural validation
        dag_res = check_dag(prompt)
        if not dag_res["is_dag"]:
            raise ValueError(f"Graph validation failed: {dag_res}")

        # 2. Identify required active nodes via reverse DFS from output_nodes
        required_nodes = set()
        stack = list(output_nodes)
        while stack:
            nid = str(stack.pop())
            if nid not in prompt:
                raise KeyError(f"Terminal node {nid} not found in prompt")
            required_nodes.add(nid)
            node_data = prompt[nid]
            for input_val in (node_data.get("inputs") or {}).values():
                if isinstance(input_val, list) and len(input_val) == 2 and isinstance(input_val[0], (str, int)):
                    parent_id = str(input_val[0])
                    if parent_id not in required_nodes:
                        stack.append(parent_id)

        # 3. Topological sort of required nodes
        # in-degree of required subgraph
        in_degree = {nid: 0 for nid in required_nodes}
        adj = {nid: [] for nid in required_nodes}
        for nid in required_nodes:
            for input_val in (prompt[nid].get("inputs") or {}).values():
                if isinstance(input_val, list) and len(input_val) == 2:
                    p_id = str(input_val[0])
                    if p_id in required_nodes:
                        adj[p_id].append(nid)
                        in_degree[nid] += 1

        zero_in = [nid for nid, deg in in_degree.items() if deg == 0]
        topo_order = []
        while zero_in:
            curr = zero_in.pop(0)
            topo_order.append(curr)
            for child in adj[curr]:
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    zero_in.append(child)

        node_status = {nid: "not_evaluated" for nid in prompt}
        for nid in required_nodes:
            node_status[nid] = "pending"

        node_outputs = {}
        execution_log = []

        # 4. Execute in topological order
        for nid in topo_order:
            if node_status[nid] == "blocked":
                continue

            node_def = prompt[nid]
            class_type = node_def["class_type"]
            node_cls = self.node_classes.get(class_type)
            if not node_cls:
                node_status[nid] = "failed"
                execution_log.append(f"Node {nid} failed: class '{class_type}' not found in registry")
                self._block_descendants(nid, adj, node_status)
                continue

            # Resolve input values
            resolved_inputs = {}
            input_spec = node_def.get("inputs", {})
            has_blocked_input = False
            for k, v in input_spec.items():
                if isinstance(v, list) and len(v) == 2 and isinstance(v[0], (str, int)):
                    src_id = str(v[0])
                    src_slot = v[1]
                    if node_status.get(src_id) in ["failed", "blocked"]:
                        has_blocked_input = True
                        break
                    src_outs = node_outputs.get(src_id)
                    if src_outs is not None:
                        if isinstance(src_outs, (tuple, list)):
                            resolved_inputs[k] = src_outs[src_slot] if src_slot < len(src_outs) else src_outs[0]
                        else:
                            resolved_inputs[k] = src_outs
                else:
                    resolved_inputs[k] = v

            if has_blocked_input:
                node_status[nid] = "blocked"
                self._block_descendants(nid, adj, node_status)
                continue

            # Compute cache key (simulating IS_CHANGED)
            cache_key = None
            if hasattr(node_cls, "IS_CHANGED"):
                try:
                    # Filter inputs matching IS_CHANGED parameters
                    cache_key = node_cls.IS_CHANGED(**{k: v for k, v in resolved_inputs.items() if not isinstance(v, Path)})
                except Exception:
                    cache_key = None
            if cache_key is None:
                # Fallback to tuple of inputs
                cache_key = str(sorted(resolved_inputs.items()))

            # Check cache
            cached_entry = self.cache_store.get(nid)
            if cached_entry and cached_entry[0] == cache_key:
                node_status[nid] = "cached"
                node_outputs[nid] = cached_entry[1]
                execution_log.append(f"Node {nid} ({class_type}): CACHED")
                continue

            # Execute node
            func_name = getattr(node_cls, "FUNCTION", "run")
            inst = node_cls()
            func = getattr(inst, func_name)

            try:
                # Only pass args accepted by func or required
                outs = func(**resolved_inputs)
                node_outputs[nid] = outs
                node_status[nid] = "executed"
                self.cache_store[nid] = (cache_key, outs)
                execution_log.append(f"Node {nid} ({class_type}): EXECUTED successfully")
            except Exception as e:
                node_status[nid] = "failed"
                execution_log.append(f"Node {nid} ({class_type}): FAILED with {type(e).__name__}: {e}")
                self._block_descendants(nid, adj, node_status)

        return {
            "node_status": node_status,
            "node_outputs": node_outputs,
            "execution_log": execution_log,
        }

    def _block_descendants(self, failed_nid: str, adj: dict, node_status: dict):
        stack = list(adj.get(failed_nid, []))
        while stack:
            curr = stack.pop()
            if node_status[curr] != "failed":
                node_status[curr] = "blocked"
            stack.extend(adj.get(curr, []))


# -----------------------------------------------------------------------------
# Test Fixtures & Classes
# -----------------------------------------------------------------------------

class MockFastpNode:
    FUNCTION = "run"
    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {"read1": ("STRING",), "threads": ("INT", {"default": 2})}}
    @classmethod
    def IS_CHANGED(cls, read1, threads=2):
        p = Path(read1)
        mtime = p.stat().st_mtime_ns if p.exists() else 0
        size = p.stat().st_size if p.exists() else 0
        return f"{p}:{mtime}:{size}:{threads}"
    def run(self, read1, threads=2):
        p = Path(read1)
        if not p.exists():
            raise FileNotFoundError(f"Input file not found: {read1}")
        out_f = p.parent / f"trimmed_{p.name}"
        out_f.write_text(f"TRIMMED({threads}):\n" + p.read_text())
        return (str(out_f),)


class MockTximportNode:
    FUNCTION = "run"
    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {"quant_file": ("STRING",), "scaling": ("STRING", {"default": "lengthScaledTPM"})}}
    @classmethod
    def IS_CHANGED(cls, quant_file, scaling="lengthScaledTPM"):
        p = Path(quant_file)
        mtime = p.stat().st_mtime_ns if p.exists() else 0
        size = p.stat().st_size if p.exists() else 0
        return f"{p}:{mtime}:{size}:{scaling}"
    def run(self, quant_file, scaling="lengthScaledTPM"):
        p = Path(quant_file)
        if not p.exists():
            raise FileNotFoundError(f"Quant file not found: {quant_file}")
        out_mat = p.parent / "counts_matrix.tsv"
        out_mat.write_text(f"gene\tcount\nGENE_1\t100\nGENE_2\t200\n# scaling={scaling}\n")
        return (str(out_mat),)


class MockPlotNode:
    FUNCTION = "run"
    @classmethod
    def INPUT_TYPES(cls):
        return {"required": {"count_table": ("STRING",), "cutoff": ("FLOAT", {"default": 0.05})}}
    @classmethod
    def IS_CHANGED(cls, count_table, cutoff=0.05):
        p = Path(count_table)
        mtime = p.stat().st_mtime_ns if p.exists() else 0
        size = p.stat().st_size if p.exists() else 0
        return f"{p}:{mtime}:{size}:{cutoff}"
    def run(self, count_table, cutoff=0.05):
        p = Path(count_table)
        if not p.exists():
            raise FileNotFoundError(f"Table file not found: {count_table}")
        out_plot = p.parent / "volcano.png"
        out_plot.write_bytes(b"FAKE_PNG_BYTES")
        return (str(out_plot),)


MOCK_REGISTRY = {
    "Fastp": MockFastpNode,
    "Tximport": MockTximportNode,
    "VolcanoPlot": MockPlotNode,
}


@pytest.fixture
def sample_pipeline(tmp_path):
    raw_fq = tmp_path / "sample.fastq"
    raw_fq.write_text("@READ_1\nACTG\n+\nIIII\n")

    prompt = {
        "1": {"class_type": "Fastp", "inputs": {"read1": str(raw_fq), "threads": 2}},
        "2": {"class_type": "Tximport", "inputs": {"quant_file": ["1", 0], "scaling": "lengthScaledTPM"}},
        "3": {"class_type": "VolcanoPlot", "inputs": {"count_table": ["2", 0], "cutoff": 0.05}},
        "4": {"class_type": "Fastp", "inputs": {"read1": str(raw_fq), "threads": 4}},  # Disconnected / independent branch
    }
    return prompt, raw_fq, tmp_path


def test_execution_01_normal_run(sample_pipeline):
    """Normal execution of linear sub-graph to target output node '3'."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    res = executor.execute_prompt(prompt, output_nodes=["3"])
    status = res["node_status"]

    assert status["1"] == "executed"
    assert status["2"] == "executed"
    assert status["3"] == "executed"
    assert status["4"] == "not_evaluated", "Disconnected node not in target subgraph must remain not_evaluated"


def test_execution_02_resubmit_identical_cache_hit(sample_pipeline):
    """Resubmission under identical conditions -> upstream and downstream nodes must be CACHED."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    # 1st run
    res1 = executor.execute_prompt(prompt, output_nodes=["3"])
    assert res1["node_status"]["1"] == "executed"

    # 2nd run
    res2 = executor.execute_prompt(prompt, output_nodes=["3"])
    status2 = res2["node_status"]
    assert status2["1"] == "cached"
    assert status2["2"] == "cached"
    assert status2["3"] == "cached"


def test_execution_03_downstream_option_change(sample_pipeline):
    """Modifying downstream option (cutoff 0.05 -> 0.01) must reuse upstream cache and re-execute downstream."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    # 1st run
    executor.execute_prompt(prompt, output_nodes=["3"])

    # Modify downstream cutoff only
    p2 = copy.deepcopy(prompt)
    p2["3"]["inputs"]["cutoff"] = 0.01

    res2 = executor.execute_prompt(p2, output_nodes=["3"])
    status2 = res2["node_status"]
    assert status2["1"] == "cached", "Upstream node 1 must be cached"
    assert status2["2"] == "cached", "Upstream node 2 must be cached"
    assert status2["3"] == "executed", "Downstream node 3 must be re-executed"


def test_execution_04_upstream_option_change(sample_pipeline):
    """Modifying upstream option (threads 2 -> 8) must re-execute upstream and invalidate downstream."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    # 1st run
    executor.execute_prompt(prompt, output_nodes=["3"])

    # Modify upstream parameter
    p2 = copy.deepcopy(prompt)
    p2["1"]["inputs"]["threads"] = 8

    res2 = executor.execute_prompt(p2, output_nodes=["3"])
    status2 = res2["node_status"]
    assert status2["1"] == "executed", "Upstream node 1 must re-execute"
    assert status2["2"] == "executed", "Downstream node 2 must re-execute"
    assert status2["3"] == "executed", "Downstream node 3 must re-execute"


def test_execution_05_input_file_content_modification(sample_pipeline):
    """Modifying file content (changing mtime and size) must trigger cache invalidation."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    # 1st run
    executor.execute_prompt(prompt, output_nodes=["3"])

    # Sleep slightly to guarantee different mtime
    time.sleep(0.01)
    raw_fq.write_text("@READ_1\nACTG\n+\nIIII\n@READ_2\nGGCC\n+\nIIII\n")

    res2 = executor.execute_prompt(prompt, output_nodes=["3"])
    status2 = res2["node_status"]
    assert status2["1"] == "executed", "Input change must invalidate upstream node cache"
    assert status2["2"] == "executed"
    assert status2["3"] == "executed"


def test_execution_06_fault_missing_input_blocks_downstream(sample_pipeline):
    """Missing input file causes node 1 to FAIL and marks downstream nodes as BLOCKED."""
    prompt, raw_fq, _ = sample_pipeline
    executor = SimpleComfyExecutor(MOCK_REGISTRY)

    p_fail = copy.deepcopy(prompt)
    p_fail["1"]["inputs"]["read1"] = "/nonexistent/path/missing.fastq"

    res = executor.execute_prompt(p_fail, output_nodes=["3"])
    status = res["node_status"]

    assert status["1"] == "failed"
    assert status["2"] == "blocked", "Downstream node 2 must be blocked due to upstream failure"
    assert status["3"] == "blocked", "Downstream node 3 must be blocked due to upstream failure"


def test_execution_07_invalid_graph_rejected():
    """Graph with cycle must be rejected during graph validation before execution."""
    cyclic_prompt = {
        "1": {"class_type": "Fastp", "inputs": {"read1": ["2", 0]}},
        "2": {"class_type": "Tximport", "inputs": {"quant_file": ["1", 0]}},
    }
    executor = SimpleComfyExecutor(MOCK_REGISTRY)
    with pytest.raises(ValueError, match="Graph validation failed"):
        executor.execute_prompt(cyclic_prompt, output_nodes=["1"])


def test_execution_08_cache_boundary_audit_mtime_size_limitation(tmp_path):
    """
    Audits the boundary of the 'mtime:size' caching strategy.
    Demonstrates that same-size same-mtime content manipulation is not detected by mtime:size alone,
    verifying why the manuscript must document this as an explicit metadata-based cache boundary.
    """
    test_f = tmp_path / "test.txt"
    test_f.write_text("ABCD")  # 4 bytes
    stat1 = test_f.stat()

    key1 = MockFastpNode.IS_CHANGED(str(test_f))

    # Overwrite with same byte length and preserve mtime
    test_f.write_text("WXYZ")  # 4 bytes
    os.utime(test_f, ns=(stat1.st_atime_ns, stat1.st_mtime_ns))

    key2 = MockFastpNode.IS_CHANGED(str(test_f))

    # IS_CHANGED produces identical key because size and mtime are identical!
    assert key1 == key2, "Demonstrates boundary: mtime:size does not detect byte-level mutation if size and mtime are preserved"
