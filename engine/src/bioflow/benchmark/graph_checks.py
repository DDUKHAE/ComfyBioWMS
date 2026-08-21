from __future__ import annotations

import ast
from typing import Any


def _links(workflow: dict) -> list[tuple[str, str, str]]:
    """Yields (source_node_id, target_node_id, input_name) for every [node_id, slot] link."""
    out = []
    for node_id, node in workflow.items():
        if not isinstance(node, dict):
            continue
        for input_name, value in (node.get("inputs") or {}).items():
            if isinstance(value, list) and len(value) == 2 and isinstance(value[0], (str, int)):
                out.append((str(value[0]), str(node_id), str(input_name)))
    return out


def check_dag(workflow: Any) -> dict:
    """Structural check of a ComfyUI API-format workflow: acyclicity and link resolution."""
    if not isinstance(workflow, dict) or not workflow:
        return {"is_dag": False, "has_cycle": False, "dangling_links": [], "node_count": 0}

    node_ids = {str(node_id) for node_id in workflow}
    edges = _links(workflow)
    dangling = [f"{target}.{name} -> {source}" for source, target, name in edges if source not in node_ids]

    adjacency: dict[str, list[str]] = {node_id: [] for node_id in node_ids}
    for source, target, _name in edges:
        if source in node_ids:
            adjacency[source].append(target)

    # Iterative DFS with colouring: 0 unvisited, 1 on stack, 2 done. Avoids recursion limits on
    # long pipelines and reports a cycle the moment a grey node is re-entered.
    colour = dict.fromkeys(node_ids, 0)
    has_cycle = False
    for root in sorted(node_ids):
        if colour[root] != 0:
            continue
        stack = [(root, iter(adjacency[root]))]
        colour[root] = 1
        while stack:
            node, children = stack[-1]
            nxt = next(children, None)
            if nxt is None:
                colour[node] = 2
                stack.pop()
                continue
            if colour.get(nxt) == 1:
                has_cycle = True
                stack.clear()
                break
            if colour.get(nxt) == 0:
                colour[nxt] = 1
                stack.append((nxt, iter(adjacency[nxt])))
        if has_cycle:
            break

    return {"is_dag": not has_cycle and not dangling, "has_cycle": has_cycle, "dangling_links": dangling, "node_count": len(node_ids)}


def check_port_types(workflow: Any, tool_types: dict[str, dict]) -> dict:
    """Compares each [node_id, slot] link against the registry's declared port types.

    Links whose endpoints are not in `tool_types` are excluded from the denominator rather than
    counted as passing -- an unknown tool is unverified, not verified-correct.
    """
    if not isinstance(workflow, dict) or not workflow:
        return {"checked_links": 0, "mismatches": [], "fidelity_pct": 0.0}

    classes = {str(node_id): str((node or {}).get("class_type", "")).lower() for node_id, node in workflow.items() if isinstance(node, dict)}
    checked, mismatches = 0, []
    for source, target, input_name in _links(workflow):
        source_spec = tool_types.get(classes.get(source, ""))
        target_spec = tool_types.get(classes.get(target, ""))
        if not source_spec or not target_spec:
            continue
        expected = (target_spec.get("inputs") or {}).get(input_name)
        emitted = (source_spec.get("outputs") or [None])[0]
        if expected is None or emitted is None:
            continue
        checked += 1
        if expected != emitted:
            mismatches.append(f"{target}.{input_name} expects {expected} but {source} ({classes[source]}) emits {emitted}")

    fidelity = ((checked - len(mismatches)) / checked * 100.0) if checked else 0.0
    return {"checked_links": checked, "mismatches": mismatches, "fidelity_pct": fidelity}


def load_tool_types(catalog_path: Path | str | None = None) -> dict[str, dict]:
    """Extracts {class_type: {'inputs': {name: type}, 'outputs': [type, ...]}} from catalog or registered nodes."""
    import json
    from pathlib import Path
    if catalog_path is None:
        return {}
    path = Path(catalog_path)
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, dict] = {}
    if isinstance(data, list):
        for item in data:
            if not isinstance(item, dict):
                continue
            name = str(item.get("node_type", "") or item.get("title", "")).lower()
            if not name:
                continue
            inputs = {inp["name"]: inp.get("type", "ANY") for inp in item.get("inputs", []) if isinstance(inp, dict) and "name" in inp}
            outputs = [outp.get("type", "ANY") for outp in item.get("outputs", []) if isinstance(outp, dict)]
            out[name] = {"inputs": inputs, "outputs": outputs}
            clean = name.replace("node", "")
            if clean:
                out[clean] = {"inputs": inputs, "outputs": outputs}
    return out


def check_custom_node_source(source: Any) -> dict:
    blank = {"parses": False, "class_name": None, "has_input_types": False, "has_return_types": False, "has_function": False, "api_compliant": bool(False)}
    if not isinstance(source, str) or not source.strip():
        return blank
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return blank

    best = {**blank, "parses": True}
    for node in ast.walk(tree):
        if not isinstance(node, ast.ClassDef):
            continue
        names = set()
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                names.add(item.name)
            elif isinstance(item, ast.Assign):
                names.update(target.id for target in item.targets if isinstance(target, ast.Name))
        report = {
            "parses": True,
            "class_name": node.name,
            "has_input_types": "INPUT_TYPES" in names,
            "has_return_types": "RETURN_TYPES" in names,
            "has_function": "FUNCTION" in names,
        }
        report["api_compliant"] = all((report["has_input_types"], report["has_return_types"], report["has_function"]))
        if report["api_compliant"]:
            return report
        best = report
    return best

