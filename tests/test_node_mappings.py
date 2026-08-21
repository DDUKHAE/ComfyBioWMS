import inspect
from pathlib import Path

def test_root_init_exports():
    """Verify that root __init__.py exports the required ComfyUI mappings and web directory."""
    import sys
    root_dir = Path(__file__).resolve().parents[1]
    if str(root_dir) not in sys.path:
        sys.path.insert(0, str(root_dir))

    import __init__ as root_pkg

    assert hasattr(root_pkg, "NODE_CLASS_MAPPINGS"), "Missing NODE_CLASS_MAPPINGS in root __init__.py"
    assert hasattr(root_pkg, "NODE_DISPLAY_NAME_MAPPINGS"), "Missing NODE_DISPLAY_NAME_MAPPINGS in root __init__.py"
    assert hasattr(root_pkg, "WEB_DIRECTORY"), "Missing WEB_DIRECTORY in root __init__.py"
    assert isinstance(root_pkg.NODE_CLASS_MAPPINGS, dict)
    assert isinstance(root_pkg.NODE_DISPLAY_NAME_MAPPINGS, dict)
    assert len(root_pkg.NODE_CLASS_MAPPINGS) > 0, "No nodes registered in NODE_CLASS_MAPPINGS"


def test_all_registered_nodes_comply_with_comfyui_spec():
    """Verify that all 170+ registered node classes meet ComfyUI custom node requirements."""
    from nodes.registry import NODE_CLASS_MAPPINGS, NODE_DISPLAY_NAME_MAPPINGS

    assert len(NODE_CLASS_MAPPINGS) >= 150, f"Expected 150+ nodes, found {len(NODE_CLASS_MAPPINGS)}"
    assert len(NODE_DISPLAY_NAME_MAPPINGS) == len(NODE_CLASS_MAPPINGS)

    for node_name, node_cls in NODE_CLASS_MAPPINGS.items():
        assert inspect.isclass(node_cls), f"{node_name} is not a class"

        # Check INPUT_TYPES
        assert hasattr(node_cls, "INPUT_TYPES"), f"{node_name} is missing INPUT_TYPES"
        input_types = node_cls.INPUT_TYPES() if callable(node_cls.INPUT_TYPES) else node_cls.INPUT_TYPES
        assert isinstance(input_types, dict), f"{node_name}.INPUT_TYPES must return a dict"
        assert "required" in input_types, f"{node_name}.INPUT_TYPES missing 'required' key"

        # Check RETURN_TYPES
        assert hasattr(node_cls, "RETURN_TYPES"), f"{node_name} is missing RETURN_TYPES"
        assert isinstance(node_cls.RETURN_TYPES, tuple), f"{node_name}.RETURN_TYPES must be a tuple"

        # Check FUNCTION
        assert hasattr(node_cls, "FUNCTION"), f"{node_name} is missing FUNCTION attribute"
        func_name = node_cls.FUNCTION
        assert isinstance(func_name, str), f"{node_name}.FUNCTION must be a string"
        assert hasattr(node_cls, func_name), f"{node_name} has no method named '{func_name}'"
        assert callable(getattr(node_cls, func_name)), f"{node_name}.{func_name} is not callable"

        # Check CATEGORY
        assert hasattr(node_cls, "CATEGORY"), f"{node_name} is missing CATEGORY attribute"
        assert isinstance(node_cls.CATEGORY, str), f"{node_name}.CATEGORY must be a string"
