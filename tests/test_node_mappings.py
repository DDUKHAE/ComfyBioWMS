import inspect


def test_root_init_exports_exact_verified_set():
    import __init__ as root_package

    assert root_package.WEB_DIRECTORY == "./web"
    assert len(root_package.NODE_CLASS_MAPPINGS) == 90
    assert set(root_package.NODE_CLASS_MAPPINGS) == set(
        root_package.NODE_DISPLAY_NAME_MAPPINGS
    )


def test_all_registered_nodes_comply_with_comfyui_spec():
    from nodes.registry import NODE_CLASS_MAPPINGS

    for name, node in NODE_CLASS_MAPPINGS.items():
        assert inspect.isclass(node), name
        assert "required" in node.INPUT_TYPES()
        assert isinstance(node.RETURN_TYPES, tuple)
        assert callable(getattr(node, node.FUNCTION))
        assert isinstance(node.CATEGORY, str)
