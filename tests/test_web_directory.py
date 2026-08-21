import json
from pathlib import Path

def test_web_directory_exists_and_valid():
    """Ensure WEB_DIRECTORY points to an existing directory with valid frontend scripts."""
    root_dir = Path(__file__).resolve().parents[1]
    web_dir = root_dir / "web"

    assert web_dir.is_dir(), "web/ directory does not exist"
    assert (web_dir / "js" / "comfybio_main.js").is_file(), "web/js/comfybio_main.js is missing"
    assert (web_dir / "css" / "comfybio.css").is_file(), "web/css/comfybio.css is missing"


def test_example_workflows_valid_json():
    """Ensure all bundled example workflows are valid JSON and have required ComfyUI keys."""
    root_dir = Path(__file__).resolve().parents[1]
    workflows_dir = root_dir / "workflows"

    assert workflows_dir.is_dir(), "workflows/ directory does not exist"
    workflow_files = list(workflows_dir.glob("*.json"))
    assert len(workflow_files) >= 3, f"Expected at least 3 example workflows, found {len(workflow_files)}"

    for wf_path in workflow_files:
        content = json.loads(wf_path.read_text(encoding="utf-8"))
        assert "nodes" in content, f"{wf_path.name} missing 'nodes' key"
        assert "links" in content, f"{wf_path.name} missing 'links' key"
        assert isinstance(content["nodes"], list)
        assert isinstance(content["links"], list)
