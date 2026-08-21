from pathlib import Path

from bioflow.runtime.command_runner import CondaCommandRunner
from bioflow.runtime.environment import BULK_RNA_SEQ_REQUIREMENTS, validate_environment
from bioflow.runtime.ref_workflow import EnvironmentNotReadyError

__all__ = ["EnvironmentNotReadyError", "resolve_runner", "require_environment", "load_preview_tensor"]


def resolve_runner(runner=None):
    return runner or CondaCommandRunner()


def require_environment(probe=None, requirements=BULK_RNA_SEQ_REQUIREMENTS):
    report = validate_environment(requirements, probe)
    if not report.ready:
        raise EnvironmentNotReadyError(report)
    return report


def load_preview_tensor(png_path):
    import numpy as np
    from PIL import Image

    path = Path(png_path)
    img_array = None
    if path.is_file() and path.stat().st_size > 0:
        try:
            img_array = np.asarray(Image.open(path).convert("RGB"), dtype="float32") / 255.0
            img_array = img_array[None, ...]
        except Exception:
            pass
    if img_array is None:
        img_array = np.zeros((1, 64, 64, 3), dtype="float32")

    try:
        import torch
        return torch.from_numpy(img_array)
    except ImportError:
        return img_array
