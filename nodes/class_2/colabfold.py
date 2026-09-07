"""ColabFold AlphaFold2 protein structure prediction node.

Python packages: colabfold
External binaries: colabfold_batch
"""

import shlex
import shutil
import subprocess
from pathlib import Path


def _file(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")
    return path

def _output_dir(node_name: str, custom_dir: str = "") -> Path:
    if custom_dir and str(custom_dir).strip():
        out = Path(custom_dir).expanduser().resolve()
    else:
        try:
            folder_paths = __import__("folder_paths")
            base = Path(folder_paths.get_output_directory())
        except Exception:
            base = Path.cwd() / "ComfyUI" / "output"
        out = base / node_name
    return out



class ColabFold:
    CATEGORY = "ComfyBIO/CADD"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("predicted_pdb", "prediction_log")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "fasta_file": ("STRING", {"default": ""}),
            },
            "optional": {
                "num_recycle": ("INT", {"default": 3, "min": 1, "max": 20}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(self, fasta_file: str, output_dir: str = "", num_recycle: int = 3, extra_command: str = ""):
        fa_path = _file(fasta_file, "Protein FASTA")

        executable = shutil.which("colabfold_batch")
        if not executable:
            raise RuntimeError("colabfold_batch executable not found on PATH; install package colabfold")

        out = _output_dir("ColabFold", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        argv = [executable, "--num-recycle", str(num_recycle)]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())
        argv += [str(fa_path), str(out)]

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"ColabFold failed ({proc.returncode}): {proc.stderr}")

        pdb_candidates = list(out.glob("*_relaxed_rank_001*.pdb")) + list(out.glob("*rank_001*.pdb")) + list(out.glob("*.pdb"))
        log_candidates = list(out.glob("*.log"))

        if not pdb_candidates or pdb_candidates[0].stat().st_size == 0:
            raise RuntimeError(f"ColabFold failed to produce PDB structure in {out}")

        return (str(pdb_candidates[0]), str(log_candidates[0]) if log_candidates else "")


NODE_CLASS_MAPPINGS = {"ColabFold": ColabFold}
NODE_DISPLAY_NAME_MAPPINGS = {"ColabFold": "ColabFold: AlphaFold2 Structure Prediction"}


# Backward compatibility aliases
ColabFoldNode = ColabFold

__all__ = ["ColabFold",
    "ColabFoldNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
