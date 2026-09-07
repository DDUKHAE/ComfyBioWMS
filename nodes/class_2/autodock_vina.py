"""AutoDock Vina molecular docking node.

Python packages: none
External binaries: vina
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



class AutoDockVina:
    CATEGORY = "ComfyBIO/CADD"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("docked_poses_pdbqt", "docking_log")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "receptor_pdbqt": ("STRING", {"default": ""}),
                "ligand_pdbqt": ("STRING", {"default": ""}),
            },
            "optional": {
                "center_x": ("FLOAT", {"default": 0.0}),
                "center_y": ("FLOAT", {"default": 0.0}),
                "center_z": ("FLOAT", {"default": 0.0}),
                "size_x": ("FLOAT", {"default": 20.0}),
                "size_y": ("FLOAT", {"default": 20.0}),
                "size_z": ("FLOAT", {"default": 20.0}),
                "exhaustiveness": ("INT", {"default": 8, "min": 1, "max": 64}),
                "cpu": ("INT", {"default": 4, "min": 1, "max": 64}),
                "extra_command": ("STRING", {"default": "", "multiline": True}),
            },
        }

    def run(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        output_dir: str = "",
        center_x: float = 0.0,
        center_y: float = 0.0,
        center_z: float = 0.0,
        size_x: float = 20.0,
        size_y: float = 20.0,
        size_z: float = 20.0,
        exhaustiveness: int = 8,
        cpu: int = 4,
        extra_command: str = "",
    ):
        rec_path = _file(receptor_pdbqt, "Receptor PDBQT")
        lig_path = _file(ligand_pdbqt, "Ligand PDBQT")

        executable = shutil.which("vina")
        if not executable:
            raise RuntimeError("vina executable not found on PATH; install bioconda package autodock-vina")

        out = _output_dir("AutoDockVina", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        out_pdbqt = out / f"{lig_path.stem}_docked.pdbqt"
        out_log = out / f"{lig_path.stem}_vina.log"

        argv = [
            executable,
            "--receptor",
            str(rec_path),
            "--ligand",
            str(lig_path),
            "--out",
            str(out_pdbqt),
            "--log",
            str(out_log),
            "--center_x",
            str(center_x),
            "--center_y",
            str(center_y),
            "--center_z",
            str(center_z),
            "--size_x",
            str(size_x),
            "--size_y",
            str(size_y),
            "--size_z",
            str(size_z),
            "--exhaustiveness",
            str(exhaustiveness),
            "--cpu",
            str(cpu),
        ]
        if extra_command.strip():
            argv += shlex.split(extra_command.strip())

        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"AutoDock Vina docking failed ({proc.returncode}): {proc.stderr}")

        if not out_pdbqt.is_file() or out_pdbqt.stat().st_size == 0:
            raise RuntimeError(f"AutoDock Vina failed to generate output poses at {out_pdbqt}")

        return (str(out_pdbqt), str(out_log) if out_log.is_file() else "")


NODE_CLASS_MAPPINGS = {"AutoDockVina": AutoDockVina}
NODE_DISPLAY_NAME_MAPPINGS = {"AutoDockVina": "AutoDock Vina: Molecular Docking"}


# Backward compatibility aliases
AutoDockVinaNode = AutoDockVina

__all__ = ["AutoDockVina",
    "AutoDockVinaNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
