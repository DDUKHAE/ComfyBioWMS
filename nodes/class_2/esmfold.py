"""ESMFold ultra-fast single-sequence protein folding node.

Python packages: esm, torch
External binaries: esm-fold (or python script)
"""

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



class ESMFold:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/CADD"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "FLOAT")
    RETURN_NAMES = ("predicted_pdb", "mean_plddt")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "fasta_file": ("STRING", {"default": ""}),
            },
        }

    def run(self, fasta_file: str, output_dir: str = ""):
        import torch
        import esm

        fa_path = _file(fasta_file, "Protein FASTA")
        out = _output_dir("ESMFold", output_dir)
        out.mkdir(parents=True, exist_ok=True)

        lines = [l.strip() for l in open(fa_path).readlines() if l.strip()]
        seq = "".join([l for l in lines if not l.startswith(">")]).upper()

        model = esm.pretrained.esmfold_v1()
        model = model.eval().cuda() if torch.cuda.is_available() else model.eval()

        with torch.no_grad():
            output = model.infer_pdb(seq)

        out_pdb = out / f"{fa_path.stem}_esmfold.pdb"
        with open(out_pdb, "w") as f:
            f.write(output)

        # Parse B-factor column for pLDDT
        plddts = []
        for line in output.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    plddts.append(float(line[60:66].strip()))
                except ValueError:
                    pass

        mean_plddt = round(sum(plddts) / len(plddts), 2) if plddts else 0.0
        return (str(out_pdb), mean_plddt)


NODE_CLASS_MAPPINGS = {"ESMFold": ESMFold}
NODE_DISPLAY_NAME_MAPPINGS = {"ESMFold": "ESMFold: Fast Protein Structure Prediction"}


# Backward compatibility aliases
ESMFoldNode = ESMFold

__all__ = ["ESMFold",
    "ESMFoldNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
