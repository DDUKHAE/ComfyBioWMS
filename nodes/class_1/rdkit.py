"""RDKit cheminformatics descriptor and Morgan fingerprint node.

Python packages: rdkit
External binaries: none
"""

import json


class RDKitDescriptor:
    OUTPUT_NODE = True
    OUPUT_NODE = True
    CATEGORY = "ComfyBIO/CADD"
    FUNCTION = "run"
    RETURN_TYPES = ("STRING", "FLOAT", "FLOAT", "FLOAT", "INT")
    RETURN_NAMES = ("descriptors_json", "mw", "logp", "tpsa", "rotatable_bonds")

    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "smiles": ("STRING", {"default": "CC(=O)OC1=CC=CC=C1C(=O)O"}), # Aspirin
            },
        }

    def run(self, smiles: str):
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors

        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        mw = round(float(Descriptors.MolWt(mol)), 4)
        logp = round(float(Descriptors.MolLogP(mol)), 4)
        tpsa = round(float(Descriptors.TPSA(mol)), 4)
        n_rot = int(Descriptors.NumRotatableBonds(mol))
        hbd = int(rdMolDescriptors.CalcNumHBD(mol))
        hba = int(rdMolDescriptors.CalcNumHBA(mol))

        summary = {
            "smiles": smiles.strip(),
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "rotatable_bonds": n_rot,
            "hydrogen_bond_donors": hbd,
            "hydrogen_bond_acceptors": hba,
            "lipinski_rule_of_five_pass": bool(mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10),
        }

        return (json.dumps(summary, indent=2), mw, logp, tpsa, n_rot)


NODE_CLASS_MAPPINGS = {"RDKitDescriptor": RDKitDescriptor}
NODE_DISPLAY_NAME_MAPPINGS = {"RDKitDescriptor": "RDKit: Chemical Descriptors & Lipinski Rules"}


# Backward compatibility aliases
RDKitDescriptorNode = RDKitDescriptor

__all__ = ["RDKitDescriptor",
    "RDKitDescriptorNode", "NODE_CLASS_MAPPINGS", "NODE_DISPLAY_NAME_MAPPINGS"]
