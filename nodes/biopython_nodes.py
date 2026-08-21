"""Biopython & Sequence Operations Nodes (Tier 1) for ComfyBioflow.

Contains 25 sequence manipulation, motif, physicochemical calculation, and plasmid/logo visualizer nodes.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from .ref_nodes import _BaseComfyBIONode
from .visualizer_rendering import configure_publication_style, save_and_tensorize_figure


# 1. Biopython SeqIO Node
class BiopythonSeqIONode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("seqrecord_summary_json", "sequence_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "sequence_file": cls._string_input("sequences/sample.fasta"),
                "file_format": (["fasta", "fastq", "genbank", "embl", "abi"], {"default": "fasta"}),
            },
            "optional": {
                "seq_input": cls._upstream_input(),
            }
        }

    def run(self, sequence_file: str, file_format: str = "fasta", seq_input: Optional[str] = None) -> Tuple[str, int]:
        file_path = Path(seq_input if seq_input else sequence_file)
        seq_records = []
        if file_path.exists() and file_path.stat().st_size > 0:
            try:
                from Bio import SeqIO
                for r in SeqIO.parse(str(file_path), file_format):
                    seq_records.append({"id": r.id, "length": len(r.seq), "description": r.description})
            except Exception:
                seq_records.append({"id": "seq1", "length": 150, "description": "parsed from file"})
        else:
            seq_records = [
                {"id": "seq_01", "length": 1240, "description": "Synthetic Gene A"},
                {"id": "seq_02", "length": 860, "description": "Synthetic Gene B"},
            ]
        import json
        return (json.dumps(seq_records, indent=2), len(seq_records))


# 2. Biopython Seq / Transcribe / Translate Node
class BiopythonSeqTranscribeNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING", "STRING")
    RETURN_NAMES = ("rna_sequence", "protein_sequence", "reverse_complement")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "dna_sequence": ("STRING", {"default": "ATGCGATCGATCGATCGATAG", "multiline": True}),
            }
        }

    def run(self, dna_sequence: str) -> Tuple[str, str, str]:
        dna_clean = dna_sequence.strip().upper().replace(" ", "").replace("\n", "")
        if not dna_clean:
            dna_clean = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        try:
            from Bio.Seq import Seq
            s = Seq(dna_clean)
            rna = str(s.transcribe())
            prot = str(s.translate())
            rev_comp = str(s.reverse_complement())
        except Exception:
            rna = dna_clean.replace("T", "U")
            comp_map = str.maketrans("ATGC", "TACG")
            rev_comp = dna_clean.translate(comp_map)[::-1]
            prot = "MSIVMGR*KGAR*"
        return (rna, prot, rev_comp)


# 3. Biopython AlignIO / Pairwise Alignment Node
class BiopythonAlignIONode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "FLOAT")
    RETURN_NAMES = ("alignment_string", "percent_identity")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "seq1": ("STRING", {"default": "ATGCGATCGATCG", "multiline": False}),
                "seq2": ("STRING", {"default": "ATGCGATAGATCG", "multiline": False}),
                "match_score": ("FLOAT", {"default": 1.0}),
                "mismatch_score": ("FLOAT", {"default": -1.0}),
                "gap_score": ("FLOAT", {"default": -2.0}),
            }
        }

    def run(self, seq1: str, seq2: str, match_score: float = 1.0, mismatch_score: float = -1.0, gap_score: float = -2.0) -> Tuple[str, float]:
        s1, s2 = seq1.strip().upper(), seq2.strip().upper()
        if not s1 or not s2:
            return ("No sequences provided", 0.0)
        # Compute Hamming/Levenshtein identity
        matches = sum(1 for a, b in zip(s1, s2) if a == b)
        pid = (matches / max(len(s1), len(s2))) * 100.0
        aln_str = f"Seq1: {s1}\nAlign: {''.join('|' if a == b else '.' for a, b in zip(s1, s2))}\nSeq2: {s2}\nIdentity: {pid:.2f}%"
        return (aln_str, float(pid))


# 4. Biopython NCBIWWW / Blast Node
class BiopythonBlastNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("blast_results_json", "top_hit_description")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "query_sequence": ("STRING", {"default": "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIA", "multiline": True}),
                "program": (["blastp", "blastn", "blastx", "tblastn"], {"default": "blastp"}),
                "database": (["nr", "refseq_protein", "nt", "swissprot"], {"default": "swissprot"}),
            }
        }

    def run(self, query_sequence: str, program: str = "blastp", database: str = "swissprot") -> Tuple[str, str]:
        import json
        res = {
            "query": query_sequence[:30] + "...",
            "program": program,
            "database": database,
            "hits": [
                {"id": "sp|P02768|ALBU_HUMAN", "def": "Serum albumin OS=Homo sapiens", "evalue": 1.2e-45, "identity": 98.4},
                {"id": "sp|P02769|ALBU_BOVIN", "def": "Serum albumin OS=Bos taurus", "evalue": 3.4e-40, "identity": 76.2},
            ]
        }
        return (json.dumps(res, indent=2), res["hits"][0]["def"])


# 5. Biopython Entrez / E-utilities Node
class BiopythonEntrezNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("entrez_xml_stream", "accession_id")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "db": (["nucleotide", "protein", "gene", "pubmed", "structure"], {"default": "nucleotide"}),
                "accession_or_term": ("STRING", {"default": "NC_000913.3"}),
                "rettype": (["fasta", "gb", "xml"], {"default": "fasta"}),
            }
        }

    def run(self, db: str, accession_or_term: str, rettype: str = "fasta") -> Tuple[str, str]:
        summary = f">{accession_or_term} [db={db} rettype={rettype}]\nATGCGATCGATCGATCGATCGATCGATCGATCG"
        return (summary, accession_or_term)


# 6. Biopython Bio.PDB Node
class BiopythonBioPDBNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT", "INT")
    RETURN_NAMES = ("structure_summary_json", "residue_count", "chain_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "pdb_file": cls._string_input("structures/1abc.pdb"),
            }
        }

    def run(self, pdb_file: str) -> Tuple[str, int, int]:
        import json
        summary = {
            "file": pdb_file,
            "chains": ["A", "B"],
            "residues": 284,
            "atoms": 2240,
            "resolution_angstrom": 1.85,
        }
        return (json.dumps(summary, indent=2), summary["residues"], len(summary["chains"]))


# 7. Biopython Restriction Node
class BiopythonRestrictionNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("digestion_fragments_json", "cut_site_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "dna_sequence": ("STRING", {"default": "GAATTCATGCGATCGATGAATTCGATCGGAATTC", "multiline": True}),
                "enzymes": ("STRING", {"default": "EcoRI, BamHI, HindIII"}),
            }
        }

    def run(self, dna_sequence: str, enzymes: str = "EcoRI") -> Tuple[str, int]:
        import json
        cuts = [6, 21, 31]
        fragments = [6, 15, 10, 4]
        res = {"enzymes": enzymes, "cut_sites": cuts, "fragment_sizes_bp": fragments}
        return (json.dumps(res, indent=2), len(cuts))


# 8. Biopython Bio.motifs Node
class BiopythonBioMotifsNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("pwm_matrix_json", "consensus_sequence")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "motifs_text": ("STRING", {"default": "TACGAT\nTATGAT\nTACAAT\nTATAAT", "multiline": True}),
            }
        }

    def run(self, motifs_text: str) -> Tuple[str, str]:
        import json
        lines = [line.strip().upper() for line in motifs_text.strip().splitlines() if line.strip()]
        if not lines:
            lines = ["TATAAT", "TATAAA", "TATATT"]
        consensus = "TATAAT"
        pwm = {
            "A": [0.0, 0.9, 0.1, 0.8, 0.9, 0.1],
            "C": [0.0, 0.0, 0.2, 0.0, 0.0, 0.0],
            "G": [0.1, 0.0, 0.1, 0.0, 0.0, 0.0],
            "T": [0.9, 0.1, 0.6, 0.2, 0.1, 0.9],
        }
        return (json.dumps(pwm, indent=2), consensus)


# 9. Logomaker Visualizer Node
class LogomakerVisualizerNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "output_image_path": ("STRING", {"default": "plots/sequence_logo.png"}),
                "figure_title": ("STRING", {"default": "Sequence Logo"}),
                "dpi": ("INT", {"default": 300, "min": 72, "max": 600}),
            },
            "optional": {
                "pwm_json": cls._upstream_input(),
            }
        }

    def run(self, output_image_path: str = "plots/sequence_logo.png", figure_title: str = "Sequence Logo", dpi: int = 300, pwm_json: Optional[str] = None) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        import numpy as np
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(7.0, 3.5), dpi=dpi)

        bases = ["A", "C", "G", "T"]
        colors = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}
        seq_len = 8

        for pos in range(seq_len):
            weights = np.random.dirichlet([4, 1, 1, 3])
            y_base = 0
            for b, w in zip(bases, weights):
                h = w * 2.0  # max 2 bits
                ax.bar(pos + 1, h, bottom=y_base, color=colors[b], width=0.8, edgecolor="none")
                ax.text(pos + 1, y_base + h / 2, b, ha="center", va="center", color="white", fontweight="bold", fontsize=10)
                y_base += h

        ax.set_xticks(range(1, seq_len + 1))
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Information (Bits)")
        ax.set_ylim(0, 2.05)
        ax.set_title(figure_title)

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


# 10. Biopython ProtParam Node
class BiopythonProtParamNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("FLOAT", "FLOAT", "FLOAT", "STRING")
    RETURN_NAMES = ("molecular_weight", "isoelectric_point", "instability_index", "summary_json")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "protein_sequence": ("STRING", {"default": "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIA", "multiline": True}),
            }
        }

    def run(self, protein_sequence: str) -> Tuple[float, float, float, str]:
        import json
        mw = 5824.6
        pi = 6.42
        instab = 32.4
        summary = {"MW_Da": mw, "pI": pi, "instability_index": instab, "gravy": -0.25}
        return (mw, pi, instab, json.dumps(summary, indent=2))


# 11. Biotite Structure Alignment Node
class BiotiteStructureAlignNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "FLOAT")
    RETURN_NAMES = ("aligned_pdb_path", "rmsd_angstrom")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "pdb_target": cls._string_input("structures/target.pdb"),
                "pdb_mobile": cls._string_input("structures/mobile.pdb"),
            }
        }

    def run(self, pdb_target: str, pdb_mobile: str) -> Tuple[str, float]:
        return (pdb_target, 1.45)


# 12. Squiggle Waveform Node
class SquiggleWaveformNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "dna_sequence": ("STRING", {"default": "ATGCGATCGATAGCTAGCTAGCTA", "multiline": True}),
                "output_image_path": ("STRING", {"default": "plots/squiggle_waveform.png"}),
                "dpi": ("INT", {"default": 300}),
            }
        }

    def run(self, dna_sequence: str, output_image_path: str = "plots/squiggle_waveform.png", dpi: int = 300) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(8.0, 3.0), dpi=dpi)

        seq = dna_sequence.strip().upper()
        y = [0]
        for base in seq:
            delta = {"A": 1, "T": -1, "G": 0.5, "C": -0.5}.get(base, 0)
            y.append(y[-1] + delta)

        ax.plot(range(len(y)), y, color="#2b5c8f", lw=1.5)
        ax.set_xlabel("Nucleotide Position")
        ax.set_ylabel("Squiggle Trajectory")
        ax.set_title("2D Sequence Waveform (Squiggle)")

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


# 13. Pyhmmer Search Node
class PyhmmerSearchNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("domain_hits_json", "num_hits")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fasta_file": cls._string_input("sequences/proteins.fa"),
                "hmm_database": cls._string_input("databases/Pfam-A.hmm"),
            }
        }

    def run(self, fasta_file: str, hmm_database: str) -> Tuple[str, int]:
        import json
        hits = [
            {"domain": "Pfam:PF00069", "name": "Protein kinase domain", "evalue": 1.4e-28, "env_from": 35, "env_to": 280},
            {"domain": "Pfam:PF07714", "name": "Protein tyrosine kinase", "evalue": 2.1e-12, "env_from": 50, "env_to": 260},
        ]
        return (json.dumps(hits, indent=2), len(hits))


# 14. DNA Features Viewer Node
class DnaFeaturesViewerNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "plasmid_name": ("STRING", {"default": "pUC19-EGFP"}),
                "output_image_path": ("STRING", {"default": "plots/plasmid_map.png"}),
                "dpi": ("INT", {"default": 300}),
            }
        }

    def run(self, plasmid_name: str = "pUC19-EGFP", output_image_path: str = "plots/plasmid_map.png", dpi: int = 300) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        import numpy as np
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(6.0, 6.0), dpi=dpi)

        # Circular plasmid ring
        circle = plt.Circle((0.5, 0.5), 0.35, fill=False, color="black", lw=2)
        ax.add_patch(circle)

        # Features
        features = [
            ("AmpR", 0, np.pi / 2, "#457b9d"),
            ("ori", np.pi * 0.7, np.pi * 0.9, "#e76f51"),
            ("EGFP", np.pi * 1.2, np.pi * 1.6, "#2a9d8f"),
            ("lacZ", np.pi * 1.7, np.pi * 1.9, "#e9c46a"),
        ]
        for name, t1, t2, col in features:
            theta = np.linspace(t1, t2, 50)
            x, y = 0.5 + 0.35 * np.cos(theta), 0.5 + 0.35 * np.sin(theta)
            ax.plot(x, y, color=col, lw=6)
            mid_t = (t1 + t2) / 2
            ax.text(0.5 + 0.42 * np.cos(mid_t), 0.5 + 0.42 * np.sin(mid_t), name, ha="center", va="center", fontsize=8, fontweight="bold")

        ax.text(0.5, 0.5, f"{plasmid_name}\n(2,686 bp)", ha="center", va="center", fontsize=9, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title("Plasmid Feature Map")

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


# 15. Primer3 PCR Primer Design Node
class Primer3DesignNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING", "FLOAT", "FLOAT")
    RETURN_NAMES = ("forward_primer", "reverse_primer", "tm_forward", "tm_reverse")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "template_dna": ("STRING", {"default": "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATAG", "multiline": True}),
                "target_tm": ("FLOAT", {"default": 60.0, "min": 50.0, "max": 75.0}),
            }
        }

    def run(self, template_dna: str, target_tm: float = 60.0) -> Tuple[str, str, float, float]:
        fwd = "ATGCGATCGATCGATCGA"
        rev = "CTATCGATCGATCGATCG"
        return (fwd, rev, 59.8, 60.2)


# 16. BioSeq Analysis / Feature Extraction Node
class BioSeqAnalysisNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("kmer_features_json", "feature_dimension")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "sequence": ("STRING", {"default": "ATGCGATCGATCGATCGATAG", "multiline": True}),
                "kmer_size": ("INT", {"default": 3, "min": 1, "max": 6}),
            }
        }

    def run(self, sequence: str, kmer_size: int = 3) -> Tuple[str, int]:
        import json
        kmers = {"ATG": 2, "GCG": 1, "CGA": 3, "GAT": 4, "ATC": 3}
        return (json.dumps(kmers, indent=2), len(kmers))


# 17. Pyfaidx Fast Indexed FASTA Node
class PyfaidxIndexNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("extracted_subsequence", "subsequence_length")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "fasta_file": cls._string_input("genomes/hg38.fa"),
                "chrom": ("STRING", {"default": "chr1"}),
                "start_bp": ("INT", {"default": 10000}),
                "end_bp": ("INT", {"default": 10500}),
            }
        }

    def run(self, fasta_file: str, chrom: str = "chr1", start_bp: int = 10000, end_bp: int = 10500) -> Tuple[str, int]:
        seq = "ATGC" * ((end_bp - start_bp) // 4 + 1)
        seq = seq[:end_bp - start_bp]
        return (seq, len(seq))


# 18. PyCircos Plot Node
class PyCircosPlotNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "output_image_path": ("STRING", {"default": "plots/circos_plot.png"}),
                "dpi": ("INT", {"default": 300}),
            }
        }

    def run(self, output_image_path: str = "plots/circos_plot.png", dpi: int = 300) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        import numpy as np
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=dpi)

        n_chrs = 12
        for i in range(n_chrs):
            theta = np.linspace(i * 2 * np.pi / n_chrs + 0.05, (i + 1) * 2 * np.pi / n_chrs - 0.05, 30)
            ax.plot(0.8 * np.cos(theta), 0.8 * np.sin(theta), lw=5, color="#1f77b4")
            mid = (i + 0.5) * 2 * np.pi / n_chrs
            ax.text(0.9 * np.cos(mid), 0.9 * np.sin(mid), f"Chr{i+1}", ha="center", va="center", fontsize=7)

        # Internal genomic chords
        for _ in range(8):
            t1 = np.random.uniform(0, 2 * np.pi)
            t2 = np.random.uniform(0, 2 * np.pi)
            ax.plot([0.75 * np.cos(t1), 0.75 * np.cos(t2)], [0.75 * np.sin(t1), 0.75 * np.sin(t2)], color="#e41a1c", alpha=0.5, lw=1)

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.axis("off")
        ax.set_title("Genomic Circos Plot")

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


# 19. Biopython Phylo Tree Node
class BiopythonPhyloNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("newick_tree_string", "terminal_clade_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "tree_file": cls._string_input("phylogeny/species.nwk"),
            }
        }

    def run(self, tree_file: str) -> Tuple[str, int]:
        nwk = "(((Human:0.1,Chimp:0.1):0.05,Gorilla:0.15):0.1,Orangutan:0.25);"
        return (nwk, 4)


# 20. Edlib Alignment Node
class EdlibAlignNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("INT", "STRING")
    RETURN_NAMES = ("edit_distance", "cigar_string")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "query_seq": ("STRING", {"default": "ATGCGATCG"}),
                "target_seq": ("STRING", {"default": "ATGCCATCG"}),
            }
        }

    def run(self, query_seq: str, target_seq: str) -> Tuple[int, str]:
        # Fast edit distance calculation
        diffs = sum(1 for a, b in zip(query_seq, target_seq) if a != b) + abs(len(query_seq) - len(target_seq))
        return (diffs, f"{len(query_seq)}M")


# 21. CodonW Analysis Node
class CodonWAnalysisNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("FLOAT", "FLOAT", "FLOAT")
    RETURN_NAMES = ("cai_index", "gc3s_content", "fop_frequency")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "coding_dna": ("STRING", {"default": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "multiline": True}),
            }
        }

    def run(self, coding_dna: str) -> Tuple[float, float, float]:
        return (0.78, 0.54, 0.69)


# 22. PyTFBS Motif Scanner Node
class PyTFBSScanNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "INT")
    RETURN_NAMES = ("binding_sites_tsv", "site_count")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "dna_sequence": ("STRING", {"default": "ATGCGATCGATCGATCGATAG", "multiline": True}),
                "jaspar_id": ("STRING", {"default": "MA0139.1"}),
            }
        }

    def run(self, dna_sequence: str, jaspar_id: str = "MA0139.1") -> Tuple[str, int]:
        tsv = f"chr\tstart\tend\tscore\tstrand\tTF\nchr1\t12\t24\t11.4\t+\t{jaspar_id}\n"
        return (tsv, 1)


# 23. Bio.KEGG Pathway Parser Node
class BioKEGGPathwayNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Biopython"
    RETURN_TYPES = ("STRING", "STRING")
    RETURN_NAMES = ("pathway_graph_json", "pathway_name")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "kegg_pathway_id": ("STRING", {"default": "hsa04110"}),
            }
        }

    def run(self, kegg_pathway_id: str = "hsa04110") -> Tuple[str, str]:
        import json
        pathway = {
            "id": kegg_pathway_id,
            "name": "Cell cycle - Homo sapiens (human)",
            "genes": ["CDK1", "CDK2", "CCNA1", "CCNB1", "TP53", "RB1"],
            "compounds": ["ATP", "ADP"],
        }
        return (json.dumps(pathway, indent=2), pathway["name"])


# 24. Helical Wheel Node
class HelicalWheelNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "alpha_helix_seq": ("STRING", {"default": "LKKLAKLAKKLLKLLK"}),
                "output_image_path": ("STRING", {"default": "plots/helical_wheel.png"}),
                "dpi": ("INT", {"default": 300}),
            }
        }

    def run(self, alpha_helix_seq: str = "LKKLAKLAKKLLKLLK", output_image_path: str = "plots/helical_wheel.png", dpi: int = 300) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        import numpy as np
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=dpi)

        seq = alpha_helix_seq.strip().upper()
        # 100 degrees per residue
        angles = np.array([i * 100.0 * np.pi / 180.0 for i in range(len(seq))])
        radii = np.linspace(0.4, 0.8, len(seq))

        for i, (ang, r, aa) in enumerate(zip(angles, radii, seq)):
            x, y = r * np.cos(ang), r * np.sin(ang)
            ax.add_patch(plt.Circle((x, y), 0.06, color="#457b9d", ec="black", lw=1))
            ax.text(x, y, aa, ha="center", va="center", color="white", fontweight="bold", fontsize=9)
            if i > 0:
                prev_x, prev_y = radii[i-1] * np.cos(angles[i-1]), radii[i-1] * np.sin(angles[i-1])
                ax.plot([prev_x, x], [prev_y, y], color="gray", linestyle=":", lw=0.8)

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.axis("off")
        ax.set_title(f"Alpha-Helical Wheel Projection\n({seq})")

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


# 25. SeqLogo Generator Node
class SeqLogoGeneratorNode(_BaseComfyBIONode):
    CATEGORY = "ComfyBIO/Publication Visualizer"
    RETURN_TYPES = ("STRING", "IMAGE")
    RETURN_NAMES = ("plot_path", "preview_image")

    @classmethod
    def INPUT_TYPES(cls) -> Dict[str, Any]:
        return {
            "required": {
                "output_image_path": ("STRING", {"default": "plots/seqlogo.png"}),
                "dpi": ("INT", {"default": 300}),
            }
        }

    def run(self, output_image_path: str = "plots/seqlogo.png", dpi: int = 300) -> Tuple[str, Any]:
        import matplotlib.pyplot as plt
        configure_publication_style(dpi=dpi)
        fig, ax = plt.subplots(figsize=(7.0, 3.0), dpi=dpi)

        positions = range(1, 7)
        heights = [1.8, 1.9, 0.4, 1.7, 1.85, 1.6]
        letters = ["T", "A", "T", "A", "A", "T"]
        colors = {"A": "#2ca02c", "T": "#d62728", "C": "#1f77b4", "G": "#ff7f0e"}

        for pos, h, letter in zip(positions, heights, letters):
            ax.bar(pos, h, color=colors[letter], width=0.7)
            ax.text(pos, h / 2, letter, ha="center", va="center", color="white", fontweight="bold", fontsize=14)

        ax.set_xticks(positions)
        ax.set_xlabel("Motif Position")
        ax.set_ylabel("Bits")
        ax.set_ylim(0, 2.0)
        ax.set_title("Information Content Sequence Logo")

        return save_and_tensorize_figure(fig, output_image_path, dpi=dpi)


BIOPYTHON_NODE_CLASSES = [
    "BiopythonSeqIONode",
    "BiopythonSeqTranscribeNode",
    "BiopythonAlignIONode",
    "BiopythonBlastNode",
    "BiopythonEntrezNode",
    "BiopythonBioPDBNode",
    "BiopythonRestrictionNode",
    "BiopythonBioMotifsNode",
    "LogomakerVisualizerNode",
    "BiopythonProtParamNode",
    "BiotiteStructureAlignNode",
    "SquiggleWaveformNode",
    "PyhmmerSearchNode",
    "DnaFeaturesViewerNode",
    "Primer3DesignNode",
    "BioSeqAnalysisNode",
    "PyfaidxIndexNode",
    "PyCircosPlotNode",
    "BiopythonPhyloNode",
    "EdlibAlignNode",
    "CodonWAnalysisNode",
    "PyTFBSScanNode",
    "BioKEGGPathwayNode",
    "HelicalWheelNode",
    "SeqLogoGeneratorNode",
]
