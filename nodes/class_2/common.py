"""Common utilities for bioinformatics nodes: sample discovery and batch processing."""

import re
from pathlib import Path
from typing import List, Optional, Tuple


def _clean_stem(filename: str) -> str:
    stem = filename
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".sam", ".cram"):
        if stem.lower().endswith(ext):
            stem = stem[:-len(ext)]
            break
    for sfx in ("_val_1", "_val_2", "_R1", "_R2", "_1", "_2", ".1", ".2", "-R1", "-R2"):
        if stem.endswith(sfx):
            stem = stem[:-len(sfx)]
            break
    return stem


def discover_samples(fwd_input: str, rev_input: str = "") -> List[Tuple[str, Path, Optional[Path]]]:
    """Discover sample file pairs or single files from a file path, directory path, or comma-separated paths.

    Returns:
        List of (sample_id, fwd_path, rev_path_or_None)
    """
    fwd_raw = str(fwd_input).strip()
    rev_raw = str(rev_input).strip() if rev_input else ""
    if not fwd_raw:
        raise ValueError("Input reads / file path cannot be empty")

    # 1. Comma-separated list of files or directories
    if "," in fwd_raw:
        fwds = [Path(x.strip()).expanduser().resolve() for x in fwd_raw.split(",") if x.strip()]
        revs = [Path(x.strip()).expanduser().resolve() for x in rev_raw.split(",") if x.strip()] if rev_raw else []
        if rev_raw and len(revs) != len(fwds):
            raise ValueError(f"Mismatched forward ({len(fwds)}) and reverse ({len(revs)}) file counts.")
        samples = []
        for idx, f in enumerate(fwds):
            if not f.is_file() and not f.is_dir():
                raise FileNotFoundError(f"Input path not found: {f}")
            if f.is_dir():
                sub = discover_samples(str(f), str(revs[idx]) if idx < len(revs) else "")
                samples.extend(sub)
                continue
            r = revs[idx] if idx < len(revs) else None
            if r is not None and not r.is_file():
                raise FileNotFoundError(f"Reverse mate file not found: {r}")
            stem = _clean_stem(f.name)
            samples.append((stem, f, r))
        return samples

    fwd_path = Path(fwd_raw).expanduser().resolve()
    rev_path = Path(rev_raw).expanduser().resolve() if rev_raw else None

    # 2. Single file input
    if fwd_path.is_file():
        if rev_raw:
            if not rev_path.is_file():
                raise FileNotFoundError(f"Reverse mate file not found: {rev_path}")
        stem = _clean_stem(fwd_path.name)
        return [(stem, fwd_path, rev_path)]

    # 3. Directory input
    if not fwd_path.exists():
        raise FileNotFoundError(f"Forward input file or directory not found: {fwd_path}")
    if not fwd_path.is_dir():
        raise FileNotFoundError(f"Input path is neither a file nor a directory: {fwd_path}")
    if rev_raw and not rev_path.is_dir():
        raise FileNotFoundError(f"Reverse input directory not found: {rev_path}")

    exts = (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".sam", ".cram")
    all_files = sorted([p for p in fwd_path.rglob("*") if p.is_file() and any(p.name.lower().endswith(e) for e in exts)])

    if not all_files:
        raise FileNotFoundError(f"No sequencing files (FASTQ/BAM/SAM) found in directory: {fwd_path}")

    r1_pattern = re.compile(r"^(.*?)(?:[._-](?:R1|1|val_1))([._-][^._-]+)*(\.fastq.*|\.fq.*)$", re.IGNORECASE)
    r2_pattern = re.compile(r"^(.*?)(?:[._-](?:R2|2|val_2))([._-][^._-]+)*(\.fastq.*|\.fq.*)$", re.IGNORECASE)

    r1_map = {}
    r2_map = {}
    unpaired = []

    for f in all_files:
        name = f.name
        m1 = r1_pattern.match(name)
        m2 = r2_pattern.match(name)
        if m1 and not m2:
            key = m1.group(1)
            r1_map[key] = f
        elif m2 and not m1:
            key = m2.group(1)
            r2_map[key] = f
        else:
            unpaired.append(f)

    samples = []
    for key, fwd_f in r1_map.items():
        if key in r2_map:
            samples.append((key, fwd_f, r2_map[key]))
        else:
            samples.append((key, fwd_f, None))

    unmatched_r2 = [f for k, f in r2_map.items() if k not in r1_map]
    for f in unpaired + unmatched_r2:
        stem = _clean_stem(f.name)
        samples.append((stem, f, None))

    samples.sort(key=lambda s: s[0])
    return samples
