import hashlib
import os
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class RuntimeArtifact:
    artifact_id: str
    artifact_type: str
    path: Path
    producer_stage_id: str


def write_text_artifact(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def ensure_dir_artifact(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_project_root() -> Path:
    """Return the repository root directory."""
    if "COMFYBIO_ROOT" in os.environ:
        return Path(os.environ["COMFYBIO_ROOT"]).resolve()
    # Fallback to repository root relative to this file
    return Path(__file__).resolve().parents[4]


def resolve_input_path(path_str: str) -> Path:
    """Resolve an input file or directory path against the project root if relative."""
    p = Path(path_str).expanduser()
    if not p.is_absolute():
        p = get_project_root() / p
    return p.resolve()


COMMON_BIO_EXTENSIONS = {
    ".fa", ".fasta", ".fna", ".fastq", ".fq", ".gz", ".bam", ".sam", ".cram",
    ".vcf", ".bcf", ".gff", ".gtf", ".bed", ".narrowpeak", ".broadpeak",
    ".tsv", ".csv", ".txt", ".json", ".h5ad", ".sf", ".h5", ".loom",
    ".tab", ".counts", ".matrix", ".mtx", ".sh", ".html", ".log",
}


def is_probable_path(p_str: str) -> bool:
    """Determine if a string represents a file or directory path."""
    if not p_str or len(p_str) > 4096:
        return False
    if any(sep in p_str for sep in ("/", "\\")):
        return True
    if p_str.startswith("."):
        return True
    p = Path(p_str)
    if p.suffix.lower() in COMMON_BIO_EXTENSIONS:
        return True
    try:
        if (Path.cwd() / p_str).exists():
            return True
        if resolve_input_path(p_str).exists():
            return True
    except Exception:
        pass
    return False


# Upper bound on the number of files hashed per directory input. Directory scans that
# reach this ceiling cannot prove the inputs are unchanged, so they force re-execution.
MAX_DIR_SCAN_FILES = 20000


def compute_input_fingerprint(*args: Any, **kwargs: Any) -> str:
    """Compute a deterministic input fingerprint for ComfyUI IS_CHANGED caching.

    Inspects all input parameters. For file/directory paths, captures normalized path,
    file size, and modification time (st_mtime_ns). For directories, recursively
    inspects child file sizes and mtimes to detect internal modifications (e.g. quant.sf).
    If any referenced input file is missing, or a directory holds more than
    MAX_DIR_SCAN_FILES files so the scan cannot cover them all, returns a unique token to
    force re-execution.
    """
    hasher = hashlib.sha256()

    def process_item(key: str, val: Any):
        if val is None:
            hasher.update(f"{key}:None;".encode("utf-8"))
            return

        if isinstance(val, (int, float, bool)):
            hasher.update(f"{key}:{type(val).__name__}:{val};".encode("utf-8"))
            return

        if isinstance(val, str):
            val_str = val.strip()
            # Check if this string looks like a path or comma-separated list of paths
            tokens = [p.strip() for p in val_str.split(",") if p.strip()]
            if tokens and all(is_probable_path(t) for t in tokens):
                for p_str in tokens:
                    try:
                        resolved = (Path.cwd() / p_str).resolve() if (Path.cwd() / p_str).exists() else resolve_input_path(p_str)
                        if resolved.is_file():
                            st = resolved.stat()
                            hasher.update(f"{key}:file:{resolved}:{st.st_size}:{st.st_mtime_ns};".encode("utf-8"))
                        elif resolved.is_dir():
                            st = resolved.stat()
                            hasher.update(f"{key}:dir:{resolved}:{st.st_mtime_ns};".encode("utf-8"))
                            # Recursively inspect files inside directory to track content/size/mtime changes
                            file_count = 0
                            for child in sorted(resolved.rglob("*")):
                                if any(part.startswith(".") for part in child.parts):
                                    continue
                                if child.is_file():
                                    c_st = child.stat()
                                    rel = child.relative_to(resolved)
                                    hasher.update(f"f:{rel}:{c_st.st_size}:{c_st.st_mtime_ns};".encode("utf-8"))
                                    file_count += 1
                                    if file_count >= MAX_DIR_SCAN_FILES:
                                        # Scan ceiling reached: the fingerprint no longer covers every
                                        # input file, so force re-execution instead of silently
                                        # reusing a cache entry that may be stale.
                                        hasher.update(
                                            f"{key}:dir_scan_truncated:{resolved}:{os.urandom(8).hex()};".encode("utf-8")
                                        )
                                        break
                        else:
                            # Missing file: mark non-existent to avoid false caching
                            hasher.update(f"{key}:missing:{p_str}:{os.urandom(8).hex()};".encode("utf-8"))
                    except Exception:
                        hasher.update(f"{key}:error:{p_str}:{os.urandom(8).hex()};".encode("utf-8"))
                return

            hasher.update(f"{key}:str:{val_str};".encode("utf-8"))
            return

        if isinstance(val, (list, tuple, set)):
            for i, elem in enumerate(val):
                process_item(f"{key}[{i}]", elem)
            return

        if isinstance(val, dict):
            for k, v in sorted(val.items()):
                process_item(f"{key}.{k}", v)
            return

        hasher.update(f"{key}:{repr(val)};".encode("utf-8"))

    for idx, arg in enumerate(args):
        process_item(f"arg_{idx}", arg)

    for k, v in sorted(kwargs.items()):
        process_item(k, v)

    return hasher.hexdigest()


def get_run_output_dir(base_dir: str | Path, node_name: str, sample_id: str = "", run_id: str = "") -> Path:
    """Construct an isolated output directory path: base / [run_id] / node_name / [sample_id]."""
    out = Path(base_dir).expanduser().resolve()
    if run_id.strip():
        out = out / run_id.strip()
    out = out / node_name
    if sample_id.strip():
        out = out / sample_id.strip()
    out.mkdir(parents=True, exist_ok=True)
    return out


def validate_artifact(
    path: Path | str,
    artifact_type: str,
    allow_empty: bool = False,
) -> tuple[bool, str]:
    """Validate that an output artifact exists and meets type-specific integrity criteria.

    Distinguishes legitimate biological empty results (e.g. 0 significant DEGs) from process failures.
    Returns (is_valid, reason).
    """
    p = Path(path).resolve()
    if not p.exists():
        return False, f"Artifact does not exist: {p}"

    if p.is_dir():
        children = list(p.iterdir())
        if not children and not allow_empty:
            return False, f"Directory artifact is empty: {p}"
        return True, "Valid directory artifact"

    size = p.stat().st_size

    if size == 0:
        if allow_empty:
            return True, f"Empty artifact allowed: {p}"
        if artifact_type == "deg_csv" or artifact_type == "deg_table":
            return False, f"DEG table is 0 bytes (missing header): {p}"
        if artifact_type == "fastq":
            return False, f"FASTQ file is 0 bytes: {p}"
        return False, f"Artifact is 0 bytes: {p}"

    # Type-specific validation
    if artifact_type in ("deg_csv", "deg_table"):
        # Check if CSV/TSV has header and validate column structure
        try:
            with open(p, "r", encoding="utf-8") as f:
                header = f.readline().strip()
            if not header:
                return False, f"DEG table has empty header: {p}"
            # Even if row count == 0 (no DEGs), having a valid header is a valid biological empty result!
            return True, "Valid DEG table artifact"
        except Exception as e:
            return False, f"Failed to parse DEG table: {e}"

    if artifact_type == "fastq":
        # Check if file has data
        return True, "Valid FASTQ artifact"

    return True, "Valid artifact"


