import hashlib
import json
import urllib.request
from pathlib import Path


_MANIFEST = json.loads(Path(__file__).with_name("official_data.json").read_text())


def fetch_official_data(name: str, cache_dir: Path) -> Path:
    item = _MANIFEST[name]
    target = cache_dir / name
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        with urllib.request.urlopen(item["url"], timeout=60) as response, target.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    if actual != item["sha256"]:
        target.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {name}: {actual}")
    return target
