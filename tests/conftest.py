import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
parent = root.parent

# Remove any parent paths that might collide with 'nodes'
sys.path = [p for p in sys.path if Path(p).resolve() != parent.resolve()]

if str(root) not in sys.path:
    sys.path.insert(0, str(root))

engine_src = str(root / "engine" / "src")
if engine_src not in sys.path:
    sys.path.insert(1, engine_src)
