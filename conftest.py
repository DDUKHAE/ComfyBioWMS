import sys
from pathlib import Path

root = Path(__file__).resolve().parent
parent = root.parent

# Remove parent from sys.path to avoid collisions with sibling/parent directories
sys.path = [p for p in sys.path if Path(p).resolve() != parent.resolve()]

import nodes

if str(root) not in sys.path:
    sys.path.insert(0, str(root))

engine_src = str(root / "engine" / "src")
if engine_src not in sys.path:
    sys.path.insert(1, engine_src)

# Purge any incorrectly resolved 'nodes' modules from sys.modules
for mod_name in list(sys.modules.keys()):
    if mod_name == "nodes" or mod_name.startswith("nodes."):
        mod = sys.modules[mod_name]
        mod_file = getattr(mod, "__file__", "")
        if mod_file and str(parent) in mod_file and str(root) not in mod_file:
            del sys.modules[mod_name]
