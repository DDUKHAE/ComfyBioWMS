"""Compatibility shims for upstream library version discrepancies.

Handles legacy AnnData / Pickle deserialization across Pandas 1.x, 2.x, and 3.x.
"""

import sys


def ensure_pandas_compat() -> None:
    """Ensure backwards compatibility for legacy H5AD/pickle files serialized with Pandas < 2.0.

    Pandas 2.0+ removed 'pandas.core.index' and 'pandas.core.indexes.numeric', which causes
    'ModuleNotFoundError: No module named pandas.core.index' when unpickling legacy AnnData metadata
    (such as cell/gene index objects or uns annotations).
    This shim aliases the removed module paths to 'pandas.core.indexes.base' and maps legacy
    index classes (Int64Index, Float64Index, UInt64Index, NumericIndex) to 'pandas.Index'.
    """
    try:
        import pandas as pd
        import pandas.core.indexes.base as _base

        # Alias pandas.core.index -> pandas.core.indexes.base
        if not hasattr(pd.core, "index"):
            pd.core.index = _base
        if "pandas.core.index" not in sys.modules:
            sys.modules["pandas.core.index"] = _base

        # Alias pandas.core.indexes.numeric -> pandas.core.indexes.base
        if not hasattr(pd.core.indexes, "numeric"):
            pd.core.indexes.numeric = _base
        if "pandas.core.indexes.numeric" not in sys.modules:
            sys.modules["pandas.core.indexes.numeric"] = _base

        # Ensure legacy index class names resolve to Index
        for idx_name in ("Int64Index", "Float64Index", "UInt64Index", "NumericIndex"):
            if not hasattr(_base, idx_name):
                setattr(_base, idx_name, pd.Index)

        for attr in ("MultiIndex", "RangeIndex", "CategoricalIndex"):
            if not hasattr(_base, attr) and hasattr(pd, attr):
                setattr(_base, attr, getattr(pd, attr))
    except Exception:
        pass


# Run shim immediately on import
ensure_pandas_compat()
