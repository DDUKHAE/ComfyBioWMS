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
        pd.core.index = _base
        sys.modules["pandas.core.index"] = _base

        # Alias pandas.core.indexes.numeric -> pandas.core.indexes.base
        if hasattr(pd.core, "indexes"):
            pd.core.indexes.numeric = _base
        sys.modules["pandas.core.indexes.numeric"] = _base

        # Ensure legacy index class names resolve to Index
        for idx_name in ("Int64Index", "Float64Index", "UInt64Index", "NumericIndex"):
            if not hasattr(_base, idx_name):
                setattr(_base, idx_name, pd.Index)

        for attr in dir(pd):
            if attr.endswith("Index") and hasattr(pd, attr):
                if not hasattr(_base, attr):
                    setattr(_base, attr, getattr(pd, attr))

        # Ensure is_categorical exists on pandas.api.types and pandas.core.dtypes.common (removed in Pandas 2.0+)
        def _is_categorical(arr_or_dtype) -> bool:
            if arr_or_dtype is None:
                return False
            if isinstance(arr_or_dtype, pd.CategoricalDtype):
                return True
            if hasattr(arr_or_dtype, "dtype") and isinstance(arr_or_dtype.dtype, pd.CategoricalDtype):
                return True
            if isinstance(arr_or_dtype, pd.Categorical):
                return True
            if hasattr(pd, "CategoricalIndex") and isinstance(arr_or_dtype, pd.CategoricalIndex):
                return True
            if hasattr(pd.api.types, "is_categorical_dtype"):
                try:
                    return bool(pd.api.types.is_categorical_dtype(arr_or_dtype))
                except Exception:
                    pass
            return getattr(arr_or_dtype, "name", None) == "category" or str(arr_or_dtype) == "category"

        if not hasattr(pd.api.types, "is_categorical"):
            setattr(pd.api.types, "is_categorical", _is_categorical)
        if "pandas.api.types" in sys.modules:
            setattr(sys.modules["pandas.api.types"], "is_categorical", _is_categorical)

        try:
            import pandas.core.dtypes.common as _common
            if not hasattr(_common, "is_categorical"):
                setattr(_common, "is_categorical", _is_categorical)
            if "pandas.core.dtypes.common" in sys.modules:
                setattr(sys.modules["pandas.core.dtypes.common"], "is_categorical", _is_categorical)
        except Exception:
            pass
    except Exception:
        pass


# Run shim immediately on import
ensure_pandas_compat()
