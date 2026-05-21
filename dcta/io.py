"""I/O helpers for DCTA."""
from __future__ import annotations
import os
import anndata as ad
from dcta import DCTADataError


def load_h5ad(path: str) -> ad.AnnData:
    """Load an .h5ad file from disk.

    Parameters
    ----------
    path : str
        Path to .h5ad file. Must exist.

    Returns
    -------
    AnnData
        The loaded annotated data matrix.

    Raises
    ------
    DCTADataError
        If the file does not exist or cannot be read as h5ad.
    """
    if not os.path.exists(path):
        raise DCTADataError(f"AnnData file not found: {path}")
    try:
        adata = ad.read_h5ad(path)
    except Exception as e:  # noqa: BLE001
        raise DCTADataError(f"Failed to read h5ad {path}: {e}") from e

    # CELLxGENE-formatted h5ad uses integer var_names with gene symbols in
    # var['feature_name']. Promote symbols to var_names so marker lookup works.
    if "feature_name" in adata.var.columns:
        looks_integer = all(str(s).isdigit() for s in list(adata.var_names[:5]))
        if looks_integer:
            adata.var_names = adata.var["feature_name"].astype(str).values
            adata.var_names_make_unique()

    return adata
