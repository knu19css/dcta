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
    return adata
