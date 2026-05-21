"""I/O helpers for DCTA."""
from __future__ import annotations
import os
import numpy as np
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


def attach_paper_embedding(
    adata: ad.AnnData,
    source_h5ad: str,
    *,
    match_column: str = "observation_joinid",
    embeddings: tuple[str, ...] = ("X_umap", "X_pca", "X_harmony"),
) -> ad.AnnData:
    """Attach paper-original embeddings from a source h5ad onto adata.

    Matches cells on a shared metadata column (default 'observation_joinid'
    used by CELLxGENE Census). Only the requested embeddings present in the
    source are copied. Cells in adata not found in source receive NaN rows
    in the merged embedding.

    Implements feedback_use_original_embedding: prefer published embedding
    over locally-recomputed UMAP/PCA.

    Parameters
    ----------
    adata : AnnData
        The pipeline's working AnnData. Must contain `match_column` in obs.
    source_h5ad : str
        Path to the publisher's source h5ad (downloaded e.g. via
        cellxgene_census.download_source_h5ad). Read in backed mode.
    match_column : str
        adata.obs column used as the cell identity key. Must exist in both
        adata.obs and source.obs.
    embeddings : tuple of str
        Names of obsm keys to copy from source if present.

    Returns
    -------
    AnnData
        Same `adata` (mutated) with `obsm[k]` populated for each matched key.

    Raises
    ------
    DCTADataError
        If source file missing, match_column missing in either side, or
        zero cells match.
    """
    if not os.path.exists(source_h5ad):
        raise DCTADataError(f"source h5ad not found: {source_h5ad}")
    if match_column not in adata.obs.columns:
        raise DCTADataError(
            f"adata.obs is missing match_column={match_column!r}"
        )
    src = ad.read_h5ad(source_h5ad, backed="r")
    if match_column not in src.obs.columns:
        raise DCTADataError(
            f"source h5ad obs is missing match_column={match_column!r}"
        )

    src_ids = src.obs[match_column].astype(str).to_numpy()
    target_ids = adata.obs[match_column].astype(str).to_numpy()

    # Build src_index lookup: id -> row position
    src_pos = {x: i for i, x in enumerate(src_ids)}
    target_pos = np.array([src_pos.get(t, -1) for t in target_ids], dtype=int)
    n_matched = int((target_pos >= 0).sum())
    if n_matched == 0:
        raise DCTADataError(
            f"No cells matched between adata and source on {match_column!r}"
        )

    for key in embeddings:
        if key not in src.obsm:
            continue
        src_mat = src.obsm[key]
        # backed mode: pull only matched rows
        if hasattr(src_mat, "shape"):
            ncols = src_mat.shape[1]
        else:
            continue
        out = np.full((adata.n_obs, ncols), np.nan, dtype=np.float32)
        valid = target_pos >= 0
        # Read needed rows (may require sorted indices for some backends)
        idx = target_pos[valid]
        order = np.argsort(idx)
        sorted_idx = idx[order]
        sorted_rows = np.asarray(src_mat[sorted_idx, :])
        # Restore original order
        unsorted_rows = np.empty_like(sorted_rows)
        unsorted_rows[order] = sorted_rows
        out[valid] = unsorted_rows
        adata.obsm[key] = out

    adata.uns.setdefault("dcta", {})["paper_embedding_source"] = {
        "path": source_h5ad,
        "match_column": match_column,
        "n_matched": n_matched,
        "n_target": int(adata.n_obs),
        "embeddings_attached": [k for k in embeddings if k in src.obsm],
    }
    return adata
