"""DCTA diffusion stage: PCA, KNN, diffusion map with eigval-gap auto-n_dcs."""
from __future__ import annotations
import numpy as np
import scanpy as sc
import anndata as ad
from dcta import defaults as D


def _select_n_dcs_by_gap(eigvals: np.ndarray, min_dcs: int = 2,
                         max_dcs: int | None = None) -> int:
    """Pick n_dcs at the largest relative drop (eigenvalue gap).

    Skips index 0 (trivial eigenvalue ~1.0 of the diffusion operator).
    """
    e = np.asarray(eigvals, dtype=float)
    if max_dcs is None:
        max_dcs = len(e) - 1
    max_dcs = min(max_dcs, len(e) - 1)
    if max_dcs <= min_dcs:
        return max_dcs
    # Consider gaps from min_dcs..max_dcs (between e[i] and e[i+1])
    gaps = e[min_dcs:max_dcs] - e[min_dcs + 1:max_dcs + 1]
    rel_gaps = gaps / np.maximum(np.abs(e[min_dcs:max_dcs]), 1e-9)
    best = int(np.argmax(rel_gaps)) + min_dcs + 1
    return max(min_dcs, min(best, max_dcs))


def run(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    n_comps_pca: int = D.N_COMPS_PCA,
    n_comps: int = D.N_COMPS_DIFFMAP,
    knn_k: int = D.KNN_PCA_K,
    random_state: int = D.RANDOM_STATE,
) -> ad.AnnData:
    """Run PCA -> KNN -> diffusion map, then pick n_dcs from eigenvalue gap.

    Mutates adata (or copy) and returns it.
    """
    a = adata.copy() if copy else adata

    sc.pp.pca(a, n_comps=min(n_comps_pca, a.n_vars - 1, a.n_obs - 1),
              random_state=random_state)
    sc.pp.neighbors(a, n_neighbors=knn_k, n_pcs=a.obsm["X_pca"].shape[1],
                    random_state=random_state)
    sc.tl.diffmap(a, n_comps=min(n_comps, a.n_obs - 2))

    eigvals = np.asarray(a.uns["diffmap_evals"])
    n_dcs = _select_n_dcs_by_gap(eigvals, min_dcs=2, max_dcs=len(eigvals) - 1)

    uns = a.uns.setdefault("dcta", {})
    uns["diffmap_eigvals"] = eigvals.tolist()
    uns["n_dcs"] = int(n_dcs)
    return a
