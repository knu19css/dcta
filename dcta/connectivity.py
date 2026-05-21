"""DCTA connectivity stage: diffusion-based KNN -> Leiden -> PAGA."""
from __future__ import annotations
import numpy as np
import scanpy as sc
import anndata as ad
from dcta import defaults as D


def run(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    knn_k: int = D.KNN_DIFFMAP_K,
    resolution: float = D.LEIDEN_RESOLUTION,
    random_state: int = D.RANDOM_STATE,
) -> ad.AnnData:
    """Build diffusion-based KNN, run Leiden clustering and PAGA.

    Stores connectivities into obsp['dcta_knn'], cluster labels into
    obs['dcta_leiden'], and PAGA outputs into uns['dcta']['paga'].
    """
    a = adata.copy() if copy else adata

    n_dcs = int(a.uns.get("dcta", {}).get("n_dcs", a.obsm["X_diffmap"].shape[1] - 1))
    # Build neighbors on diffusion components
    sc.pp.neighbors(a, n_neighbors=knn_k, use_rep="X_diffmap",
                    n_pcs=n_dcs, random_state=random_state, key_added="dcta")
    # Move the resulting connectivities into obsp['dcta_knn'] for clarity
    a.obsp["dcta_knn"] = a.obsp["dcta_connectivities"]

    sc.tl.leiden(a, neighbors_key="dcta", resolution=resolution,
                 random_state=random_state, key_added="dcta_leiden")

    sc.tl.paga(a, groups="dcta_leiden", neighbors_key="dcta")
    paga_conn = a.uns["paga"]["connectivities"]
    # Cluster centroids in diffusion-map space (for downstream viz / metrics)
    dm = a.obsm["X_diffmap"][:, :n_dcs]
    labels = a.obs["dcta_leiden"].astype(str).values
    uniq = sorted(np.unique(labels), key=int)
    centroids = np.vstack([dm[labels == k].mean(axis=0) for k in uniq])
    cluster_sizes = np.array([int((labels == k).sum()) for k in uniq])

    uns = a.uns.setdefault("dcta", {})
    uns["paga"] = {
        "connectivities": paga_conn.toarray() if hasattr(paga_conn, "toarray") else np.asarray(paga_conn),
        "cluster_sizes": cluster_sizes,
        "cluster_centroids": centroids,
        "cluster_labels": uniq,
    }
    return a
