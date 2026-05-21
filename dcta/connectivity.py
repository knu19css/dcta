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
    cluster_from_obs: str | None = None,
) -> ad.AnnData:
    """Build diffusion-based KNN, cluster, and PAGA.

    cluster_from_obs: if given, SKIP Leiden and use that obs column as
      the cluster assignment (copied into 'dcta_leiden'). Useful when the
      input dataset already carries publisher-provided cluster labels
      (e.g. paper_AuthorCellType_Broad attached via io.attach_paper_embedding).
      Cells with NaN/None in that column are assigned cluster '__NA__'.

    Stores connectivities into obsp['dcta_knn'], cluster labels into
    obs['dcta_leiden'], and PAGA outputs into uns['dcta']['paga'].
    """
    a = adata.copy() if copy else adata

    n_dcs = int(a.uns.get("dcta", {}).get("n_dcs", a.obsm["X_diffmap"].shape[1] - 1))
    sc.pp.neighbors(a, n_neighbors=knn_k, use_rep="X_diffmap",
                    n_pcs=n_dcs, random_state=random_state, key_added="dcta_nbr")
    a.obsp["dcta_knn"] = a.obsp["dcta_nbr_connectivities"]

    if cluster_from_obs is not None:
        if cluster_from_obs not in a.obs.columns:
            raise KeyError(
                f"cluster_from_obs={cluster_from_obs!r} not in adata.obs"
            )
        raw = a.obs[cluster_from_obs].astype(object).to_numpy()
        labels_str = np.array([
            "__NA__" if (v is None or (isinstance(v, float) and np.isnan(v))
                         or str(v).lower() in ("nan", "none", ""))
            else str(v)
            for v in raw
        ])
        a.obs["dcta_leiden"] = labels_str
        a.obs["dcta_leiden"] = a.obs["dcta_leiden"].astype("category")
        cluster_source = f"obs:{cluster_from_obs}"
    else:
        sc.tl.leiden(a, neighbors_key="dcta_nbr", resolution=resolution,
                     random_state=random_state, key_added="dcta_leiden")
        cluster_source = f"leiden(resolution={resolution})"

    sc.tl.paga(a, groups="dcta_leiden", neighbors_key="dcta_nbr")
    paga_conn = a.uns["paga"]["connectivities"]
    # Cluster centroids in diffusion-map space (for downstream viz / metrics)
    dm = a.obsm["X_diffmap"][:, :n_dcs]
    labels = a.obs["dcta_leiden"].astype(str).values
    # When labels are strings (paper categories), sort lexicographically;
    # when leiden ids (all digits), sort numerically.
    uniq_raw = list(np.unique(labels))
    if all(str(u).isdigit() for u in uniq_raw):
        uniq = sorted(uniq_raw, key=int)
    else:
        uniq = sorted(uniq_raw)
    centroids = np.vstack([dm[labels == k].mean(axis=0) for k in uniq])
    cluster_sizes = np.array([int((labels == k).sum()) for k in uniq])

    uns = a.uns.setdefault("dcta", {})
    uns["paga"] = {
        "connectivities": paga_conn.toarray() if hasattr(paga_conn, "toarray") else np.asarray(paga_conn),
        "cluster_sizes": cluster_sizes,
        "cluster_centroids": centroids,
        "cluster_labels": uniq,
    }
    uns["cluster_source"] = cluster_source
    return a
