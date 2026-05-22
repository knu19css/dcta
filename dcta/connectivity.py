"""DCTA connectivity stage: KNN on a chosen embedding -> cluster -> PAGA.

Phase 2.8: rebuilds KNN on a publisher-provided embedding (X_pca, X_harmony,
or X_umap) rather than a locally-computed diffusion map. Preprocessing /
diffusion stages are no longer required upstream; the input AnnData only
needs the chosen `obsm[use_rep]` already attached (e.g. via
`dcta.io.attach_paper_embedding`).
"""
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
    use_rep: str = "X_harmony",
    anonymize_cluster_labels: bool = True,
) -> ad.AnnData:
    """Build KNN on `use_rep`, cluster, and PAGA.

    use_rep: name of an `obsm` key (e.g. 'X_harmony', 'X_pca', 'X_umap').
      Must already exist in adata.obsm. No PCA/diffusion is recomputed.
    cluster_from_obs: if given, SKIP Leiden and use that obs column as the
      cluster assignment (copied into 'dcta_leiden').

    Stores connectivities into obsp['dcta_knn'], cluster labels into
    obs['dcta_leiden'], and PAGA outputs into uns['dcta']['paga'].
    """
    a = adata.copy() if copy else adata
    if use_rep not in a.obsm:
        raise KeyError(
            f"use_rep={use_rep!r} not in adata.obsm "
            f"(have: {list(a.obsm.keys())}). Pass an embedding name that "
            f"exists (e.g. X_harmony/X_pca/X_umap from paper)."
        )
    rep_shape = a.obsm[use_rep].shape
    n_pcs = int(rep_shape[1])

    sc.pp.neighbors(a, n_neighbors=knn_k, use_rep=use_rep,
                    n_pcs=n_pcs, random_state=random_state, key_added="dcta_nbr")
    a.obsp["dcta_knn"] = a.obsp["dcta_nbr_connectivities"]

    if cluster_from_obs is not None:
        if cluster_from_obs not in a.obs.columns:
            raise KeyError(
                f"cluster_from_obs={cluster_from_obs!r} not in adata.obs"
            )
        raw = a.obs[cluster_from_obs].astype(object).to_numpy()
        raw_str = np.array([
            "__NA__" if (v is None or (isinstance(v, float) and np.isnan(v))
                         or str(v).lower() in ("nan", "none", ""))
            else str(v)
            for v in raw
        ])
        if anonymize_cluster_labels:
            # Anonymize: preserve cluster membership but strip semantic
            # labels (avoid label leakage when the source column is a
            # paper-curated cell-type annotation).
            uniq_raw = sorted(set(raw_str) - {"__NA__"})
            mapping = {lbl: f"C{i:02d}" for i, lbl in enumerate(uniq_raw)}
            mapping["__NA__"] = "__NA__"
            labels_str = np.array([mapping[s] for s in raw_str])
            cluster_label_map = {v: k for k, v in mapping.items() if k != "__NA__"}
            cluster_source = (f"obs:{cluster_from_obs}(anonymized -> "
                              f"C00..C{len(uniq_raw)-1:02d})")
        else:
            labels_str = raw_str
            cluster_label_map = None
            cluster_source = f"obs:{cluster_from_obs}"
        a.obs["dcta_leiden"] = labels_str
        a.obs["dcta_leiden"] = a.obs["dcta_leiden"].astype("category")
    else:
        sc.tl.leiden(a, neighbors_key="dcta_nbr", resolution=resolution,
                     random_state=random_state, key_added="dcta_leiden")
        cluster_source = f"leiden(resolution={resolution})"

    sc.tl.paga(a, groups="dcta_leiden", neighbors_key="dcta_nbr")
    paga_conn = a.uns["paga"]["connectivities"]

    # Cluster centroids on the embedding (used by downstream viz / pseudotime)
    rep = a.obsm[use_rep]
    labels = a.obs["dcta_leiden"].astype(str).values
    uniq_raw = list(np.unique(labels))
    uniq = (sorted(uniq_raw, key=int) if all(str(u).isdigit() for u in uniq_raw)
            else sorted(uniq_raw))
    centroids = np.vstack([rep[labels == k].mean(axis=0) for k in uniq])
    cluster_sizes = np.array([int((labels == k).sum()) for k in uniq])

    uns = a.uns.setdefault("dcta", {})
    uns["paga"] = {
        "connectivities": paga_conn.toarray() if hasattr(paga_conn, "toarray") else np.asarray(paga_conn),
        "cluster_sizes": cluster_sizes,
        "cluster_centroids": centroids,
        "cluster_labels": uniq,
    }
    uns["cluster_source"] = cluster_source
    uns["use_rep"] = use_rep
    uns["n_rep_dims"] = n_pcs
    if cluster_from_obs is not None and anonymize_cluster_labels:
        # Persist the C## -> original-label mapping so downstream
        # debugging/viz can still see what each anonymized cluster maps to.
        uns["cluster_label_map"] = cluster_label_map
    return a
