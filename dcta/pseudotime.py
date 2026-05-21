"""DCTA pseudotime stage: multi-root-aware DPT.

Default strategy: HSC marker medoid per connected component, top-K
robustness re-sampling, virtual_hsc_ref hook deferred to Phase 2.
Also stores per-cluster HSC marker score (used by dag.build for marker-direction).
"""
from __future__ import annotations
import numpy as np
import scanpy as sc
import anndata as ad
import yaml
from scipy.sparse.csgraph import connected_components
from dcta import DCTARootError
from dcta import defaults as D


def _load_hsc_markers(path: str) -> list[str]:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return list(cfg.get("hsc_positive", []))


def _score_hsc(adata: ad.AnnData, markers: list[str]) -> np.ndarray:
    """Per-cell mean expression of HSC markers (from log-normalized X)."""
    present = [g for g in markers if g in adata.var_names]
    if not present:
        raise DCTARootError(
            f"None of the HSC marker genes are present in adata.var_names: {markers}"
        )
    idx = [adata.var_names.get_loc(g) for g in present]
    X = adata.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X).mean(axis=1)


def _medoid_cell(coords: np.ndarray) -> int:
    """Index of the cell minimizing summed distance to others (in coords subset)."""
    if coords.shape[0] == 1:
        return 0
    diff = coords[:, None, :] - coords[None, :, :]
    d = np.sqrt((diff ** 2).sum(axis=-1))
    return int(np.argmin(d.sum(axis=1)))


def _connected_components_from_paga(adata: ad.AnnData) -> tuple[np.ndarray, int]:
    """Return (cell_component_id, n_components) using PAGA cluster connectivity."""
    paga = adata.uns["dcta"]["paga"]
    conn = np.asarray(paga["connectivities"])
    labels = list(paga["cluster_labels"])
    adj = (conn > 0).astype(int)
    n_comp, cluster_comp = connected_components(adj, directed=False)
    label_to_comp = dict(zip(labels, cluster_comp))
    cell_comp = adata.obs["dcta_leiden"].astype(str).map(label_to_comp).to_numpy()
    return cell_comp.astype(int), int(n_comp)


def _dpt_from_iroot(adata: ad.AnnData, iroot: int) -> np.ndarray:
    """Run scanpy DPT with given iroot, return per-cell pseudotime array."""
    a = adata
    a.uns["iroot"] = int(iroot)
    sc.tl.dpt(a, neighbors_key="dcta_nbr")
    return a.obs["dpt_pseudotime"].to_numpy().copy()


def run(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    root_strategy: str = "hsc_marker_medoid",
    hsc_marker_yaml: str | None = None,
    robustness_k: int = D.ROBUSTNESS_K,
    ref_adata: ad.AnnData | None = None,
) -> ad.AnnData:
    """Multi-root-aware DPT.

    Strategies:
      - "hsc_marker_medoid" (default): medoid of HSC cluster per connected
        component as iroot. Run DPT once per component. Then re-run K times
        using top-K HSC candidates and store cell-wise std.
      - "virtual_hsc_ref": NotImplementedError (Phase 2).

    Also stores per-cluster HSC marker score and maturation score
    (= negative HSC score for Phase 1; mature-marker YAML is Phase 2).
    """
    if root_strategy == "virtual_hsc_ref":
        raise NotImplementedError("virtual_hsc_ref root strategy is Phase 2")
    if root_strategy != "hsc_marker_medoid":
        raise ValueError(f"Unknown root_strategy: {root_strategy}")
    if hsc_marker_yaml is None:
        raise DCTARootError("hsc_marker_yaml is required for hsc_marker_medoid")

    a = adata.copy() if copy else adata
    markers = _load_hsc_markers(hsc_marker_yaml)
    hsc_score = _score_hsc(a, markers)
    if hsc_score.max() <= 0:
        raise DCTARootError("HSC marker score is non-positive for all cells")

    cell_comp, n_comp = _connected_components_from_paga(a)
    a.obs["dcta_component"] = cell_comp

    dpt_combined = np.full(a.n_obs, np.nan, dtype=float)
    root_cells: list[int] = []
    root_clusters: list[str] = []

    n_dcs = int(a.uns.get("dcta", {}).get("n_dcs", a.obsm["X_diffmap"].shape[1] - 1))
    dm = a.obsm["X_diffmap"][:, :n_dcs]
    leiden = a.obs["dcta_leiden"].astype(str).to_numpy()

    for comp_id in range(n_comp):
        mask_comp = (cell_comp == comp_id)
        if mask_comp.sum() < 3:
            continue
        comp_clusters = np.unique(leiden[mask_comp])
        cluster_scores = {c: hsc_score[(leiden == c) & mask_comp].mean()
                          for c in comp_clusters}
        hsc_cluster = max(cluster_scores, key=cluster_scores.get)
        cluster_mask = (leiden == hsc_cluster) & mask_comp
        cluster_idx = np.where(cluster_mask)[0]
        medoid_local = _medoid_cell(dm[cluster_idx])
        iroot = int(cluster_idx[medoid_local])
        root_cells.append(iroot)
        root_clusters.append(hsc_cluster)

        dpt_full = _dpt_from_iroot(a, iroot)
        dpt_combined[mask_comp] = dpt_full[mask_comp]

    a.obs["dcta_dpt"] = dpt_combined

    # Robustness: top-K HSC candidates (global), re-run DPT, compute std
    top_k_idx = np.argsort(-hsc_score)[:robustness_k]
    dpts = np.zeros((robustness_k, a.n_obs), dtype=float)
    for i, root in enumerate(top_k_idx):
        dpts[i, :] = _dpt_from_iroot(a, int(root))
    dpt_root_std = np.nanstd(dpts, axis=0)

    # Per-cluster mean HSC marker score (used by dag.build for marker-direction)
    cluster_hsc_score = {
        c: float(hsc_score[leiden == c].mean())
        for c in np.unique(leiden)
    }

    uns = a.uns.setdefault("dcta", {})
    uns["root_strategy"] = root_strategy
    uns["root_cells"] = [int(x) for x in root_cells]
    uns["root_cluster"] = root_clusters[0] if root_clusters else None
    uns["component_count"] = int(n_comp)
    uns["dpt_root_std"] = dpt_root_std
    uns["cluster_hsc_score"] = cluster_hsc_score
    # Phase 1: maturation = inverse of HSC score (no mature-marker YAML yet)
    uns["cluster_maturation_score"] = {c: -v for c, v in cluster_hsc_score.items()}
    return a
