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
    """Per-cell z-scored mean expression of HSC markers."""
    present = [g for g in markers if g in adata.var_names]
    if not present:
        raise DCTARootError(
            f"None of the HSC marker genes are present in adata.var_names: {markers}"
        )
    return _score_panel(adata, present)


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


def _load_lineage_markers(path: str | None) -> dict:
    """Load lineage_mature_markers.yaml -> {cell_type: {markers, stage, lineage, ...}}."""
    if path is None:
        return {}
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return dict(cfg or {})


def _score_panel(adata: ad.AnnData, panel: list[str]) -> np.ndarray:
    """Per-cell z-scored mean expression of a marker panel; 0.0 if none present.

    Z-scoring per gene before averaging makes the score robust to differences
    in absolute expression magnitude between markers (e.g. IGHG3 ~10 vs CCR7
    ~1 in log-norm). This is the same effect that sc.pp.scale(max_value=10)
    provided when preprocessing was part of the pipeline.
    """
    present = [g for g in panel if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    idx = [adata.var_names.get_loc(g) for g in present]
    X = adata.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    Xz = (X - mu) / sd
    return Xz.mean(axis=1)


def run(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    root_strategy: str = "hsc_marker_medoid",
    hsc_marker_yaml: str | None = None,
    lineage_marker_yaml: str | None = None,
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

    # Phase 2.8: medoid computed on whichever embedding the connectivity
    # stage used (paper X_harmony / X_pca / X_umap), not a local diffusion map.
    use_rep = a.uns.get("dcta", {}).get("use_rep", "X_harmony")
    if use_rep not in a.obsm:
        raise KeyError(
            f"pseudotime: embedding {use_rep!r} not in adata.obsm "
            f"(have {list(a.obsm.keys())})"
        )
    dm = a.obsm[use_rep]
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

    # Phase 2.1: lineage-specific maturation score.
    # If lineage_marker_yaml is provided, for each cluster compute:
    #   - per-cell-type score = mean(panel) per cell, mean over cluster
    #   - cluster_lineage_scores: {cluster: {cell_type: score}}
    #   - cluster_mature_score = max over `mature`/`terminal`-stage cell types
    #     (i.e., "how mature this cluster looks regardless of which lineage")
    #   - cluster_best_lineage = the mature cell type with highest score
    #   - cluster_maturation_score = mature_score - hsc_score
    # If no YAML provided, fall back to Phase 1 behavior (- hsc_score).
    if lineage_marker_yaml is not None:
        types = _load_lineage_markers(lineage_marker_yaml)
        # Pre-compute per-cell score for every cell type with a non-empty panel
        type_cell_scores: dict[str, np.ndarray] = {}
        for tname, tdata in types.items():
            panel = list(tdata.get("markers") or [])
            if not panel:
                continue
            type_cell_scores[tname] = _score_panel(a, panel)

        # Aggregate to cluster level
        cluster_lineage_scores: dict[str, dict[str, float]] = {}
        cluster_mature_score: dict[str, float] = {}
        cluster_best_lineage: dict[str, str | None] = {}
        cluster_maturation: dict[str, float] = {}
        mature_types = [
            t for t, td in types.items()
            if td.get("stage") in ("mature", "terminal")
        ]
        for c in np.unique(leiden):
            mask = leiden == c
            per_type = {
                t: float(s[mask].mean()) for t, s in type_cell_scores.items()
            }
            cluster_lineage_scores[c] = per_type
            mature_only = {t: per_type[t] for t in mature_types if t in per_type}
            if mature_only:
                best = max(mature_only, key=mature_only.get)
                cluster_best_lineage[c] = best
                cluster_mature_score[c] = mature_only[best]
            else:
                cluster_best_lineage[c] = None
                cluster_mature_score[c] = 0.0
            cluster_maturation[c] = (
                cluster_mature_score[c] - cluster_hsc_score.get(c, 0.0)
            )

        uns["cluster_lineage_scores"] = cluster_lineage_scores
        uns["cluster_mature_score"] = cluster_mature_score
        uns["cluster_best_lineage"] = cluster_best_lineage
        uns["cluster_maturation_score"] = cluster_maturation
        uns["lineage_marker_yaml"] = lineage_marker_yaml
    else:
        # Phase 1 fallback
        uns["cluster_maturation_score"] = {
            c: -v for c, v in cluster_hsc_score.items()
        }
        uns["lineage_marker_yaml"] = None
    return a
