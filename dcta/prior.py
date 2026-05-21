"""Hematopoietic prior DAG loading + cluster-to-prior label mapping.

Phase 2.2. Reads `data/priors/dag_prior_bm.yaml` (30 nodes, 31 edges parsed
from the user's DAG_editable.pptx slide 2) and exposes:

- `load_dag_prior(path)` -> networkx.DiGraph with nodes + edges
- `map_clusters_to_prior(adata, ...)` -> attaches cluster->prior labels
  by scoring each cluster against every cell-type panel in the lineage
  marker YAML, then taking argmax. Stores results under
  `adata.uns['dcta']['cluster_prior_label']` and
  `adata.obs['dcta_prior_label']`.
- `prior_path_reachable(G_prior, label_u, label_v, max_hops)` -> bool
- `edge_prior_consistency(G_data, cluster_prior_label, G_prior)` ->
  dict per (u,v): {consistent: bool, reason: str}
"""
from __future__ import annotations
import numpy as np
import networkx as nx
import yaml
import anndata as ad


def load_dag_prior(path: str | list[str]) -> nx.DiGraph:
    """Load DAG prior YAML(s) -> nx.DiGraph.

    Accepts either a single path or a list of paths. Multiple YAMLs are
    merged: union of nodes and edges. Useful for stitching the broad BM
    prior with a maturation-subtree prior (Phase 2.5 T cell tree).

    YAML schema (see data/priors/dag_prior_bm.yaml):
      nodes: list[str]
      edges: list[[parent, child]]
    """
    paths = [path] if isinstance(path, str) else list(path)
    G = nx.DiGraph()
    for p in paths:
        with open(p) as f:
            cfg = yaml.safe_load(f) or {}
        for n in cfg.get("nodes", []) or []:
            G.add_node(str(n))
        for src, dst in cfg.get("edges", []) or []:
            G.add_edge(str(src), str(dst))
    return G


def _score_panel(adata: ad.AnnData, panel: list[str]) -> np.ndarray:
    present = [g for g in panel if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    idx = [adata.var_names.get_loc(g) for g in present]
    X = adata.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X).mean(axis=1)


def map_clusters_to_prior(
    adata: ad.AnnData,
    *,
    lineage_marker_yaml: str,
    g_prior: nx.DiGraph,
    leiden_key: str = "dcta_leiden",
    min_score_quantile: float = 0.0,
) -> ad.AnnData:
    """Attach prior labels to clusters by argmax over ALL cell-type panels
    that map to a node in g_prior.

    For each cluster:
      - score against every cell-type panel whose canonical name OR any
        alias matches a node in g_prior
      - cluster_prior_label = panel with highest cluster-mean score
      - cluster_prior_score = the corresponding score
      - cluster_prior_alternatives = sorted (label, score) top-3 for inspection

    Stores into:
      - adata.obs['dcta_prior_label']  (per-cell, same label for all cells in cluster)
      - adata.uns['dcta']['cluster_prior_label']
      - adata.uns['dcta']['cluster_prior_score']
      - adata.uns['dcta']['cluster_prior_alternatives']
      - adata.uns['dcta']['prior_mapping_unmatched']  (cell types in YAML
        whose canonical name isn't a g_prior node — informational)
    """
    with open(lineage_marker_yaml) as f:
        types = yaml.safe_load(f) or {}

    prior_node_set = set(g_prior.nodes)
    # Build {prior_node: panel_markers}
    # A panel "matches" a prior node if its canonical name == node OR any
    # alias == node.
    panel_for_node: dict[str, list[str]] = {}
    unmatched: list[str] = []
    for tname, tdata in types.items():
        markers = list((tdata or {}).get("markers") or [])
        if not markers:
            continue
        aliases = list((tdata or {}).get("aliases") or [])
        candidates = [tname] + aliases
        hit = next((c for c in candidates if c in prior_node_set), None)
        if hit is None:
            unmatched.append(tname)
            continue
        # If multiple yaml panels map to the same prior node, keep the
        # longest panel (most evidence).
        if hit not in panel_for_node or len(markers) > len(panel_for_node[hit]):
            panel_for_node[hit] = markers

    # Score every panel once per cell
    panel_cell_scores: dict[str, np.ndarray] = {
        node: _score_panel(adata, panel) for node, panel in panel_for_node.items()
    }

    leiden = adata.obs[leiden_key].astype(str).to_numpy()
    cluster_prior_label: dict[str, str | None] = {}
    cluster_prior_score: dict[str, float] = {}
    cluster_prior_alts: dict[str, list] = {}

    for c in np.unique(leiden):
        mask = leiden == c
        if not mask.any():
            continue
        per_node = {
            node: float(scores[mask].mean())
            for node, scores in panel_cell_scores.items()
        }
        if not per_node:
            cluster_prior_label[c] = None
            cluster_prior_score[c] = 0.0
            cluster_prior_alts[c] = []
            continue
        # Optional quantile floor — if best score below floor, label = None
        all_scores = np.array(list(per_node.values()))
        floor = float(np.quantile(all_scores, min_score_quantile)) if min_score_quantile > 0 else -np.inf
        best_node = max(per_node, key=per_node.get)
        if per_node[best_node] < floor:
            cluster_prior_label[c] = None
            cluster_prior_score[c] = float(per_node[best_node])
        else:
            cluster_prior_label[c] = best_node
            cluster_prior_score[c] = float(per_node[best_node])
        # Top-3 alternatives
        cluster_prior_alts[c] = sorted(
            per_node.items(), key=lambda kv: -kv[1]
        )[:3]

    adata.obs["dcta_prior_label"] = [cluster_prior_label.get(c) for c in leiden]
    uns = adata.uns.setdefault("dcta", {})
    uns["cluster_prior_label"] = cluster_prior_label
    uns["cluster_prior_score"] = cluster_prior_score
    uns["cluster_prior_alternatives"] = {
        c: [[lbl, float(s)] for lbl, s in alts]
        for c, alts in cluster_prior_alts.items()
    }
    uns["prior_mapping_unmatched"] = sorted(unmatched)
    return adata


def prior_path_reachable(
    g_prior: nx.DiGraph, u_label: str | None, v_label: str | None,
    max_hops: int | None = None,
) -> bool:
    """Is there a directed path u_label ->* v_label in g_prior (within max_hops)?

    Returns False if either label is None or not a node in g_prior.
    """
    if u_label is None or v_label is None:
        return False
    if u_label == v_label:
        return True
    if u_label not in g_prior or v_label not in g_prior:
        return False
    if max_hops is None:
        return nx.has_path(g_prior, u_label, v_label)
    # Bounded BFS
    try:
        path_len = nx.shortest_path_length(g_prior, u_label, v_label)
        return path_len <= max_hops
    except nx.NetworkXNoPath:
        return False


def edge_prior_consistency(
    g_data: nx.DiGraph,
    cluster_prior_label: dict[str, str | None],
    g_prior: nx.DiGraph,
    *,
    max_hops: int | None = None,
) -> dict:
    """For each edge of g_data, check whether its prior labels are reachable
    in g_prior in the same direction.

    Returns:
      {
        "per_edge": { "u|v": {"u_label","v_label","consistent","reason"} },
        "consistent_count": int,
        "scored_count": int,    # edges where both endpoints have a prior label
        "unlabeled_count": int,
        "consistency_rate": float | None,
      }
    """
    per_edge = {}
    consistent = 0
    scored = 0
    unlabeled = 0
    for u, v in g_data.edges:
        lu = cluster_prior_label.get(u)
        lv = cluster_prior_label.get(v)
        key = f"{u}|{v}"
        if lu is None or lv is None:
            unlabeled += 1
            per_edge[key] = {
                "u_label": lu, "v_label": lv,
                "consistent": None, "reason": "endpoint without prior label",
            }
            continue
        # u_label and v_label may be SAME (cluster labels collapsed) — count as consistent
        if lu == lv:
            consistent += 1
            scored += 1
            per_edge[key] = {
                "u_label": lu, "v_label": lv,
                "consistent": True, "reason": "same prior label",
            }
            continue
        # Forward path lu -> lv in g_prior?
        if prior_path_reachable(g_prior, lu, lv, max_hops=max_hops):
            consistent += 1
            scored += 1
            per_edge[key] = {
                "u_label": lu, "v_label": lv,
                "consistent": True, "reason": "prior-forward reachable",
            }
            continue
        # Reverse path? = G_data edge is anti-prior
        if prior_path_reachable(g_prior, lv, lu, max_hops=max_hops):
            scored += 1
            per_edge[key] = {
                "u_label": lu, "v_label": lv,
                "consistent": False, "reason": "reverse direction (anti-prior)",
            }
        else:
            scored += 1
            per_edge[key] = {
                "u_label": lu, "v_label": lv,
                "consistent": False, "reason": "labels unrelated in prior",
            }
    rate = (consistent / scored) if scored > 0 else None
    return {
        "per_edge": per_edge,
        "consistent_count": consistent,
        "scored_count": scored,
        "unlabeled_count": unlabeled,
        "consistency_rate": rate,
    }
