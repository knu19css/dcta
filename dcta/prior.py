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
    """Per-cell z-scored mean of a marker panel (same as pseudotime._score_panel)."""
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


def _is_progenitor(g_prior: nx.DiGraph, node: str) -> bool:
    """A node is a progenitor in g_prior if it has at least one successor."""
    return node in g_prior and g_prior.out_degree(node) > 0


def _lca_in_dag(g_prior: nx.DiGraph, a: str, b: str,
                max_lca_hops: int = 2) -> str | None:
    """Lowest (closest) common ancestor of a and b in the directed g_prior.

    Returns the common ancestor whose MAX shortest-path-length to (a, b) is
    SMALLEST -- i.e. the most recent shared progenitor. Returns None if:
      - either node is not in g_prior
      - no common ancestor exists
      - the closest LCA is more than `max_lca_hops` away from BOTH a and b
        (= too far up to be a meaningful "shared progenitor" of these
        terminals; usually it's HSC which is uninformative)
    """
    if a not in g_prior or b not in g_prior:
        return None
    anc_a = nx.ancestors(g_prior, a) | {a}
    anc_b = nx.ancestors(g_prior, b) | {b}
    common = anc_a & anc_b
    if not common:
        return None

    def far(n):
        try:
            return max(nx.shortest_path_length(g_prior, n, a),
                       nx.shortest_path_length(g_prior, n, b))
        except nx.NetworkXNoPath:
            return 10**9
    best = min(common, key=far)
    if far(best) > max_lca_hops:
        return None
    return best


def map_clusters_to_prior(
    adata: ad.AnnData,
    *,
    lineage_marker_yaml: str,
    g_prior: nx.DiGraph,
    leiden_key: str = "dcta_leiden",
    min_score_quantile: float = 0.0,
    prior_fusion_weight: float = 0.0,
    prior_fusion_hops: int = 1,
    ambiguity_margin: float = 0.15,
    max_lca_hops: int = 2,
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
    node_stage: dict[str, str] = {}
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
        if hit not in panel_for_node or len(markers) > len(panel_for_node[hit]):
            panel_for_node[hit] = markers
        st = (tdata or {}).get("stage")
        if st and hit not in node_stage:
            node_stage[hit] = str(st)

    # Score every panel once per cell
    panel_cell_scores: dict[str, np.ndarray] = {
        node: _score_panel(adata, panel) for node, panel in panel_for_node.items()
    }

    leiden = adata.obs[leiden_key].astype(str).to_numpy()

    # Compute per-cluster x per-prior-node MARKER score (cluster mean)
    cluster_marker_score: dict[str, dict[str, float]] = {}
    for c in np.unique(leiden):
        mask = leiden == c
        if not mask.any():
            continue
        cluster_marker_score[c] = {
            node: float(scores[mask].mean())
            for node, scores in panel_cell_scores.items()
        }

    # -------- Marker-only baseline (DAG ignored) --------
    cluster_marker_label: dict[str, str | None] = {}
    cluster_marker_alts: dict[str, list] = {}
    for c, per_node in cluster_marker_score.items():
        if not per_node:
            cluster_marker_label[c] = None
            cluster_marker_alts[c] = []
            continue
        all_scores = np.array(list(per_node.values()))
        floor = float(np.quantile(all_scores, min_score_quantile)) if min_score_quantile > 0 else -np.inf
        best_node = max(per_node, key=per_node.get)
        cluster_marker_label[c] = best_node if per_node[best_node] >= floor else None
        cluster_marker_alts[c] = sorted(per_node.items(), key=lambda kv: -kv[1])[:3]

    # -------- Prior-fused score (Phase 2.6) --------
    # For each cluster c and each candidate prior node `label`:
    #   neighbor_support(c, label) = mean over G_prior in/out neighbors n of
    #     label, of mean marker_score of `n` over PAGA-adjacent clusters of c.
    #
    # final_score(c, label) =
    #   (1 - alpha) * marker_score(c, label) + alpha * neighbor_support(c, label)
    cluster_fused_label: dict[str, str | None] = {}
    cluster_fused_alts: dict[str, list] = {}
    cluster_fused_score: dict[str, dict[str, float]] = {}

    alpha = float(prior_fusion_weight)
    hops = max(1, int(prior_fusion_hops))
    # PAGA neighbor structure (cluster->cluster)
    paga = adata.uns.get("dcta", {}).get("paga", {})
    paga_labels = list(paga.get("cluster_labels", []))
    paga_conn = paga.get("connectivities")
    paga_idx = {str(c): i for i, c in enumerate(paga_labels)}
    if paga_conn is not None and len(paga_labels) > 0:
        import numpy as _np
        conn_arr = _np.asarray(paga_conn)
    else:
        conn_arr = None

    def paga_neighbors(c):
        if conn_arr is None or c not in paga_idx:
            return []
        i = paga_idx[c]
        return [paga_labels[j] for j in range(len(paga_labels))
                if j != i and conn_arr[i, j] > 0]

    def prior_neighbors(label, hops):
        if label not in prior_node_set:
            return set()
        out = set()
        # k-hop in- and out-neighbors (excluding self)
        frontier = {label}
        for _ in range(hops):
            new = set()
            for n in frontier:
                if n in g_prior:
                    new.update(g_prior.predecessors(n))
                    new.update(g_prior.successors(n))
            frontier = new - out - {label}
            out |= frontier
        return out

    # -------- Per-cluster disambiguation reasons (for trace/inspection) -----
    cluster_disambig_reason: dict[str, str] = {}

    for c, per_node in cluster_marker_score.items():
        if not per_node:
            cluster_fused_label[c] = None
            cluster_fused_score[c] = {}
            cluster_fused_alts[c] = []
            cluster_disambig_reason[c] = "no panel scored"
            continue
        # Sort by marker score descending
        sorted_labels = sorted(per_node.items(), key=lambda kv: -kv[1])
        best_label, best_score = sorted_labels[0]
        runner_label, runner_score = (sorted_labels[1] if len(sorted_labels) > 1
                                      else (None, 0.0))
        margin = (best_score - runner_score) / max(best_score, 1e-9)

        # Default: marker argmax wins
        chosen = best_label
        reason = f"marker argmax (margin {margin:.3f})"
        fused = dict(per_node)

        # AMBIGUOUS only when margin small AND runner_label exists
        if (runner_label is not None and margin < ambiguity_margin
                and best_label in g_prior and runner_label in g_prior):
            best_is_prog = _is_progenitor(g_prior, best_label)
            runner_is_prog = _is_progenitor(g_prior, runner_label)
            both_terminal = (not best_is_prog) and (not runner_is_prog)
            both_progenitor = best_is_prog and runner_is_prog
            # Phase 2.8.1: if BOTH candidates are mature/terminal (per
            # yaml `stage`), skip LCA-upshift even when both have
            # out-degree>0 in G_prior. Picking a thymic LCA (e.g. DP) as
            # the label for a mature peripheral cluster is wrong.
            best_stage = node_stage.get(best_label, "")
            runner_stage = node_stage.get(runner_label, "")
            both_mature = (best_stage in ("mature", "terminal")
                           and runner_stage in ("mature", "terminal"))

            if both_terminal or both_mature:
                # === Case B: both terminals -> PAGA neighbor tie-break ===
                # Which candidate's G_prior neighborhood is more enriched
                # in cluster c's PAGA-adjacent clusters' marker scores?
                # Example: NK vs CD8_effector — if neighbor clusters
                # express CD4, that supports CD8 (T family); if not, NK.
                nbrs_c = paga_neighbors(c)
                support = {}
                for cand in (best_label, runner_label):
                    p_nbrs = prior_neighbors(cand, hops)
                    vals = []
                    for nbr_c in nbrs_c:
                        if nbr_c not in cluster_marker_score:
                            continue
                        nbr_scores = cluster_marker_score[nbr_c]
                        v = [nbr_scores.get(pn, 0.0) for pn in p_nbrs if pn in nbr_scores]
                        if v:
                            vals.append(float(np.mean(v)))
                    support[cand] = float(np.mean(vals)) if vals else 0.0
                if support.get(runner_label, 0) > support.get(best_label, 0):
                    chosen = runner_label
                    reason = (f"terminal-tiebreak: runner {runner_label} "
                              f"(neighbor support {support[runner_label]:.3f} "
                              f"> best {support[best_label]:.3f})")
                else:
                    chosen = best_label
                    reason = (f"terminal-tiebreak: best {best_label} "
                              f"(neighbor support {support.get(best_label,0):.3f} "
                              f">= runner {support.get(runner_label,0):.3f})")
            else:
                # === Case A: at least one is a progenitor -> consider LCA ===
                # If the closest LCA is within max_lca_hops, use it. Else
                # fall back to whichever of the top-2 IS the progenitor
                # (multipotent expression of itself + one of its descendants).
                lca = _lca_in_dag(g_prior, best_label, runner_label,
                                  max_lca_hops=max_lca_hops)
                if lca is not None and lca not in (best_label, runner_label):
                    if _is_progenitor(g_prior, lca):
                        chosen = lca
                        reason = (f"LCA-progenitor ({lca}) of "
                                  f"({best_label}, {runner_label}) within "
                                  f"{max_lca_hops} hops")
                    else:
                        chosen = best_label
                        reason = (f"LCA ({lca}) is not progenitor; keep best")
                elif lca in (best_label, runner_label):
                    # One IS the parent of the other -> use the parent
                    chosen = lca
                    reason = (f"parent-child top-2; use parent ({lca})")
                elif both_progenitor:
                    # Both progenitors but no close LCA -> prefer the one
                    # higher up (more in-degree-0 or fewer ancestors)
                    cand_anc = {
                        best_label: len(nx.ancestors(g_prior, best_label)),
                        runner_label: len(nx.ancestors(g_prior, runner_label)),
                    }
                    chosen = min(cand_anc, key=cand_anc.get)
                    reason = (f"both progenitors no close LCA; "
                              f"choose more-upstream ({chosen})")
                else:
                    # One progenitor + one terminal, no close LCA -> keep best
                    chosen = best_label
                    reason = f"mixed progenitor+terminal, no close LCA; keep best"

        # Optional global fusion (legacy): blend marker scores w/ neighbor
        # support — kept for backward compat when alpha > 0
        if alpha > 0:
            nbrs_c = paga_neighbors(c)
            blended = {}
            for label, m_score in per_node.items():
                p_nbrs = prior_neighbors(label, hops)
                if not p_nbrs or not nbrs_c:
                    blended[label] = m_score
                    continue
                neighbor_supports = []
                for nbr_c in nbrs_c:
                    if nbr_c not in cluster_marker_score:
                        continue
                    nbr_scores = cluster_marker_score[nbr_c]
                    vals = [nbr_scores.get(pn, 0.0) for pn in p_nbrs if pn in nbr_scores]
                    if vals:
                        neighbor_supports.append(float(np.mean(vals)))
                if neighbor_supports:
                    blended[label] = (
                        (1 - alpha) * m_score
                        + alpha * float(np.mean(neighbor_supports))
                    )
                else:
                    blended[label] = m_score
            fused = blended
            # Only override `chosen` if marker margin was not ambiguous
            if margin >= ambiguity_margin:
                chosen = max(fused, key=fused.get)
                reason = f"fused argmax (alpha={alpha})"

        cluster_fused_score[c] = fused
        cluster_fused_label[c] = chosen
        cluster_fused_alts[c] = sorted(fused.items(), key=lambda kv: -kv[1])[:3]
        cluster_disambig_reason[c] = reason

    # The "official" cluster_prior_label = fused result (or marker-only when
    # alpha=0). Both are kept side-by-side for inspection.
    cluster_prior_label = cluster_fused_label
    cluster_prior_score = {
        c: (max(cluster_fused_score[c].values()) if cluster_fused_score.get(c) else 0.0)
        for c in cluster_marker_score
    }

    adata.obs["dcta_prior_label"] = [cluster_prior_label.get(c) for c in leiden]
    uns = adata.uns.setdefault("dcta", {})
    uns["cluster_prior_label"] = cluster_prior_label
    uns["cluster_prior_score"] = cluster_prior_score
    # Baseline (marker-only) labels for comparison
    uns["cluster_marker_only_label"] = cluster_marker_label
    uns["cluster_marker_only_alternatives"] = {
        c: [[lbl, float(s)] for lbl, s in alts]
        for c, alts in cluster_marker_alts.items()
    }
    uns["cluster_prior_alternatives"] = {
        c: [[lbl, float(s)] for lbl, s in alts]
        for c, alts in cluster_fused_alts.items()
    }
    uns["prior_fusion_weight"] = alpha
    uns["prior_fusion_hops"] = hops
    uns["ambiguity_margin"] = ambiguity_margin
    uns["cluster_disambig_reason"] = cluster_disambig_reason
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
