"""DCTA G_data construction.

For each candidate PAGA edge above `edge_threshold`, compute TWO directions:
  - "dpt":    low cluster-mean DPT -> high cluster-mean DPT
  - "marker": high cluster HSC-marker score -> low (immature -> mature)

If a label column with seniority mapping is available, score each method's
agreement with biology-derived expected direction and pick the winner. Else
default to "marker". Build G_data from the chosen method's directions, then
FAS-remove cycles.
"""
from __future__ import annotations
import numpy as np
import networkx as nx
import anndata as ad
import yaml
from dcta import DCTAGraphError
from dcta import defaults as D


def _greedy_fas_remove(G: nx.DiGraph) -> list[tuple[str, str, dict]]:
    """Greedy feedback arc set: while G has a cycle, remove the weakest edge
    on a found cycle. Returns the list of removed edges with their data.
    Modifies G in place.
    """
    removed: list[tuple[str, str, dict]] = []
    safety = G.number_of_edges() + 10
    for _ in range(safety):
        try:
            cyc = nx.find_cycle(G, orientation="original")
        except nx.NetworkXNoCycle:
            return removed
        weakest = min(cyc, key=lambda e: G[e[0]][e[1]].get("weight", 0.0))
        u, v = weakest[0], weakest[1]
        data = dict(G[u][v])
        removed.append((u, v, data))
        G.remove_edge(u, v)
    raise DCTAGraphError("FAS greedy did not converge")


def _load_seniority(path: str | None) -> dict[str, int]:
    if path is None:
        return {}
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return dict(cfg.get("seniority", {}))


def _cluster_dominant_label(adata: ad.AnnData, label_column: str | None) -> dict[str, str]:
    """Return {cluster_id: dominant_label} using mode of label_column per cluster."""
    if label_column is None or label_column not in adata.obs.columns:
        return {}
    leiden = adata.obs["dcta_leiden"].astype(str).to_numpy()
    labels = adata.obs[label_column].astype(str).to_numpy()
    out: dict[str, str] = {}
    for c in np.unique(leiden):
        mask = leiden == c
        if not mask.any():
            continue
        vals, counts = np.unique(labels[mask], return_counts=True)
        out[str(c)] = str(vals[int(np.argmax(counts))])
    return out


def _direction_dpt(cluster_mean_dpt: dict[str, float],
                   ci: str, cj: str) -> tuple[str, str] | None:
    """Return (source, target) for DPT direction, or None if undecidable."""
    di, dj = cluster_mean_dpt.get(ci), cluster_mean_dpt.get(cj)
    if di is None or dj is None or np.isnan(di) or np.isnan(dj):
        return None
    return (ci, cj) if di <= dj else (cj, ci)


def _direction_marker(cluster_hsc_score: dict[str, float],
                      ci: str, cj: str) -> tuple[str, str] | None:
    """Return (source, target) for marker direction.

    Higher HSC score = more immature = source.
    """
    hi, hj = cluster_hsc_score.get(ci), cluster_hsc_score.get(cj)
    if hi is None or hj is None or np.isnan(hi) or np.isnan(hj):
        return None
    return (ci, cj) if hi >= hj else (cj, ci)


def _score_method(direction_map: dict[tuple[str, str], tuple[str, str]],
                  cluster_seniority: dict[str, int]) -> float | None:
    """Fraction of edges where method's direction agrees with seniority.

    Only edges whose BOTH endpoints have a seniority and seniorities differ
    contribute. Returns None if no scoreable edges.
    """
    correct, total = 0, 0
    for (ci, cj), (src, dst) in direction_map.items():
        sa = cluster_seniority.get(src)
        sb = cluster_seniority.get(dst)
        if sa is None or sb is None or sa == sb:
            continue
        total += 1
        if sa < sb:  # expected: source has lower seniority
            correct += 1
    return None if total == 0 else correct / total


def build(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    edge_threshold: float = D.EDGE_THRESHOLD,
    fas_method: str = D.FAS_METHOD,
    transitive_reduction: bool = D.TRANSITIVE_REDUCTION,
    direction_strategy: str = "auto_select",
    seniority_yaml: str | None = None,
    label_column: str | None = "cell_type",
) -> ad.AnnData:
    """Build G_data from PAGA connectivities with dual-direction selection.

    Direction strategy:
      - "auto_select" (default): compute both `dpt` and `marker` directions,
        score each against seniority labels if available, choose higher.
        Falls back to "marker" when labels missing.
      - "marker": force marker direction.
      - "dpt": force DPT direction.
    """
    a = adata.copy() if copy else adata
    if fas_method != "greedy":
        raise NotImplementedError(f"fas_method={fas_method!r} not implemented")
    if direction_strategy not in ("auto_select", "marker", "dpt"):
        raise ValueError(f"direction_strategy={direction_strategy!r} not in "
                         "('auto_select', 'marker', 'dpt')")

    paga = a.uns["dcta"]["paga"]
    conn = np.asarray(paga["connectivities"])
    labels = [str(l) for l in paga["cluster_labels"]]
    dpt = a.obs["dcta_dpt"].to_numpy()
    leiden = a.obs["dcta_leiden"].astype(str).to_numpy()
    cluster_mean_dpt = {c: float(np.nanmean(dpt[leiden == c])) for c in labels}
    cluster_hsc_score = dict(a.uns["dcta"]["cluster_hsc_score"])

    # Compute both direction methods per candidate edge
    candidate_pairs: list[tuple[str, str, float]] = []
    n = len(labels)
    for i in range(n):
        for j in range(i + 1, n):
            w = float(conn[i, j])
            if w < edge_threshold:
                continue
            candidate_pairs.append((labels[i], labels[j], w))

    dpt_dirs: dict[tuple[str, str], tuple[str, str]] = {}
    marker_dirs: dict[tuple[str, str], tuple[str, str]] = {}
    for ci, cj, _w in candidate_pairs:
        d = _direction_dpt(cluster_mean_dpt, ci, cj)
        m = _direction_marker(cluster_hsc_score, ci, cj)
        if d is not None:
            dpt_dirs[(ci, cj)] = d
        if m is not None:
            marker_dirs[(ci, cj)] = m

    # Per-pair agreement
    pairs_in_both = set(dpt_dirs) & set(marker_dirs)
    agreement = sum(1 for k in pairs_in_both if dpt_dirs[k] == marker_dirs[k])
    agreement_rate = (agreement / len(pairs_in_both)) if pairs_in_both else None

    # Score methods against seniority if labels available
    seniority_map = _load_seniority(seniority_yaml)
    cluster_labels_map = _cluster_dominant_label(a, label_column)
    cluster_seniority = {
        c: seniority_map.get(lbl) for c, lbl in cluster_labels_map.items()
        if seniority_map.get(lbl) is not None
    }
    if cluster_seniority:
        acc_dpt = _score_method(dpt_dirs, cluster_seniority)
        acc_marker = _score_method(marker_dirs, cluster_seniority)
        accuracy = {"dpt": acc_dpt, "marker": acc_marker}
    else:
        accuracy = {}

    # Choose method
    if direction_strategy in ("marker", "dpt"):
        chosen = direction_strategy
    elif accuracy and accuracy.get("dpt") is not None and accuracy.get("marker") is not None:
        chosen = "marker" if accuracy["marker"] >= accuracy["dpt"] else "dpt"
    else:
        chosen = "marker"  # default: biology-grounded

    chosen_dirs = marker_dirs if chosen == "marker" else dpt_dirs

    # Build G with chosen directions
    G = nx.DiGraph()
    for c in labels:
        G.add_node(c, mean_dpt=cluster_mean_dpt.get(c, float("nan")),
                   hsc_score=cluster_hsc_score.get(c, float("nan")),
                   seniority=cluster_seniority.get(c))
    pair_weight = {(ci, cj): w for ci, cj, w in candidate_pairs}
    for (ci, cj), (src, dst) in chosen_dirs.items():
        w = pair_weight[(ci, cj)]
        G.add_edge(src, dst, weight=w)

    removed = _greedy_fas_remove(G)

    if transitive_reduction:
        attrs_backup = {n_: dict(G.nodes[n_]) for n_ in G.nodes}
        G = nx.transitive_reduction(G)
        for n_, data in attrs_backup.items():
            if n_ in G.nodes:
                G.nodes[n_].update(data)

    uns = a.uns.setdefault("dcta", {})
    uns["G_data"] = G
    uns["cycle_removal_log"] = [
        {"u": u, "v": v, "weight": d.get("weight", 0.0)} for u, v, d in removed
    ]
    uns["direction_methods"] = {
        "dpt": {f"{k[0]}|{k[1]}": list(v) for k, v in dpt_dirs.items()},
        "marker": {f"{k[0]}|{k[1]}": list(v) for k, v in marker_dirs.items()},
        "agreement_rate": agreement_rate,
        "accuracy": accuracy,
        "n_pairs_compared": len(pairs_in_both),
    }
    uns["direction_method_chosen"] = chosen
    uns["biological_seniority_map"] = cluster_seniority

    # Backbone-vs-hematopoiesis validation (per user feedback 2026-05-21)
    uns["backbone_validation"] = _validate_backbone(G, cluster_seniority)
    return a


def _longest_path_length(G: nx.DiGraph, source: str) -> tuple[list[str], int]:
    """Return (longest_path_from_source, length) in a DAG. Empty list if no
    reachable nodes."""
    if source not in G.nodes:
        return [], 0
    # Among descendants of source, find the topological-order path with max length
    desc = nx.descendants(G, source) | {source}
    sub = G.subgraph(desc)
    if sub.number_of_nodes() == 0:
        return [], 0
    # Longest path in DAG via DP on topo order
    topo = list(nx.topological_sort(sub))
    dist = {n: 0 for n in topo}
    parent = {n: None for n in topo}
    for u in topo:
        for v in sub.successors(u):
            if dist[u] + 1 > dist[v]:
                dist[v] = dist[u] + 1
                parent[v] = u
    end = max(dist, key=dist.get)
    path = []
    cur = end
    while cur is not None:
        path.append(cur)
        cur = parent[cur]
    path.reverse()
    return path, dist[end]


def _validate_backbone(G: nx.DiGraph,
                       cluster_seniority: dict[str, int]) -> dict:
    """Check whether G_data's backbone chain matches expected hematopoietic order.

    The "backbone" is defined as the longest directed path starting from the
    HSC root cluster (= node whose seniority is the minimum, defaults to the
    first node if no seniority known).

    Returns a dict suitable for adata.uns['dcta']['backbone_validation'].

    Metrics:
      - backbone_path: list of cluster ids on the chain
      - backbone_length: number of edges in the chain
      - backbone_seniority: list of (cluster, seniority) along the path
      - seniority_monotonic: True if seniority is non-decreasing along path
      - seniority_strictly_increasing: True if strictly increasing (no ties)
      - kendall_tau: Kendall's tau between topological order on path
        and seniority order (only for clusters with known seniority)
      - n_lineages_reachable_from_root: count of distinct mature labels reachable
      - n_disconnected_mature: count of mature clusters NOT reached from root
      - root_cluster: which node was used as start
    """
    out: dict = {
        "backbone_path": [],
        "backbone_length": 0,
        "backbone_seniority": [],
        "seniority_monotonic": None,
        "seniority_strictly_increasing": None,
        "kendall_tau": None,
        "n_lineages_reachable_from_root": 0,
        "n_disconnected_mature": 0,
        "root_cluster": None,
    }
    if G.number_of_nodes() == 0:
        return out

    # Root = node with minimum known seniority. Tie-break: lexicographic.
    if cluster_seniority:
        min_sen = min(cluster_seniority.values())
        roots = sorted(c for c, s in cluster_seniority.items() if s == min_sen)
        root = roots[0] if roots else None
    else:
        root = None
    if root is None:
        # Fallback: node with no incoming edges
        sources = [n for n in G.nodes if G.in_degree(n) == 0]
        root = sources[0] if sources else next(iter(G.nodes))
    out["root_cluster"] = root

    path, length = _longest_path_length(G, root)
    out["backbone_path"] = path
    out["backbone_length"] = length
    out["backbone_seniority"] = [(n, cluster_seniority.get(n)) for n in path]

    # Monotonicity over clusters whose seniority is known
    known = [s for _, s in out["backbone_seniority"] if s is not None]
    if len(known) >= 2:
        out["seniority_monotonic"] = all(
            known[i] <= known[i + 1] for i in range(len(known) - 1)
        )
        out["seniority_strictly_increasing"] = all(
            known[i] < known[i + 1] for i in range(len(known) - 1)
        )
        # Kendall's tau between path order index and seniority value
        try:
            from scipy.stats import kendalltau
            tau, _ = kendalltau(list(range(len(known))), known)
            out["kendall_tau"] = float(tau) if tau == tau else None  # NaN guard
        except Exception:
            out["kendall_tau"] = None

    # Mature-cluster reachability (seniority >= 4 = mature)
    # NOTE: Disconnected mature clusters are NOT necessarily errors — many BM
    # mature lymphoid cells (T/NK/mature B) are *tissue-resident*, having
    # arrived from thymus/periphery rather than differentiated in this BM
    # sample. The expectation is that erythroid + myeloid maturation chains
    # connect to HSC root; lymphoid lineages may legitimately disconnect.
    mature = {c for c, s in cluster_seniority.items() if s is not None and s >= 4}
    if mature:
        reachable = nx.descendants(G, root) | {root}
        reachable_mature = mature & reachable
        out["n_lineages_reachable_from_root"] = len(reachable_mature)
        out["n_disconnected_mature"] = len(mature - reachable)
        out["disconnected_mature_clusters"] = sorted(mature - reachable)

    return out
