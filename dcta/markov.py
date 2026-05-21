"""Absorbing Markov chain on G_prior for DCTA Phase 2.3 annotation.

Approach:
- Identify absorbing (terminal) states = G_prior nodes with out_degree == 0.
  All other G_prior nodes are transient (progenitors, intermediate states).
- Build a transition probability matrix on G_prior nodes using PAGA-cluster
  connectivity propagated through cluster_prior_label assignments:
    edge (n1 -> n2) weight = sum of PAGA connectivities between clusters
    mapped to n1 and n2, restricted to edges that exist in G_prior
    (forward direction). Plus a self-loop for absorbing states.
- Row-normalize to get transition probabilities P. Partition into:
    P = [[Q, R], [0, I]]
  Q = transient->transient (t x t), R = transient->absorbing (t x a).
- Absorption probability matrix B = (I - Q)^{-1} R, shape (t, a).
- For each leiden cluster c with prior label n_c:
    if n_c is absorbing => annotation = n_c, confidence = 1.0
    if n_c is transient => annotation = absorbing state argmax B[n_c, :],
                            confidence = max B[n_c, :]
- Cells inherit cluster-level annotation (Phase 2.3 simplification).
"""
from __future__ import annotations
import numpy as np
import networkx as nx
import anndata as ad


def _build_transition_matrix(
    g_prior: nx.DiGraph,
    cluster_prior_label: dict[str, str | None],
    paga_connectivities: np.ndarray,
    paga_cluster_labels: list[str],
) -> tuple[np.ndarray, list[str], list[str], list[str]]:
    """Build row-normalized transition matrix on G_prior nodes.

    Returns:
      P                : (N x N) row-stochastic transition matrix
      node_order       : list of N prior nodes (transient first, then absorbing)
      transient_nodes  : list of transient nodes (out_degree > 0)
      absorbing_nodes  : list of absorbing nodes (out_degree == 0)

    Transition probabilities aggregate PAGA connectivity between clusters
    mapped to each prior node, restricted to G_prior edges. Absorbing
    states have a 1.0 self-loop.
    """
    nodes = list(g_prior.nodes)
    absorbing = sorted([n for n in nodes if g_prior.out_degree(n) == 0])
    transient = sorted([n for n in nodes if g_prior.out_degree(n) > 0])
    node_order = transient + absorbing
    N = len(node_order)
    node_idx = {n: i for i, n in enumerate(node_order)}
    P = np.zeros((N, N), dtype=float)

    # Aggregate PAGA conn between prior nodes
    # Build cluster->prior_node lookup, ignoring None / unmapped clusters
    cl_to_node = {c: n for c, n in cluster_prior_label.items() if n in node_idx}

    # For each G_prior edge (n_from -> n_to), find total PAGA connectivity
    # between cluster sets mapped to n_from and n_to.
    cluster_idx = {c: i for i, c in enumerate(paga_cluster_labels)}
    for n_from, n_to in g_prior.edges:
        clusters_from = [c for c, n in cl_to_node.items() if n == n_from]
        clusters_to = [c for c, n in cl_to_node.items() if n == n_to]
        if not clusters_from or not clusters_to:
            continue
        idx_from = [cluster_idx[c] for c in clusters_from if c in cluster_idx]
        idx_to = [cluster_idx[c] for c in clusters_to if c in cluster_idx]
        if not idx_from or not idx_to:
            continue
        w = float(paga_connectivities[np.ix_(idx_from, idx_to)].sum())
        if w <= 0:
            continue
        P[node_idx[n_from], node_idx[n_to]] += w

    # Add a small self-loop probability for transient nodes so they don't
    # become rows of zeros (which break the chain). If a transient node has
    # no observed outgoing PAGA mass, it stays put with prob 1 (= effectively
    # treats it as absorbing for this dataset).
    for n in transient:
        i = node_idx[n]
        row_sum = P[i].sum()
        if row_sum <= 0:
            P[i, i] = 1.0  # no evidence => treat as absorbing self-loop
    # Absorbing rows: identity self-loop
    for n in absorbing:
        i = node_idx[n]
        P[i, i] = 1.0

    # Row-normalize transient rows
    for n in transient:
        i = node_idx[n]
        s = P[i].sum()
        if s > 0:
            P[i] = P[i] / s

    return P, node_order, transient, absorbing


def absorption_probabilities(
    g_prior: nx.DiGraph,
    cluster_prior_label: dict[str, str | None],
    paga_connectivities: np.ndarray,
    paga_cluster_labels: list[str],
    damping: float = 0.0,
) -> dict:
    """Compute absorption probability matrix B for the G_prior chain.

    damping (0..1) shrinks transient self-loops toward 0 and reallocates the
    mass uniformly across reachable absorbing states. Used to escape near-
    singular (I - Q) matrices when a transient node has no observed PAGA
    mass flowing to any successor (entire row was synthetic self-loop).

    Returns dict with P, B, node_order, transient_nodes, absorbing_nodes,
    Q, R, N_fund.
    """
    P, node_order, transient, absorbing = _build_transition_matrix(
        g_prior, cluster_prior_label, paga_connectivities, paga_cluster_labels
    )
    t = len(transient)
    a = len(absorbing)

    # Identify transient rows that are pure self-loop (no observed mass to
    # other states). These cause (I - Q) to be singular. If absorbing
    # states exist, redirect each such row to a uniform distribution over
    # the absorbing states reachable from that transient node in g_prior.
    node_to_idx = {n: i for i, n in enumerate(node_order)}
    for i, n in enumerate(transient):
        if P[i, i] >= 1.0 - 1e-12:
            reachable_abs = [
                m for m in absorbing if nx.has_path(g_prior, n, m)
            ]
            if reachable_abs:
                P[i] = 0.0
                for m in reachable_abs:
                    P[i, node_to_idx[m]] = 1.0 / len(reachable_abs)
            # else: leave as self-loop; transient is effectively a sink

    # Optional damping: pull a fraction of self-loop mass into absorbing
    # uniform (safety net even when row isn't strictly self-loop=1).
    if damping > 0:
        for i, n in enumerate(transient):
            if P[i, i] > 0:
                escape = damping * P[i, i]
                reachable_abs = [m for m in absorbing if nx.has_path(g_prior, n, m)]
                if reachable_abs:
                    P[i, i] -= escape
                    per = escape / len(reachable_abs)
                    for m in reachable_abs:
                        P[i, node_to_idx[m]] += per

    Q = P[:t, :t]
    R = P[:t, t:]
    try:
        I_Q = np.eye(t) - Q
        N_fund = np.linalg.inv(I_Q)
    except np.linalg.LinAlgError:
        N_fund = np.linalg.pinv(np.eye(t) - Q)
    B = N_fund @ R
    # Numerical sanity: clip tiny negatives, normalize rows that sum > 0
    B = np.where(np.isfinite(B), B, 0.0)
    B = np.clip(B, 0.0, None)
    row_sum = B.sum(axis=1, keepdims=True)
    nz = (row_sum.flatten() > 0)
    B[nz] = B[nz] / row_sum[nz]
    return {
        "P": P, "B": B, "Q": Q, "R": R, "N_fund": N_fund,
        "node_order": node_order,
        "transient_nodes": transient,
        "absorbing_nodes": absorbing,
    }


def annotate_clusters(
    adata: ad.AnnData,
    *,
    dag_prior_yaml: str | None = None,
    leiden_key: str = "dcta_leiden",
    copy: bool = False,
) -> ad.AnnData:
    """Cluster-level annotation via absorbing chain on G_prior.

    Requires (Phase 2.2 outputs):
      - adata.uns['dcta']['cluster_prior_label']
      - adata.uns['dcta']['paga'] with 'connectivities' and 'cluster_labels'
      - either dag_prior_yaml param OR existing prior_consistency['g_prior_path']

    Adds to adata.obs:
      - 'dcta_annotation' : terminal lineage label (or absorbing self-label
        if cluster's prior label is already absorbing). NaN if cluster has
        no prior label.
      - 'dcta_annotation_confidence' : max absorption probability (or 1.0
        for absorbing-label clusters).
      - 'dcta_progenitor_state' : if transient, the prior label itself
        (so the user can read 'CMP', 'MEP', etc.). NaN if absorbing.

    Adds to adata.uns['dcta']:
      - 'cluster_annotation'      : {cluster: terminal_label_or_self}
      - 'cluster_annotation_conf' : {cluster: float}
      - 'absorption_matrix'       : dict with B + node lists (JSON-safe shape)
    """
    from dcta import prior as dcta_prior
    a = adata.copy() if copy else adata

    # Resolve G_prior
    if dag_prior_yaml is None:
        pc = a.uns.get("dcta", {}).get("prior_consistency", {}) or {}
        dag_prior_yaml = pc.get("g_prior_path")
    if dag_prior_yaml is None:
        raise ValueError(
            "annotate_clusters requires dag_prior_yaml (or a prior previously "
            "loaded by dcta.dag.build with dag_prior_yaml=...)"
        )
    g_prior = dcta_prior.load_dag_prior(dag_prior_yaml)

    cluster_prior_label = dict(a.uns["dcta"]["cluster_prior_label"])
    paga = a.uns["dcta"]["paga"]
    conn = np.asarray(paga["connectivities"])
    paga_labels = [str(c) for c in paga["cluster_labels"]]

    out = absorption_probabilities(g_prior, cluster_prior_label, conn, paga_labels)
    B = out["B"]
    transient = out["transient_nodes"]
    absorbing = out["absorbing_nodes"]

    # Phase 2.3 semantics (CLAUDE.md sec 45-46):
    #   - annotation = the prior label itself (HSC, MPP, CMP, GMP, Monocyte,
    #     Erythrocyte, ...). Both progenitor (transient) and mature (absorbing)
    #     clusters keep their actual identity.
    #   - confidence = absorption max probability for transient clusters,
    #     1.0 for absorbing clusters, NaN if unmapped.
    #   - terminal_prediction = "what terminal would this transient
    #     differentiate to" = absorbing argmax. Useful for progenitor
    #     lineage-fate prediction; stored separately from annotation.
    cluster_annot: dict[str, str | None] = {}
    cluster_conf: dict[str, float] = {}
    cluster_prog: dict[str, str | None] = {}
    cluster_terminal_pred: dict[str, str | None] = {}
    transient_idx = {n: i for i, n in enumerate(transient)}
    for c, prior_node in cluster_prior_label.items():
        if prior_node is None:
            cluster_annot[c] = None
            cluster_conf[c] = float("nan")
            cluster_prog[c] = None
            cluster_terminal_pred[c] = None
        elif prior_node in absorbing:
            cluster_annot[c] = prior_node
            cluster_conf[c] = 1.0
            cluster_prog[c] = None
            cluster_terminal_pred[c] = prior_node
        elif prior_node in transient_idx:
            # ANNOTATION = the progenitor state itself (e.g. "CMP", "GMP")
            cluster_annot[c] = prior_node
            cluster_prog[c] = prior_node
            row = B[transient_idx[prior_node]]
            if len(absorbing) == 0 or row.sum() == 0:
                cluster_terminal_pred[c] = None
                cluster_conf[c] = float("nan")
            else:
                cluster_terminal_pred[c] = absorbing[int(np.argmax(row))]
                cluster_conf[c] = float(row.max())
        else:
            cluster_annot[c] = None
            cluster_conf[c] = float("nan")
            cluster_prog[c] = None
            cluster_terminal_pred[c] = None

    # Broadcast to cells
    leiden = a.obs[leiden_key].astype(str).to_numpy()
    a.obs["dcta_annotation"] = [cluster_annot.get(c) for c in leiden]
    a.obs["dcta_annotation_confidence"] = [cluster_conf.get(c, float("nan")) for c in leiden]
    a.obs["dcta_progenitor_state"] = [cluster_prog.get(c) for c in leiden]
    a.obs["dcta_terminal_prediction"] = [cluster_terminal_pred.get(c) for c in leiden]

    uns = a.uns.setdefault("dcta", {})
    uns["cluster_annotation"] = cluster_annot
    uns["cluster_annotation_conf"] = cluster_conf
    uns["cluster_progenitor_state"] = cluster_prog
    uns["cluster_terminal_prediction"] = cluster_terminal_pred
    uns["absorption_matrix"] = {
        "B": B.tolist(),
        "transient_nodes": transient,
        "absorbing_nodes": absorbing,
    }
    return a
