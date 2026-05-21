"""DCTA visualization helpers for Phase 1 artifacts."""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import anndata as ad


def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def diffmap_eigvals(adata: ad.AnnData, path: str) -> None:
    eigs = np.asarray(adata.uns["dcta"]["diffmap_eigvals"])
    n_dcs = int(adata.uns["dcta"]["n_dcs"])
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(range(len(eigs)), eigs, marker="o")
    ax.axvline(n_dcs, color="red", linestyle="--", label=f"n_dcs={n_dcs}")
    ax.set_xlabel("DC index")
    ax.set_ylabel("eigenvalue")
    ax.set_title("Diffusion map eigenvalues")
    ax.legend()
    _save(fig, path)


def paga_connectivities(adata: ad.AnnData, path: str) -> None:
    paga = adata.uns["dcta"]["paga"]
    conn = np.asarray(paga["connectivities"])
    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(conn, cmap="viridis")
    ax.set_xticks(range(len(paga["cluster_labels"])))
    ax.set_yticks(range(len(paga["cluster_labels"])))
    ax.set_xticklabels(paga["cluster_labels"], rotation=45)
    ax.set_yticklabels(paga["cluster_labels"])
    ax.set_title("PAGA connectivities")
    fig.colorbar(im, ax=ax)
    _save(fig, path)


def dpt_distribution(adata: ad.AnnData, path: str) -> None:
    dpt = adata.obs["dcta_dpt"].to_numpy()
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(dpt[np.isfinite(dpt)], bins=50, color="steelblue")
    ax.set_xlabel("DPT")
    ax.set_ylabel("cells")
    ax.set_title("DPT pseudotime distribution")
    _save(fig, path)


def dpt_root_std_hist(adata: ad.AnnData, path: str) -> None:
    std = np.asarray(adata.uns["dcta"]["dpt_root_std"])
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.hist(std[np.isfinite(std)], bins=50, color="indianred")
    ax.set_xlabel("DPT std across top-K roots")
    ax.set_ylabel("cells")
    ax.set_title("Root choice robustness")
    _save(fig, path)


def g_data_diagram(adata: ad.AnnData, path: str) -> None:
    G: nx.DiGraph = adata.uns["dcta"]["G_data"]
    # Layout by mean_dpt on x-axis
    pos = {}
    nodes = list(G.nodes)
    for n in nodes:
        x = G.nodes[n].get("mean_dpt", 0.0)
        # spread y by node index hash for stability
        y = (hash(n) % 1000) / 1000.0
        pos[n] = (x, y)
    fig, ax = plt.subplots(figsize=(7, 5))
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="lightblue",
                           node_size=600, edgecolors="black")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9)
    weights = [G[u][v]["weight"] * 4 for u, v in G.edges]
    nx.draw_networkx_edges(G, pos, ax=ax, width=weights, arrows=True,
                           arrowsize=15, edge_color="gray")
    ax.set_xlabel("cluster mean DPT")
    ax.set_title("G_data (DPT-directed DAG)")
    ax.set_yticks([])
    _save(fig, path)


def paga_vs_dag_overlay(adata: ad.AnnData, path: str) -> None:
    """Plot PAGA (undirected, grey) with G_data directed edges overlaid.

    Edge color encodes method-agreement between marker and DPT directions:
      green solid  = both methods agree
      orange dashed = methods disagree (chosen-method direction shown)
      red dotted  = edge removed by FAS cycle removal (drawn from
                    cycle_removal_log)

    Per spec Sec 11.4.
    """
    paga = adata.uns["dcta"]["paga"]
    conn = np.asarray(paga["connectivities"])
    labels = [str(l) for l in paga["cluster_labels"]]
    G_pa = nx.Graph()
    for c in labels:
        G_pa.add_node(c)
    for i, ci in enumerate(labels):
        for j in range(i + 1, len(labels)):
            if conn[i, j] > 0:
                G_pa.add_edge(ci, labels[j], weight=float(conn[i, j]))

    # Layout from PAGA edge weights
    pos = nx.spring_layout(G_pa, seed=0, weight="weight")

    G_dag: nx.DiGraph = adata.uns["dcta"]["G_data"]
    methods = adata.uns["dcta"].get("direction_methods", {})
    dpt_dirs_raw = methods.get("dpt", {}) or {}
    marker_dirs_raw = methods.get("marker", {}) or {}
    def _key(s):
        return tuple(s.split("|", 1)) if "|" in s else (s, "")
    dpt_dirs = {_key(k): tuple(v) for k, v in dpt_dirs_raw.items()}
    marker_dirs = {_key(k): tuple(v) for k, v in marker_dirs_raw.items()}

    fig, ax = plt.subplots(figsize=(8, 6))
    nx.draw_networkx_edges(G_pa, pos, ax=ax, edge_color="lightgrey", width=1.5)
    nx.draw_networkx_nodes(G_dag, pos, ax=ax, node_color="lightblue",
                           node_size=600, edgecolors="black")
    nx.draw_networkx_labels(G_dag, pos, ax=ax, font_size=9)

    # Categorize each retained edge
    green_edges, orange_edges = [], []
    for u, v in G_dag.edges:
        # Determine pair-key (ci, cj) by sorted PAGA order (labels list)
        key = (u, v) if labels.index(u) < labels.index(v) else (v, u)
        d_dir = dpt_dirs.get(key)
        m_dir = marker_dirs.get(key)
        if d_dir == m_dir and d_dir is not None:
            green_edges.append((u, v))
        else:
            orange_edges.append((u, v))

    if green_edges:
        nx.draw_networkx_edges(G_dag, pos, ax=ax, edgelist=green_edges,
                               edge_color="green", arrows=True, arrowsize=15,
                               width=2.0, style="solid")
    if orange_edges:
        nx.draw_networkx_edges(G_dag, pos, ax=ax, edgelist=orange_edges,
                               edge_color="darkorange", arrows=True,
                               arrowsize=15, width=2.0, style="dashed")
    # Removed edges (FAS) drawn as red dotted from log
    removed = adata.uns["dcta"].get("cycle_removal_log", [])
    if removed:
        red_pairs = [(r["u"], r["v"]) for r in removed]
        nx.draw_networkx_edges(G_dag, pos, ax=ax, edgelist=red_pairs,
                               edge_color="red", arrows=True, arrowsize=12,
                               width=1.5, style="dotted")

    chosen = adata.uns["dcta"].get("direction_method_chosen", "?")
    ax.set_title(f"PAGA (grey) <-> G_data overlay - chosen direction: {chosen}\n"
                 f"green=methods agree, orange=disagree, red=removed by FAS")
    ax.set_axis_off()
    _save(fig, path)


def export_all(adata: ad.AnnData, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    diffmap_eigvals(adata, os.path.join(outdir, "diffmap_eigvals.png"))
    paga_connectivities(adata, os.path.join(outdir, "paga_connectivities.png"))
    dpt_distribution(adata, os.path.join(outdir, "dpt_distribution.png"))
    dpt_root_std_hist(adata, os.path.join(outdir, "dpt_root_std_hist.png"))
    g_data_diagram(adata, os.path.join(outdir, "G_data_diagram.png"))
    paga_vs_dag_overlay(adata, os.path.join(outdir, "paga_vs_dag_overlay.png"))
