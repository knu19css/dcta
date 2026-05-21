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

    # Phase 2.2 — prior-inconsistent edges drawn purple (on top of green/orange)
    pc = adata.uns["dcta"].get("prior_consistency", {}) or {}
    pe = (pc or {}).get("per_edge", {}) or {}
    purple_edges = []
    for u, v in G_dag.edges:
        info = pe.get(f"{u}|{v}", {})
        if info.get("consistent") is False:
            purple_edges.append((u, v))
    if purple_edges:
        nx.draw_networkx_edges(G_dag, pos, ax=ax, edgelist=purple_edges,
                               edge_color="purple", arrows=True, arrowsize=15,
                               width=2.2, style="solid")

    # Backbone path: highlight in thick navy on top of everything
    bb = adata.uns["dcta"].get("backbone_validation", {})
    bb_path = bb.get("backbone_path", []) if bb else []
    if len(bb_path) >= 2:
        bb_edges = list(zip(bb_path[:-1], bb_path[1:]))
        bb_edges = [(u, v) for (u, v) in bb_edges if G_dag.has_edge(u, v)]
        if bb_edges:
            nx.draw_networkx_edges(G_dag, pos, ax=ax, edgelist=bb_edges,
                                   edge_color="navy", arrows=True, arrowsize=18,
                                   width=4.0, style="solid", alpha=0.7)

    chosen = adata.uns["dcta"].get("direction_method_chosen", "?")
    mono = bb.get("seniority_monotonic") if bb else None
    tau = bb.get("kendall_tau") if bb else None
    pcr = pc.get("consistency_rate") if pc else None
    title = (f"PAGA (grey) <-> G_data overlay - chosen={chosen}\n"
             f"green=agree, orange=disagree, red=FAS-removed, "
             f"purple=anti-prior, navy=backbone "
             f"(monotonic={mono}, tau={tau}, prior_consistency={pcr})")
    ax.set_title(title, fontsize=9)
    ax.set_axis_off()
    _save(fig, path)


def umap_paper_with_g_data(adata: ad.AnnData, path: str) -> None:
    """Plot the dataset's ORIGINAL paper UMAP (adata.obsm['X_umap']) colored
    by leiden cluster, with G_data backbone overlaid as arrows between
    cluster centroids.

    If no X_umap exists in adata.obsm, draws a notice-only placeholder.
    Per feedback_use_original_embedding: prefer published embedding over
    locally-recomputed UMAP.
    """
    fig, ax = plt.subplots(figsize=(8, 7))
    if "X_umap" not in adata.obsm:
        ax.text(0.5, 0.5, "No X_umap in adata.obsm\n"
                          "(original paper embedding not available)",
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        _save(fig, path)
        return

    um = np.asarray(adata.obsm["X_umap"])
    leiden = adata.obs["dcta_leiden"].astype(str).to_numpy()
    clusters = sorted(np.unique(leiden), key=lambda x: int(x) if x.isdigit() else x)
    cmap = plt.cm.tab20
    color_for = {c: cmap(i % 20) for i, c in enumerate(clusters)}
    for c in clusters:
        m = (leiden == c)
        ax.scatter(um[m, 0], um[m, 1], s=4, color=color_for[c],
                   label=c, alpha=0.6, linewidths=0)
    # Centroids in UMAP
    centroids = {c: um[(leiden == c)].mean(axis=0) for c in clusters}
    G: nx.DiGraph = adata.uns["dcta"]["G_data"]
    for u, v in G.edges:
        if u not in centroids or v not in centroids:
            continue
        x0, y0 = centroids[u]
        x1, y1 = centroids[v]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", color="black",
                                    lw=1.2, alpha=0.7))
    # Backbone in navy
    bb = adata.uns["dcta"].get("backbone_validation", {})
    bb_path = bb.get("backbone_path", []) if bb else []
    if len(bb_path) >= 2:
        for u, v in zip(bb_path[:-1], bb_path[1:]):
            if u not in centroids or v not in centroids:
                continue
            x0, y0 = centroids[u]
            x1, y1 = centroids[v]
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="->", color="navy",
                                        lw=3.0, alpha=0.85))
    # Cluster id labels
    for c, (cx, cy) in centroids.items():
        ax.text(cx, cy, c, fontsize=8, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          alpha=0.7, edgecolor="grey"))
    ax.set_xlabel("UMAP 1 (paper embedding)")
    ax.set_ylabel("UMAP 2 (paper embedding)")
    ax.set_title("Original paper UMAP with G_data overlay\n"
                 "(black=edges, navy=backbone)", fontsize=10)
    _save(fig, path)


def umap_annotation(adata: ad.AnnData, path: str,
                    col: str = "dcta_annotation") -> None:
    """Paper UMAP colored by annotation label (or placeholder if missing)."""
    fig, ax = plt.subplots(figsize=(8, 7))
    if "X_umap" not in adata.obsm or col not in adata.obs.columns:
        ax.text(0.5, 0.5, f"missing {col!r} or X_umap",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        _save(fig, path)
        return
    um = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs[col].astype(object).fillna("NA").to_numpy()
    uniq = sorted(set(labels))
    cmap = plt.cm.tab20
    color_for = {l: ("#cccccc" if l == "NA" else cmap(i % 20))
                 for i, l in enumerate(uniq)}
    for l in uniq:
        m = labels == l
        ax.scatter(um[m, 0], um[m, 1], s=4, color=color_for[l],
                   label=l, alpha=0.7, linewidths=0)
    ax.legend(markerscale=2, fontsize=6, loc="best", ncol=2,
              frameon=False, labelspacing=0.2)
    ax.set_xlabel("UMAP 1 (paper)")
    ax.set_ylabel("UMAP 2 (paper)")
    ax.set_title(f"Paper UMAP colored by {col}", fontsize=10)
    _save(fig, path)


def umap_annotation_confidence(adata: ad.AnnData, path: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 7))
    if "X_umap" not in adata.obsm or "dcta_annotation_confidence" not in adata.obs.columns:
        ax.text(0.5, 0.5, "missing X_umap or dcta_annotation_confidence",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        _save(fig, path)
        return
    um = np.asarray(adata.obsm["X_umap"])
    conf = np.asarray(adata.obs["dcta_annotation_confidence"], dtype=float)
    sc = ax.scatter(um[:, 0], um[:, 1], s=4, c=conf, cmap="viridis",
                    vmin=0.0, vmax=1.0, alpha=0.7, linewidths=0)
    fig.colorbar(sc, ax=ax, label="absorption confidence")
    ax.set_xlabel("UMAP 1 (paper)")
    ax.set_ylabel("UMAP 2 (paper)")
    ax.set_title("Paper UMAP colored by annotation confidence", fontsize=10)
    _save(fig, path)


def export_all(adata: ad.AnnData, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    diffmap_eigvals(adata, os.path.join(outdir, "diffmap_eigvals.png"))
    paga_connectivities(adata, os.path.join(outdir, "paga_connectivities.png"))
    dpt_distribution(adata, os.path.join(outdir, "dpt_distribution.png"))
    dpt_root_std_hist(adata, os.path.join(outdir, "dpt_root_std_hist.png"))
    g_data_diagram(adata, os.path.join(outdir, "G_data_diagram.png"))
    paga_vs_dag_overlay(adata, os.path.join(outdir, "paga_vs_dag_overlay.png"))
    umap_paper_with_g_data(adata, os.path.join(outdir, "umap_paper_with_g_data.png"))
    if "dcta_annotation" in adata.obs.columns:
        umap_annotation(adata, os.path.join(outdir, "umap_annotation.png"))
        umap_annotation_confidence(
            adata, os.path.join(outdir, "umap_annotation_confidence.png"))
