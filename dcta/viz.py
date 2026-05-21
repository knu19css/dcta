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
    # Layout by mean_dpt on x-axis; deterministic y-jitter via stable hash
    # of the string, then nudge any overlapping pairs apart so nx.draw
    # doesn't emit zero-length arrows (matplotlib bezier crashes on those).
    import hashlib
    pos = {}
    used = set()
    nodes = list(G.nodes)
    for n in nodes:
        x = float(G.nodes[n].get("mean_dpt", 0.0) or 0.0)
        h = int(hashlib.md5(str(n).encode()).hexdigest(), 16) % 10000
        y = h / 10000.0
        # Nudge if collision (same x and y as another node)
        k = 0
        while (round(x, 6), round(y, 6)) in used:
            y = (y + 0.013) % 1.0
            k += 1
            if k > 50:
                break
        used.add((round(x, 6), round(y, 6)))
        pos[n] = (x, y)
    fig, ax = plt.subplots(figsize=(7, 5))
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="lightblue",
                           node_size=600, edgecolors="black")
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=9)
    weights = [max(0.2, G[u][v].get("weight", 0.1) * 4) for u, v in G.edges]
    try:
        nx.draw_networkx_edges(G, pos, ax=ax, width=weights, arrows=True,
                               arrowsize=15, edge_color="gray")
    except (StopIteration, ValueError):
        # Fallback: draw without arrowheads (matplotlib bezier corner case
        # when two nodes overlap in display coords)
        nx.draw_networkx_edges(G, pos, ax=ax, width=weights, arrows=False,
                               edge_color="gray")
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
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    ax_raw, ax = axes
    if "X_umap" not in adata.obsm:
        for a_ in axes:
            a_.text(0.5, 0.5, "No X_umap in adata.obsm\n"
                              "(original paper embedding not available)",
                    ha="center", va="center", fontsize=12, transform=a_.transAxes)
            a_.set_axis_off()
        _save(fig, path)
        return

    um = np.asarray(adata.obsm["X_umap"])
    leiden = adata.obs["dcta_leiden"].astype(str).to_numpy()
    clusters_l = sorted(np.unique(leiden),
                        key=lambda x: int(x) if x.isdigit() else x)
    cmap_l = plt.cm.tab20
    cc_l = {c: cmap_l(i % 20) for i, c in enumerate(clusters_l)}
    for c in clusters_l:
        m = leiden == c
        ax_raw.scatter(um[m, 0], um[m, 1], s=4, color=cc_l[c],
                       alpha=0.6, linewidths=0)
    for c in clusters_l:
        m = leiden == c
        if not m.any():
            continue
        cx, cy = um[m].mean(axis=0)
        ax_raw.text(cx, cy, c, fontsize=7, ha="center", va="center",
                    weight="bold",
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="white", alpha=0.85,
                              edgecolor="grey", linewidth=0.5))
    ax_raw.set_xlabel("UMAP 1 (paper)")
    ax_raw.set_ylabel("UMAP 2 (paper)")
    ax_raw.set_title("Paper UMAP + Leiden clusters", fontsize=10)
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
    # Cluster id labels (+ annotation if available)
    cluster_annot = adata.uns.get("dcta", {}).get("cluster_annotation", {}) or {}
    for c, (cx, cy) in centroids.items():
        ann = cluster_annot.get(c)
        label = f"{c}\n{ann}" if ann else c
        ax.text(cx, cy, label, fontsize=6, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          alpha=0.8, edgecolor="grey"))
    ax.set_xlabel("UMAP 1 (paper embedding)")
    ax.set_ylabel("UMAP 2 (paper embedding)")
    ax.set_title("Paper UMAP with G_data overlay\n"
                 "(black=edges, navy=backbone, cluster id + annotation)",
                 fontsize=10)
    ax_raw.set_xlim(ax.get_xlim())
    ax_raw.set_ylim(ax.get_ylim())
    _save(fig, path)


def _draw_umap_panel(ax, um, labels_arr, title, *, with_text=True,
                     fontsize_text=7):
    """Scatter `um` colored by `labels_arr` (1-D, n_obs strings or None),
    optionally place bold centroid labels per category. Returns dict of
    label -> centroid for axis sharing.
    """
    centroids: dict = {}
    if labels_arr is None:
        ax.scatter(um[:, 0], um[:, 1], s=4, color="#777777",
                   alpha=0.6, linewidths=0)
    else:
        labels_arr = np.asarray(
            [("NA" if (v is None or (isinstance(v, float) and np.isnan(v))
                       or str(v).lower() in ("nan", "none", ""))
              else str(v)) for v in labels_arr]
        )
        uniq = sorted(np.unique(labels_arr),
                      key=lambda x: (x == "NA", int(x) if x.isdigit() else x))
        cmap = plt.cm.tab20
        color_for = {l: ("#cccccc" if l == "NA" else cmap(i % 20))
                     for i, l in enumerate(uniq)}
        for l in uniq:
            m = labels_arr == l
            ax.scatter(um[m, 0], um[m, 1], s=4, color=color_for[l],
                       alpha=0.6, linewidths=0)
            if m.any() and l != "NA":
                centroids[l] = um[m].mean(axis=0)
        if with_text:
            for l, (cx, cy) in centroids.items():
                ax.text(cx, cy, l, fontsize=fontsize_text,
                        ha="center", va="center", weight="bold",
                        bbox=dict(boxstyle="round,pad=0.2",
                                  facecolor="white", alpha=0.85,
                                  edgecolor="grey", linewidth=0.5))
    ax.set_xlabel("UMAP 1 (paper)")
    ax.set_ylabel("UMAP 2 (paper)")
    ax.set_title(title, fontsize=10)
    return centroids


def umap_annotation(adata: ad.AnnData, path: str,
                    col: str = "dcta_annotation",
                    paper_label_col: str = "paper_AuthorCellType",
                    cluster_col: str = "dcta_leiden") -> None:
    """Three-panel side-by-side comparison on the SAME paper UMAP:
        1. cluster id only         (left)   — `cluster_col`
        2. paper cell-type label   (middle) — `paper_label_col`
        3. DCTA predicted label    (right)  — `col`

    Any panel whose source column is missing falls back to a notice. All
    three panels share x/y limits so points are pixel-aligned.
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    ax_clu, ax_paper, ax_pred = axes
    if "X_umap" not in adata.obsm:
        for a_ in axes:
            a_.text(0.5, 0.5, "missing X_umap",
                    ha="center", va="center", transform=a_.transAxes)
            a_.set_axis_off()
        _save(fig, path)
        return
    um = np.asarray(adata.obsm["X_umap"])

    # ---- Panel 1: cluster id ----
    clu = adata.obs[cluster_col].astype(object).to_numpy() if cluster_col in adata.obs.columns else None
    _draw_umap_panel(
        ax_clu, um, clu,
        title=(f"1) Paper UMAP + cluster id ({cluster_col})"
               if clu is not None else "1) cluster column missing"),
    )

    # ---- Panel 2: paper annotation ----
    paper = adata.obs[paper_label_col].astype(object).to_numpy() if paper_label_col in adata.obs.columns else None
    _draw_umap_panel(
        ax_paper, um, paper,
        title=(f"2) Paper annotation ({paper_label_col})"
               if paper is not None else "2) paper label column missing"),
    )

    # ---- Panel 3: DCTA prediction ----
    pred = adata.obs[col].astype(object).to_numpy() if col in adata.obs.columns else None
    _draw_umap_panel(
        ax_pred, um, pred,
        title=(f"3) DCTA prediction ({col})"
               if pred is not None else f"3) {col} missing"),
    )

    # Share axes (same UMAP coords)
    xlim = ax_paper.get_xlim() if paper is not None else ax_clu.get_xlim()
    ylim = ax_paper.get_ylim() if paper is not None else ax_clu.get_ylim()
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
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


def umap_marker_vs_dcta(adata: ad.AnnData, path: str,
                        cluster_col: str = "dcta_leiden") -> None:
    """Three-panel UMAP demonstrating DAG impact on annotation:
      1) cluster id (`cluster_col`)
      2) marker-only baseline (no DAG) — `dcta_marker_only_annotation`
      3) DCTA full (DAG-aware disambiguation) — `dcta_annotation`

    Each cell's color/label differs between panels 2 and 3 only when
    the DAG-aware logic changed the call.
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    ax_clu, ax_mo, ax_full = axes
    if ("X_umap" not in adata.obsm
        or "dcta_marker_only_annotation" not in adata.obs.columns
        or "dcta_annotation" not in adata.obs.columns):
        for a_ in axes:
            a_.text(0.5, 0.5,
                    "marker-only baseline or DCTA annotation missing",
                    ha="center", va="center", transform=a_.transAxes)
            a_.set_axis_off()
        _save(fig, path)
        return
    um = np.asarray(adata.obsm["X_umap"])
    clu = adata.obs[cluster_col].astype(object).to_numpy() if cluster_col in adata.obs.columns else None
    mo = adata.obs["dcta_marker_only_annotation"].astype(object).to_numpy()
    full = adata.obs["dcta_annotation"].astype(object).to_numpy()
    _draw_umap_panel(ax_clu, um, clu,
                     title=(f"1) Paper UMAP + cluster id ({cluster_col})"
                            if clu is not None else "1) cluster missing"))
    _draw_umap_panel(ax_mo, um, mo,
                     title="2) Marker-only baseline (no DAG)")
    _draw_umap_panel(ax_full, um, full,
                     title="3) DCTA full (DAG-aware disambiguation)")
    xlim = ax_full.get_xlim(); ylim = ax_full.get_ylim()
    for ax in axes:
        ax.set_xlim(xlim); ax.set_ylim(ylim)
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
    if "dcta_marker_only_annotation" in adata.obs.columns:
        umap_marker_vs_dcta(
            adata, os.path.join(outdir, "umap_marker_vs_dcta.png"))
