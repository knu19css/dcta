"""Default DCTA result visualization (4 panels per dataset, single PNG).

Treating `dcta_annotation_smoothed` as the canonical DCTA output. Plot
layout (2x2 grid, same paper UMAP coordinates in every panel):

  Panel 1: cluster id only
  Panel 2: author annotation (verbose paper label)
  Panel 3: DCTA annotation (original, pre-smoothing)
  Panel 4: DCTA annotation + smoothed (default output)

Rules:
  - Every panel uses the SAME cluster -> color map (tab20/tab20b cycle by
    cluster index), so a cluster has a consistent color across panels.
  - Each panel shows ONLY cluster numbers on the plot, at cluster
    centroids. adjustText repels overlaps; leader lines stay light grey.
  - Per-panel legend on the right of that panel: `C##  ->  annotation`,
    where annotation is the panel-specific label (paper / DCTA orig /
    DCTA smoothed).

Inputs:
  results/phase1_<ds>/dcta_phase1_smoothed.h5ad
Outputs:
  results/phase1_<ds>/umap_default.png   (overwritten)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parents[1]

DATASETS = [
    ("Zeng2025_BMM", ROOT / "results/phase1_hay_paperclust",
     "AuthorCellType_Broad"),
    ("Triana2021",   ROOT / "results/phase1_triana",
     "paper_AuthorCellType_Broad"),
    ("Strati2023",   ROOT / "results/phase1_strati",
     "paper_AuthorCellType_Broad"),
]


def pick_umap_key(a: ad.AnnData) -> str | None:
    default_emb = a.uns.get("default_embedding")
    if default_emb and default_emb in a.obsm:
        return default_emb
    if "X_umap" in a.obsm:
        return "X_umap"
    return None


def cluster_sort_key(s: str):
    t = s.lstrip("C")
    try:
        return (0, int(t))
    except ValueError:
        return (1, s)


def cluster_mode(values: np.ndarray) -> str:
    """Most common value within a cluster's cells (mode)."""
    if len(values) == 0:
        return "?"
    v, c = np.unique(values, return_counts=True)
    return str(v[np.argmax(c)])


def draw_panel(ax, um, clu_arr, label_arr, color_for, title, panel_legend,
               on_plot_label: str = "cluster"):
    """Scatter UMAP with on-plot labels at cluster centroids.

    color_for : dict cluster_id -> RGB tuple (same across panels)
    panel_legend : if True, draws cluster->label legend on the right.
    on_plot_label : 'cluster' (panel 1) shows just C## ;
                    'both' (panels 2-4) shows "C##\\n<cell type>"
                    at each cluster centroid.
    """
    rng = np.random.default_rng(0)
    order = rng.permutation(len(clu_arr))
    cols = np.array([color_for[c] for c in clu_arr])
    ax.scatter(um[order, 0], um[order, 1], c=cols[order],
               s=2.0, alpha=0.55, linewidths=0)

    clusters = sorted(set(clu_arr), key=cluster_sort_key)
    centroids: dict[str, np.ndarray] = {}
    panel_label: dict[str, str] = {}
    for c in clusters:
        m = clu_arr == c
        if not m.any():
            continue
        centroids[c] = um[m].mean(axis=0)
        panel_label[c] = cluster_mode(label_arr[m]) if label_arr is not None else c

    texts = []
    for c, (cx, cy) in centroids.items():
        if on_plot_label == "cluster" or label_arr is None:
            txt = c
            fs = 6.5
        else:
            # Show "C##\n<cell type>" at the centroid; clamp very long
            # labels (Triana verbose, Strati "C##:..." style) so the on-plot
            # text stays readable.
            lab = panel_label.get(c, "?")
            if len(lab) > 22:
                lab = lab[:20] + "…"
            txt = f"{c}\n{lab}"
            fs = 5.5
        t = ax.text(cx, cy, txt, fontsize=fs, weight="bold",
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.18",
                              facecolor="white", alpha=0.9,
                              edgecolor="grey", linewidth=0.3))
        texts.append(t)
    if texts:
        adjust_text(texts, ax=ax,
                    expand=(1.05, 1.2),
                    arrowprops=dict(arrowstyle="-", color="grey",
                                    lw=0.3, alpha=0.5))

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10)

    if panel_legend:
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="",
                       color=color_for[c], markersize=4,
                       label=f"{c} → {panel_label.get(c, '?')[:30]}")
            for c in clusters
        ]
        ax.legend(handles=handles, loc="center left",
                  bbox_to_anchor=(1.01, 0.5), fontsize=5.5,
                  frameon=False, handlelength=0.8,
                  borderpad=0.2, labelspacing=0.25)


def draw(ds_name: str, outdir: Path, paper_col: str):
    h5 = outdir / "dcta_phase1_smoothed.h5ad"
    if not h5.exists():
        print(f"[skip] {ds_name}: no smoothed h5ad")
        return
    a = ad.read_h5ad(h5)
    emb = pick_umap_key(a)
    if emb is None:
        print(f"[skip] {ds_name}: no UMAP")
        return
    um = np.asarray(a.obsm[emb])
    clu = a.obs["dcta_leiden"].astype(str).to_numpy()
    paper = (a.obs[paper_col].astype(str).to_numpy()
             if paper_col in a.obs.columns else None)
    dcta_orig = (a.obs["dcta_annotation"].astype(str).to_numpy()
                 if "dcta_annotation" in a.obs.columns else None)
    dcta_sm = (a.obs["dcta_annotation_smoothed"].astype(str).to_numpy()
               if "dcta_annotation_smoothed" in a.obs.columns else dcta_orig)

    clusters = sorted(set(clu), key=cluster_sort_key)
    cmap = plt.get_cmap("tab20" if len(clusters) <= 20 else "tab20b")
    color_for = {c: cmap(i % cmap.N) for i, c in enumerate(clusters)}

    fig, axes = plt.subplots(2, 2, figsize=(18, 14),
                             gridspec_kw=dict(wspace=0.7, hspace=0.25))
    draw_panel(axes[0, 0], um, clu, None, color_for,
               title="1) Paper UMAP — cluster id only",
               panel_legend=False, on_plot_label="cluster")
    draw_panel(axes[0, 1], um, clu, paper, color_for,
               title=f"2) Author annotation ({paper_col})",
               panel_legend=True, on_plot_label="both")
    draw_panel(axes[1, 0], um, clu, dcta_orig, color_for,
               title="3) DCTA annotation (original)",
               panel_legend=True, on_plot_label="both")
    draw_panel(axes[1, 1], um, clu, dcta_sm, color_for,
               title="4) DCTA annotation + smoothed (default)",
               panel_legend=True, on_plot_label="both")

    # Share axis limits across panels (use panel 4 as ref)
    xlim = axes[1, 1].get_xlim(); ylim = axes[1, 1].get_ylim()
    for ax in axes.flat:
        ax.set_xlim(xlim); ax.set_ylim(ylim)

    fig.suptitle(f"{ds_name} — DCTA default 4-panel ({emb})",
                 fontsize=12, y=0.995)
    out = outdir / "umap_default.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    for ds, outdir, paper_col in DATASETS:
        draw(ds, outdir, paper_col)


if __name__ == "__main__":
    main()
