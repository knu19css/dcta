"""Plot the DCTA backbone path (max-weight longest path) as:
  Left  : straight chain (cluster + paper label + DCTA annotation)
  Right : paper UMAP with backbone clusters highlighted and arrows between
          consecutive backbone cluster centroids

For each dataset, writes results/phase1_<ds>/backbone_path.png.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]

DATASETS = [
    ("Zeng2025_BMM", ROOT / "results/phase1_hay_paperclust"),
    ("Triana2021",   ROOT / "results/phase1_triana"),
    ("Strati2023",   ROOT / "results/phase1_strati"),
]


def pick_umap_key(a: ad.AnnData) -> str | None:
    default_emb = a.uns.get("default_embedding")
    if default_emb and default_emb in a.obsm:
        return default_emb
    if "X_umap" in a.obsm:
        return "X_umap"
    return None


def draw(ds_name: str, outdir: Path):
    h5 = outdir / "dcta_phase1.h5ad"
    js = outdir / "summary.json"
    if not h5.exists() or not js.exists():
        print(f"[skip] {ds_name}: missing artefacts")
        return
    with open(js) as f:
        summary = json.load(f)
    bb = summary["backbone"]
    bb_path = bb["backbone_path"]

    a = ad.read_h5ad(h5)
    d = a.uns.get("dcta", {})
    cluster_annot = d.get("cluster_annotation", {})
    label_map = d.get("cluster_label_map", {})
    sen_map = dict(bb.get("backbone_seniority") or [])

    emb_key = pick_umap_key(a)
    um = np.asarray(a.obsm[emb_key]) if emb_key else None
    leiden = a.obs["dcta_leiden"].astype(str).to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(20, 7),
                             gridspec_kw=dict(width_ratios=[1.0, 1.4]))
    ax_chain, ax_umap = axes

    # ---- LEFT: chain ----
    n = len(bb_path)
    ys = np.arange(n)[::-1]  # top = root
    cmap = plt.get_cmap("viridis")
    cols = [cmap(i / max(n - 1, 1)) for i in range(n)]
    for i, (c, y) in enumerate(zip(bb_path, ys)):
        sen = sen_map.get(c)
        ann = cluster_annot.get(c, "?")
        paper = label_map.get(c, "?")
        ax_chain.scatter([0], [y], s=900, color=cols[i], edgecolor="black",
                         linewidth=1.2, zorder=3)
        ax_chain.text(0, y, c, ha="center", va="center", fontsize=8,
                      weight="bold", zorder=4)
        sen_str = f"sen={sen}" if sen is not None else "sen=?"
        ax_chain.text(0.4, y, f"{ann}", fontsize=10, va="center")
        ax_chain.text(2.5, y, f"({sen_str})", fontsize=8, va="center",
                      color="#555555")
        ax_chain.text(3.4, y, f"paper: {paper}", fontsize=8, va="center",
                      color="#444444", style="italic")
        if i < n - 1:
            ax_chain.annotate("", xy=(0, ys[i + 1] + 0.25),
                              xytext=(0, y - 0.25),
                              arrowprops=dict(arrowstyle="->", color="navy",
                                              lw=2.0))

    ax_chain.set_xlim(-0.8, 7.5)
    ax_chain.set_ylim(-0.7, n - 0.3)
    ax_chain.set_axis_off()
    ax_chain.set_title(
        f"{ds_name} — backbone chain\n"
        f"length={bb['backbone_length']}, "
        f"τ={bb.get('kendall_tau')}, "
        f"monotonic={bb.get('seniority_monotonic')}, "
        f"score={bb['backbone_score']:.2f}",
        fontsize=10, loc="left")

    # ---- RIGHT: UMAP with backbone highlighted ----
    if um is not None:
        # all cells in light grey
        ax_umap.scatter(um[:, 0], um[:, 1], s=2, color="#e0e0e0",
                        alpha=0.6, linewidths=0, zorder=1)
        # backbone clusters colored
        centroids: dict[str, np.ndarray] = {}
        for i, c in enumerate(bb_path):
            m = leiden == c
            if not m.any():
                continue
            ax_umap.scatter(um[m, 0], um[m, 1], s=10, color=cols[i],
                            alpha=0.85, linewidths=0, zorder=2,
                            label=f"{c} {cluster_annot.get(c,'?')}")
            centroids[c] = um[m].mean(axis=0)
        # arrows between consecutive centroids
        for c1, c2 in zip(bb_path[:-1], bb_path[1:]):
            if c1 not in centroids or c2 not in centroids:
                continue
            x0, y0 = centroids[c1]
            x1, y1 = centroids[c2]
            ax_umap.annotate("", xy=(x1, y1), xytext=(x0, y0),
                             arrowprops=dict(arrowstyle="->", color="navy",
                                             lw=2.5, alpha=0.9), zorder=4)
        # bold labels on centroids
        for i, c in enumerate(bb_path):
            if c not in centroids:
                continue
            cx, cy = centroids[c]
            ax_umap.text(cx, cy, c, fontsize=8, ha="center", va="center",
                         weight="bold", zorder=5,
                         bbox=dict(boxstyle="round,pad=0.2",
                                   facecolor="white", alpha=0.9,
                                   edgecolor="grey", linewidth=0.5))
        ax_umap.set_xlabel(f"UMAP1 ({emb_key})")
        ax_umap.set_ylabel(f"UMAP2 ({emb_key})")
        ax_umap.set_title(f"{ds_name} — backbone on paper UMAP",
                          fontsize=10, loc="left")
        ax_umap.set_xticks([]); ax_umap.set_yticks([])
    else:
        ax_umap.text(0.5, 0.5, "no UMAP available",
                     ha="center", va="center", transform=ax_umap.transAxes)
        ax_umap.set_axis_off()

    out = outdir / "backbone_path.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    for ds_name, outdir in DATASETS:
        draw(ds_name, outdir)


if __name__ == "__main__":
    main()
