"""Side-by-side UMAP plots: paper canonical label vs DCTA annotation.

For each dataset (Zeng / Triana / Strati), produces:
  results/phase1_<ds>/umap_paper_vs_dcta.png

Uses paper-original X_umap (from CXG download / paper_AuthorCellType_*
attach step). Falls back to obsm['X_umap'] otherwise.
"""
from __future__ import annotations
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

DATASETS = [
    ("Zeng2025_BMM", ROOT / "results/phase1_hay_paperclust/dcta_phase1.h5ad"),
    ("Triana2021",   ROOT / "results/phase1_triana/dcta_phase1.h5ad"),
    ("Strati2023",   ROOT / "results/phase1_strati/dcta_phase1.h5ad"),
]


def pick_palette(n):
    cmap = plt.get_cmap("tab20" if n <= 20 else "tab20b") if n <= 40 else plt.get_cmap("nipy_spectral")
    return [cmap(i / max(n - 1, 1)) for i in range(n)]


def scatter_by_label(ax, X, labels, title, max_label_len=25):
    uniq = [u for u in sorted(set(labels)) if u and u != "nan"]
    pal = pick_palette(len(uniq))
    color_map = dict(zip(uniq, pal))
    # plot in random order to avoid one class dominating top layer
    rng = np.random.default_rng(0)
    order = rng.permutation(len(labels))
    cols = np.array([color_map.get(l, (0.7, 0.7, 0.7, 0.5)) for l in labels])
    ax.scatter(X[order, 0], X[order, 1], c=cols[order], s=1.0, alpha=0.65,
               linewidths=0)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    # legend (truncated)
    handles = []
    for u in uniq[:40]:
        handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                   color=color_map[u], markersize=4,
                                   label=u[:max_label_len]))
    if len(uniq) > 40:
        handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                   color=(0.5, 0.5, 0.5), markersize=4,
                                   label=f"... (+{len(uniq)-40} more)"))
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              fontsize=6, frameon=False)


def main():
    for ds_name, h5 in DATASETS:
        if not h5.exists():
            print(f"[skip] {h5} not found")
            continue
        print(f"[plot] {ds_name} ← {h5}")
        a = ad.read_h5ad(h5)
        # Prefer dataset's published default_embedding (CXG Explorer default).
        # Triana 2021 has multiple UMAPs; the CXG default (X_mofaumap) is the
        # MOFA-latent UMAP, which is the dataset's canonical view.
        default_emb = a.uns.get("default_embedding")
        emb_key = None
        if default_emb and default_emb in a.obsm:
            emb_key = default_emb
        elif "X_umap" in a.obsm:
            emb_key = "X_umap"
        if emb_key is None:
            print(f"  [warn] no UMAP-like obsm; skipping")
            continue
        print(f"  using {emb_key}")
        X = a.obsm[emb_key]

        # GT label column priority: paper_AuthorCellType_canonical >
        # paper_AuthorCellType_Broad > cell_type
        gt_col = None
        for c in ("paper_AuthorCellType_canonical",
                  "paper_AuthorCellType_Broad", "cell_type"):
            if c in a.obs.columns:
                gt_col = c
                break
        if gt_col is None:
            print(f"  [warn] no GT label column; skipping")
            continue
        gt = a.obs[gt_col].astype(str).to_numpy()
        pred = a.obs["dcta_annotation"].astype(str).to_numpy() \
            if "dcta_annotation" in a.obs.columns else None

        # Layout: 1 row x 2 cols
        fig, axes = plt.subplots(1, 2, figsize=(18, 7),
                                 gridspec_kw=dict(wspace=0.55))
        scatter_by_label(axes[0], X, gt,
                         title=f"{ds_name} — paper label ({gt_col})")
        if pred is not None:
            scatter_by_label(axes[1], X, pred,
                             title=f"{ds_name} — DCTA dcta_annotation")
        else:
            axes[1].text(0.5, 0.5, "no dcta_annotation", ha='center')
            axes[1].set_xticks([]); axes[1].set_yticks([])

        fig.suptitle(f"{ds_name}: paper vs DCTA on paper UMAP", y=0.99,
                     fontsize=12)
        out = h5.parent / "umap_paper_vs_dcta.png"
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
