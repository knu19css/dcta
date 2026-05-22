"""For each dataset, dump per-cluster:
  - DCTA final annotation (dcta_annotation / cluster_annotation)
  - marker-only baseline label (cluster_marker_only_label)
  - top-3 marker-score candidates (recomputed from prior YAML on the same
    cluster centroids on the chosen embedding)
  - disambiguation reason (LCA / parent-child / PAGA tie-break) if DCTA
    deviates from marker-only

Writes:
  docs/dcta_candidates_<ds>.tsv  (TSV per dataset)
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import anndata as ad
import yaml
import csv

ROOT = Path(__file__).resolve().parents[1]


def _zscore_safe(x: np.ndarray) -> np.ndarray:
    m = float(np.nanmean(x))
    s = float(np.nanstd(x))
    return (x - m) / s if s > 0 else (x * 0)


def score_panel(a: ad.AnnData, panel: list[str]) -> np.ndarray:
    """Per-cell z-scored mean (mirrors dcta.prior._score_panel)."""
    var_idx = {g: i for i, g in enumerate(a.var_names)}
    hits = [var_idx[g] for g in panel if g in var_idx]
    if not hits:
        return np.zeros(a.n_obs, dtype=np.float32)
    X = a.X[:, hits]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    # z-score each gene column, then mean across genes
    Z = np.zeros_like(X)
    for j in range(X.shape[1]):
        Z[:, j] = _zscore_safe(X[:, j])
    return Z.mean(axis=1)


def build_panels(marker_yaml: Path) -> dict[str, list[str]]:
    with open(marker_yaml) as f:
        ydoc = yaml.safe_load(f) or {}
    return {n: list((d or {}).get("markers") or []) for n, d in ydoc.items()}


def run_for_dataset(name: str, h5_path: Path, marker_yaml: Path,
                    out_tsv: Path):
    print(f"\n=== {name} ===")
    a = ad.read_h5ad(h5_path)
    panels = build_panels(marker_yaml)
    panels = {n: p for n, p in panels.items() if p}
    print(f"  loaded {a.n_obs} cells, {len(panels)} marker panels")

    # cluster -> bool mask
    leiden = a.obs["dcta_leiden"].astype(str).to_numpy()
    def _ck(x):
        s = x.lstrip("C")
        try:
            return (0, int(s))
        except ValueError:
            return (1, x)
    clusters = sorted(np.unique(leiden), key=_ck)

    # cluster x node score matrix: average per-cell z-mean within cluster
    node_names = list(panels.keys())
    print(f"  scoring {len(node_names)} panels across cluster centroids...")
    cell_scores = {n: score_panel(a, panels[n]) for n in node_names}
    score_mat = np.full((len(clusters), len(node_names)), np.nan, dtype=np.float32)
    for ci, c in enumerate(clusters):
        m = leiden == c
        if not m.any():
            continue
        for ni, n in enumerate(node_names):
            score_mat[ci, ni] = float(cell_scores[n][m].mean())

    d = a.uns.get("dcta", {})
    cluster_annot = d.get("cluster_annotation", {})
    marker_only = d.get("cluster_marker_only_label", {})
    label_map = d.get("cluster_label_map", {})

    fields = ["cluster", "paper_label", "DCTA_decision",
              "top1", "top1_score", "top2", "top2_score",
              "top3", "top3_score",
              "DCTA_vs_top1", "margin"]
    rows = [fields]
    for ci, c in enumerate(clusters):
        order = np.argsort(-np.nan_to_num(score_mat[ci]))
        t1, t2, t3 = (node_names[order[0]] if len(order) > 0 else None,
                      node_names[order[1]] if len(order) > 1 else None,
                      node_names[order[2]] if len(order) > 2 else None)
        s1 = float(score_mat[ci, order[0]]) if len(order) > 0 else float("nan")
        s2 = float(score_mat[ci, order[1]]) if len(order) > 1 else float("nan")
        s3 = float(score_mat[ci, order[2]]) if len(order) > 2 else float("nan")
        decision = cluster_annot.get(c, "")
        paper = label_map.get(c, "")
        margin = (s1 - s2) / max(abs(s1), 1e-9) if not np.isnan(s2) else 0.0
        marker_argmax_diff = "MOVED" if (decision and t1 and decision != t1) else "SAME"
        rows.append([c, paper, decision,
                     t1, f"{s1:+.3f}",
                     t2, f"{s2:+.3f}",
                     t3, f"{s3:+.3f}",
                     marker_argmax_diff, f"{margin:+.2f}"])

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for r in rows:
            w.writerow(r)
    print(f"  wrote {out_tsv}")

    # pretty print
    widths = [max(len(str(r[i])) for r in rows) for i in range(len(fields))]
    for r in rows:
        print("  " + "  ".join(str(c).ljust(w) for c, w in zip(r, widths)))


def main():
    cases = [
        ("Zeng (Hay)", ROOT / "results/phase1_hay_paperclust/dcta_phase1.h5ad",
         ROOT / "docs/dcta_candidates_zeng.tsv"),
        ("Triana",     ROOT / "results/phase1_triana/dcta_phase1.h5ad",
         ROOT / "docs/dcta_candidates_triana.tsv"),
        ("Strati",     ROOT / "results/phase1_strati/dcta_phase1.h5ad",
         ROOT / "docs/dcta_candidates_strati.tsv"),
    ]
    marker_yaml = ROOT / "data/priors/lineage_mature_markers.yaml"
    for name, h5, out in cases:
        if not h5.exists():
            print(f"[skip] {h5} missing")
            continue
        run_for_dataset(name, h5, marker_yaml, out)


if __name__ == "__main__":
    main()
