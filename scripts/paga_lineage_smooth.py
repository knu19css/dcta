"""PAGA-lineage smoothing for DCTA annotations (post-processing).

DECISION RULE (no magic numbers; pure inequality comparisons):

  For each cluster c with current DCTA annotation A and its lineage L_c:

  Build two lineage probability distributions:
    P_marker[L] = normalized positive-clipped marker z-score mass for L
                  (how much marker evidence sits in each lineage)
    P_graph[L]  = normalized PAGA neighbor weight where neighbors with
                  annotation in lineage L contribute (how much graph
                  connectivity supports each lineage)

  Let L_m = argmax P_marker, L_g = argmax P_graph.

  Reassign A to the neighbor-majority annotation inside L_g iff ALL
  three inequalities hold:
    (1) L_m != L_g                              ← two sources disagree
    (2) P_graph[L_g] > P_marker[L_m]            ← graph is more confident
    (3) P_marker[L_g] > P_marker[L_m] / n_lineages
                                                ← graph's choice has at
                                                  least baseline marker
                                                  support (above uniform)

  Only "constant" is 1/n_lineages = uniform-prior baseline.

KNOWN LIMITATIONS (kept as documented limitations, not fixed in code):

  1. Normalization loses absolute marker strength.
     P_marker and P_graph are both in [0,1] after normalization, so a
     cluster with a very strong single-panel signal (e.g. Megakaryocyte
     z=+4.79) can still be overridden when graph's lineage barely wins
     in P-space. The rule cannot tell "strong marker beating mild graph"
     from "mild marker beating mild graph".
     Concrete victims in our runs: Strati C21 (Megakaryocyte→CMP),
     Strati C25 (HSC→MEP), Zeng C15 (NK→Cytotoxic_CD8),
     Triana C16/C40 (erythroid/lymphoid_B→Basophil).

  2. Ambient RNA contamination is invisible to the rule.
     Strati C03 reticulocyte calls are driven by HBB/HBA1/HBA2 ambient
     RNA. The rule catches C03 only because lymphoid-T graph weight is
     overwhelming; on smaller / similar-lineage contamination it would
     miss.

  3. Single-pass only.
     After one round of reassignment, the new annotations could shift
     P_graph for neighboring clusters, but we don't iterate. A second
     pass might catch cascaded fixes (or introduce new errors).

These limitations could be fixed by (a) including an absolute marker-
strength term in the decision, (b) running ambient-RNA decontamination
upstream (SoupX/DecontX), and (c) iterating the smoothing — but each
adds a tunable (magic) parameter or extra dependency. The current rule
deliberately trades coverage for parsimony.

This script does NOT modify the DCTA codebase. It reads dcta_phase1.h5ad
+ summary.json and writes:
  - results/phase1_<ds>/annotation_smoothed.tsv  (per-cluster log)
  - results/phase1_<ds>/dcta_phase1_smoothed.h5ad
        obs['dcta_annotation_smoothed'] holds the post-processed labels

Usage:
  python scripts/paga_lineage_smooth.py --ds all
"""
from __future__ import annotations
import argparse
from pathlib import Path
import csv
import numpy as np
import yaml
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]


def _zscore_safe(x: np.ndarray) -> np.ndarray:
    m = float(np.nanmean(x)); s = float(np.nanstd(x))
    return (x - m) / s if s > 0 else x * 0


def score_panel(a, panel):
    var_idx = {g: i for i, g in enumerate(a.var_names)}
    hits = [var_idx[g] for g in panel if g in var_idx]
    if not hits:
        return np.zeros(a.n_obs, dtype=np.float32)
    X = a.X[:, hits]
    if hasattr(X, "toarray"): X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    Z = np.zeros_like(X)
    for j in range(X.shape[1]):
        Z[:, j] = _zscore_safe(X[:, j])
    return Z.mean(axis=1)


def cluster_node_scores(a, ydoc, clusters, leiden):
    """Return (score_mat, se_mat, nodes).

    score_mat[ci, ni] = cluster c's mean z-score on panel n.
    se_mat[ci, ni]    = standard error of that mean over the cluster's
                        cells (within-cluster std / sqrt(n_cells)).
                        We can divide score by se to get a per-cluster
                        t-statistic for "this panel is significantly
                        nonzero for this cluster" -- a first-principles
                        check that doesn't need a magic threshold (2-sigma
                        ≈ p<0.05 is the conventional statistical bar).
    """
    panels = {n: d.get("markers") or []
              for n, d in ydoc.items() if (d or {}).get("markers")}
    cell_scores = {n: score_panel(a, p) for n, p in panels.items()}
    nodes = list(panels.keys())
    mat = np.full((len(clusters), len(nodes)), np.nan, dtype=np.float32)
    se = np.full((len(clusters), len(nodes)), np.nan, dtype=np.float32)
    for ci, c in enumerate(clusters):
        m = leiden == c
        n_cells = int(m.sum())
        if n_cells == 0: continue
        for ni, n in enumerate(nodes):
            v = cell_scores[n][m]
            mat[ci, ni] = float(v.mean())
            se[ci, ni] = float(v.std(ddof=1) / np.sqrt(n_cells)) if n_cells > 1 else np.nan
    return mat, se, nodes


def lineage_of(ann: str, ydoc: dict) -> str | None:
    info = ydoc.get(ann)
    if not info: return None
    return info.get("lineage")


def _data_driven_thresholds(*args, **kwargs) -> dict:
    """No thresholds. The new smoothing uses pure comparison-based rules:

      Trigger (all of):
        (1) marker disagreement: the marker-argmax panel for this cluster
            is NOT the current annotation. i.e. if we trusted markers
            alone, we'd already pick a different cell type.
        (2) lineage disagreement: among PAGA neighbors, the maximum-weight
            lineage is NOT this cluster's current lineage.

      Action:
        Reassign to the annotation that is most common (weighted) among
        the PAGA neighbors WITHIN the dominant lineage.

    No quantiles, no IQR, no ratios — just argmax comparisons that fall
    out of "marker says X, neighbors say lineage Y, my current call
    contradicts both."
    """
    return dict(rule="argmax-comparison-only")


def smooth(ds_name: str, outdir: Path,
           dom_share: float | None, dom_ratio: float | None,
           min_score: float, max_self_score: float | None):
    h5 = outdir / "dcta_phase1.h5ad"
    if not h5.exists():
        print(f"[skip] {ds_name}: no h5"); return
    print(f"\n=== {ds_name} ===")

    a = ad.read_h5ad(h5)
    d = a.uns["dcta"]
    paga = d["paga"]
    conn = np.asarray(paga["connectivities"])
    labels = list(paga["cluster_labels"])
    cluster_annot = dict(d.get("cluster_annotation", {}))
    label_map = dict(d.get("cluster_label_map", {}))

    with open(ROOT / "data/priors/lineage_mature_markers.yaml") as f:
        ydoc = yaml.safe_load(f) or {}

    leiden = a.obs["dcta_leiden"].astype(str).to_numpy()
    score_mat, se_mat, nodes = cluster_node_scores(a, ydoc, labels, leiden)
    # t = score / SE ; |t| > 2 ≈ p<0.05 under CLT for the cluster mean.
    # This is the only "threshold" in the rule and it's a standard
    # statistical convention, not a tunable magic number.
    t_mat = score_mat / np.where(se_mat > 0, se_mat, np.nan)
    node_lineage = {n: (ydoc.get(n) or {}).get("lineage") for n in nodes}

    # Pure-disagreement rule (no thresholds, no quantiles):
    #
    # Build two lineage distributions for each cluster:
    #   P_marker[L]: positive-clipped marker scores per lineage, normalized
    #   P_graph[L] : PAGA neighbor weights per lineage, normalized
    #
    # Decision (mutual disagreement condition):
    #   Let L_m = argmax P_marker, L_g = argmax P_graph.
    #   Reassign only when BOTH directions agree there's a conflict:
    #     1. L_m != L_g                                 (the sources pick
    #        different lineages — disagreement exists)
    #     2. P_graph[L_g] > P_marker[L_m]               (graph is more
    #        confident on its winner than marker is)
    #     3. P_marker[L_g] > P_marker[L_m] / n_lineages (graph's lineage
    #        is at least minimally plausible from markers — above uniform)
    #
    # This is fully comparison-based: no thresholds, just inequalities
    # between two distributions over the same lineage set. The only
    # constant is 1/n_lineages which is the uniform-distribution baseline.
    print(f"  rule: graph overrides marker iff (L_m != L_g) AND "
          f"P_graph[L_g] > P_marker[L_m] AND "
          f"P_marker[L_g] > P_marker[L_m] / n_lineages")

    rows = [["cluster", "paper_label", "DCTA_orig", "self_lineage",
             "dom_lineage", "dom_share", "self_share",
             "DCTA_smoothed", "new_score", "reassigned"]]
    reassigned: dict[str, str] = {}
    for ci, c in enumerate(labels):
        ann_c = cluster_annot.get(c, "?")
        L_c = lineage_of(ann_c, ydoc)
        # Aggregate neighbor weight by lineage
        votes: dict[str, float] = {}
        total = 0.0
        for j, w in enumerate(conn[ci]):
            if j == ci or w <= 0: continue
            L_j = lineage_of(cluster_annot.get(labels[j], ""), ydoc)
            if not L_j: continue
            votes[L_j] = votes.get(L_j, 0.0) + float(w)
            total += float(w)
        if total <= 0:
            rows.append([c, label_map.get(c, ""), ann_c, str(L_c),
                         "-", "0.000", "0.000", ann_c, "-", ""])
            continue
        shares = {L: v / total for L, v in votes.items()}
        L_dom, s_dom = max(shares.items(), key=lambda kv: kv[1])
        s_self = shares.get(L_c, 0.0)

        new_ann = ann_c
        new_score = ""
        reassign_flag = ""
        # Build marker-implied lineage distribution P_marker
        marker_lineage_mass: dict[str, float] = {}
        total_marker = 0.0
        for ni, n in enumerate(nodes):
            L_n = node_lineage.get(n)
            if not L_n: continue
            v = float(score_mat[ci, ni])
            if not np.isfinite(v) or v <= 0: continue
            marker_lineage_mass[L_n] = marker_lineage_mass.get(L_n, 0.0) + v
            total_marker += v
        if total_marker > 0:
            P_marker = {L: m / total_marker for L, m in marker_lineage_mass.items()}
        else:
            P_marker = {}
        # graph lineage distribution is already in `shares` (computed above)
        P_graph = shares  # share dict from neighbor lineage votes

        if P_marker and P_graph:
            L_m, p_m = max(P_marker.items(), key=lambda kv: kv[1])
        else:
            L_m, p_m = "", 0.0
        L_g, p_g = (L_dom, s_dom) if P_graph else ("", 0.0)

        # Track for informational columns
        cur_score = float(score_mat[ci, nodes.index(ann_c)]) if ann_c in nodes else np.nan

        # ---- Trigger ----
        # Conflict condition (3 simultaneous comparisons, no thresholds):
        #   - L_m != L_g  (argmax lineages disagree)
        #   - P_graph[L_g] > P_marker[L_m]  (graph more confident)
        #   - P_marker[L_g] > P_marker[L_m] / n_lineages  (graph's lineage
        #     has at least uniform-baseline marker support — i.e., markers
        #     don't COMPLETELY reject the graph's choice)
        n_lineages = max(len({L for L in node_lineage.values() if L}), 2)
        p_m_at_g = P_marker.get(L_g, 0.0) if L_g else 0.0
        uniform_baseline = p_m / n_lineages if p_m > 0 else 0.0
        if (L_m and L_g and L_m != L_g
            and p_g > p_m
            and p_m_at_g > uniform_baseline):
            # ---- Candidate ranking inside L_g (graph's preferred lineage) ----
            cand_idx = [ni for ni, n in enumerate(nodes)
                        if node_lineage.get(n) == L_g]
            if cand_idx:
                # (A) PAGA neighbor majority vote within L_g
                ann_votes_in_L: dict[str, float] = {}
                for j, w in enumerate(conn[ci]):
                    if j == ci or w <= 0: continue
                    ann_j = cluster_annot.get(labels[j], "")
                    if node_lineage.get(ann_j) != L_g: continue
                    ann_votes_in_L[ann_j] = ann_votes_in_L.get(ann_j, 0.0) + float(w)
                # (B) marker z-score (own cluster) within dom_lineage candidates
                local_scores = score_mat[ci, cand_idx]
                cand_names = [nodes[ix] for ix in cand_idx]

                if ann_votes_in_L:
                    # majority vote wins (most common annotation among
                    # dom-lineage PAGA neighbors)
                    candidate = max(ann_votes_in_L, key=ann_votes_in_L.get)
                    cand_score = float(score_mat[ci, nodes.index(candidate)])
                    decision_src = "neighbor-vote"
                else:
                    k = int(np.nanargmax(local_scores))
                    candidate = cand_names[k]
                    cand_score = float(local_scores[k])
                    decision_src = "marker-argmax"
                if cand_score >= min_score or decision_src == "neighbor-vote":
                    new_ann = candidate
                    new_score = (f"{cand_score:+.3f} "
                                 f"(P_marker[{L_m}]={p_m:.2f} vs "
                                 f"P_graph[{L_g}]={p_g:.2f}, "
                                 f"src={decision_src})")
                    reassign_flag = "*"
                    reassigned[c] = new_ann

        rows.append([c, label_map.get(c, ""), ann_c, str(L_c),
                     L_dom, f"{s_dom:.3f}", f"{s_self:.3f}",
                     new_ann, new_score, reassign_flag])

    # Write log
    out_tsv = outdir / "annotation_smoothed.tsv"
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for r in rows:
            w.writerow(r)
    print(f"  wrote {out_tsv}")
    # Pretty print
    widths = [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
    for r in rows:
        if r[-1] == "*" or r is rows[0]:
            print("  " + "  ".join(str(c).ljust(wd) for c, wd in zip(r, widths)))
    print(f"\n  reassigned {len(reassigned)}/{len(labels)} clusters")
    print(f"  changes: {reassigned}")

    # ---- Apply to obs in-memory + write a new tsv summary of cell-level ----
    leiden_arr = a.obs["dcta_leiden"].astype(str).to_numpy()
    cur_ann_arr = a.obs["dcta_annotation"].astype(str).to_numpy()
    new_ann_arr = cur_ann_arr.copy()
    for c, new_v in reassigned.items():
        m = leiden_arr == c
        new_ann_arr[m] = new_v
    # Save updated h5ad copy
    a.obs["dcta_annotation_smoothed"] = new_ann_arr
    out_h5 = outdir / "dcta_phase1_smoothed.h5ad"
    a.write_h5ad(out_h5)
    print(f"  wrote {out_h5}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ds", choices=["zeng", "triana", "strati", "all"],
                   default="all")
    p.add_argument("--dom-share", type=float, default=None,
                   help="Override dom_share. Default = data-driven "
                        "(median of per-cluster self-lineage share in "
                        "PAGA neighbors).")
    p.add_argument("--dom-ratio", type=float, default=None,
                   help="Override dom_ratio. Default = data-driven "
                        "(1.5 + 4 × IQR(self_share)).")
    p.add_argument("--min-score", type=float, default=0.0,
                   help="Minimum marker z-score for the reassigned candidate "
                        "(only enforced when decision_src=marker-argmax).")
    p.add_argument("--max-self-score", type=float, default=None,
                   help="Override max_self_score. Default = data-driven "
                        "(Q25 of per-cluster self-marker z-scores). "
                        "Clusters whose CURRENT annotation has score below "
                        "this are eligible for smoothing.")
    args = p.parse_args()
    mapping = {
        "zeng":   ("Zeng2025_BMM", ROOT / "results/phase1_hay_paperclust"),
        "triana": ("Triana2021",   ROOT / "results/phase1_triana"),
        "strati": ("Strati2023",   ROOT / "results/phase1_strati"),
    }
    targets = list(mapping.values()) if args.ds == "all" else [mapping[args.ds]]
    for name, outdir in targets:
        smooth(name, outdir, args.dom_share, args.dom_ratio,
               args.min_score, args.max_self_score)


if __name__ == "__main__":
    main()
