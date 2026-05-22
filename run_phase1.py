"""DCTA Phase 1 CLI orchestrator.

Runs the full pipeline:
  load -> preprocessing -> diffusion -> connectivity -> pseudotime -> dag -> viz

Usage:
  python run_phase1.py --h5ad /path/to/hay_hca_balanced_bm.h5ad \
                       --outdir results/phase1_hay
"""
from __future__ import annotations
import argparse
import os
import sys
import json
import networkx as nx

import dcta
from dcta import defaults as D


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="DCTA Phase 1 pipeline")
    p.add_argument("--h5ad", required=True, help="input .h5ad file")
    p.add_argument("--outdir", default="results/phase1_hay", help="output directory")
    p.add_argument("--hsc-yaml", default="data/priors/hsc_markers.yaml",
                   help="YAML file with HSC marker genes")
    p.add_argument("--use-rep", default="X_harmony",
                   help="adata.obsm key for the publisher-provided embedding "
                        "to build KNN on (e.g. X_harmony, X_pca, X_umap). "
                        "DCTA no longer recomputes PCA or diffusion map; the "
                        "embedding must already be attached.")
    p.add_argument("--knn-diffmap-k", type=int, default=D.KNN_DIFFMAP_K)
    p.add_argument("--resolution", type=float, default=D.LEIDEN_RESOLUTION)
    p.add_argument("--edge-threshold", type=float, default=D.EDGE_THRESHOLD)
    p.add_argument("--robustness-k", type=int, default=D.ROBUSTNESS_K)
    p.add_argument("--transitive-reduction", action="store_true", default=False)
    p.add_argument("--seniority-yaml", default="data/priors/lineage_seniority.yaml",
                   help="YAML mapping cell-type labels to seniority rank (for direction-accuracy)")
    p.add_argument("--label-column", default="cell_type",
                   help="adata.obs column used to look up seniority")
    p.add_argument("--direction-strategy", default="auto_select",
                   choices=["auto_select", "marker", "dpt"])
    p.add_argument("--dag-prior-yaml",
                   default="data/priors/dag_prior_bm.yaml,data/priors/dag_prior_t_maturation.yaml",
                   help="Optional prior hematopoietic DAG (comma-separated "
                        "list of yamls; merged as union of nodes/edges). "
                        "When provided together with --lineage-marker-yaml, "
                        "edges of G_data are scored for prior consistency "
                        "and backbone penalizes anti-prior edges. Pass '' "
                        "to disable.")
    p.add_argument("--prior-inconsistent-penalty", type=float, default=0.5,
                   help="Multiplicative penalty (0..1) applied to edge weight "
                        "during backbone calculation when edge is prior-"
                        "inconsistent. 0=no penalty, 1=drop entirely.")
    p.add_argument("--prior-fusion-weight", type=float, default=0.0,
                   help="Phase 2.6: per-cluster label fusion weight alpha (0..1). "
                        "0 = marker argmax only (no DAG influence on label). "
                        ">0 = each cluster's label score blends its own marker "
                        "score with the mean marker scores of G_prior neighbors "
                        "evaluated on PAGA-adjacent clusters. Use to demonstrate "
                        "that DAG topology changes cell-type calls.")
    p.add_argument("--prior-fusion-hops", type=int, default=1,
                   help="K-hop neighborhood in G_prior used for fusion (default 1).")
    p.add_argument("--ambiguity-margin", type=float, default=0.20,
                   help="Phase 2.6 disambiguation. When the cluster's top-1 "
                        "marker score - top-2 score is below this fraction of "
                        "top-1, treat as ambiguous and apply DAG-aware rules: "
                        "(A) if top-2 share an LCA progenitor in G_prior, "
                        "use the LCA label; (B) if top-2 are both terminals, "
                        "tie-break by PAGA-neighbor marker support of each "
                        "candidate's G_prior neighborhood.")
    p.add_argument("--lineage-marker-yaml", default="data/priors/lineage_mature_markers.yaml",
                   help="Optional cell-type marker panel YAML; when provided, "
                        "marker direction uses (mature_panel - hsc_score) instead "
                        "of plain HSC immaturity. Pass '' to disable.")
    p.add_argument("--paper-embedding-h5ad", default=None,
                   help="Optional source h5ad (e.g. cellxgene download) "
                        "to graft paper-original X_umap/X_pca onto adata.obsm. "
                        "Matched on observation_joinid. Also pulls "
                        "AuthorCellType/AuthorCellType_Broad/Study into obs "
                        "as paper_<col>.")
    p.add_argument("--cluster-from-obs", default=None,
                   help="If set, SKIP Leiden and use this obs column as "
                        "the cluster assignment for PAGA + downstream "
                        "(e.g. 'paper_AuthorCellType_Broad'). Requires "
                        "the column to exist in adata.obs after preprocessing.")
    p.add_argument("--keep-cluster-labels", action="store_true", default=False,
                   help="Preserve the original cluster labels from "
                        "--cluster-from-obs (default anonymizes them to "
                        "C00, C01, ... to prevent label leakage when the "
                        "source column is a paper-curated cell-type label).")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    print(f"[DCTA] loading {args.h5ad}", flush=True)
    adata = dcta.io.load_h5ad(args.h5ad)
    print(f"[DCTA] adata: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    if args.paper_embedding_h5ad:
        print(f"[DCTA] attaching paper embedding from {args.paper_embedding_h5ad}", flush=True)
        adata = dcta.io.attach_paper_embedding(adata, args.paper_embedding_h5ad)
        meta = adata.uns["dcta"]["paper_embedding_source"]
        print(f"[DCTA] paper embedding: matched {meta['n_matched']}/{meta['n_target']} "
              f"cells; attached {meta['embeddings_attached']}", flush=True)

    # Phase 2.8: skip local preprocessing + diffusion. Use the publisher's
    # log-normalized X and obsm[use_rep] (X_harmony / X_pca / X_umap) as-is.
    # If the input dataset comes with raw counts, the user must normalize
    # upstream — DCTA no longer rebuilds embeddings.
    print(f"[DCTA] using paper embedding: {args.use_rep}", flush=True)
    if args.use_rep not in adata.obsm:
        raise SystemExit(
            f"--use-rep={args.use_rep!r} not in adata.obsm "
            f"(available: {list(adata.obsm.keys())}). "
            f"Pass --paper-embedding-h5ad or change --use-rep."
        )

    print("[DCTA] connectivity", flush=True)
    cluster_from_obs = args.cluster_from_obs or None
    if cluster_from_obs and cluster_from_obs not in adata.obs.columns:
        raise SystemExit(
            f"--cluster-from-obs={cluster_from_obs!r} not found in adata.obs "
            f"(available: {list(adata.obs.columns)[:20]} ...). "
            f"Did you pass --paper-embedding-h5ad to pull paper_* columns?"
        )
    adata = dcta.connectivity.run(adata, knn_k=args.knn_diffmap_k,
                                  resolution=args.resolution,
                                  cluster_from_obs=cluster_from_obs,
                                  use_rep=args.use_rep,
                                  anonymize_cluster_labels=not args.keep_cluster_labels)
    n_clust = adata.obs["dcta_leiden"].nunique()
    src = adata.uns["dcta"].get("cluster_source", "?")
    print(f"[DCTA] clusters: {n_clust} (source: {src})", flush=True)

    print("[DCTA] pseudotime", flush=True)
    lineage_yaml = args.lineage_marker_yaml or None
    adata = dcta.pseudotime.run(adata, hsc_marker_yaml=args.hsc_yaml,
                                lineage_marker_yaml=lineage_yaml,
                                robustness_k=args.robustness_k)
    print(f"[DCTA] root cluster: {adata.uns['dcta']['root_cluster']}, "
          f"components: {adata.uns['dcta']['component_count']}", flush=True)

    print("[DCTA] dag", flush=True)
    raw = (args.dag_prior_yaml or "").strip()
    if not raw:
        dag_prior_yaml = None
    elif "," in raw:
        dag_prior_yaml = [p.strip() for p in raw.split(",") if p.strip()]
    else:
        dag_prior_yaml = raw
    adata = dcta.dag.build(adata, edge_threshold=args.edge_threshold,
                           transitive_reduction=args.transitive_reduction,
                           direction_strategy=args.direction_strategy,
                           seniority_yaml=args.seniority_yaml,
                           label_column=args.label_column,
                           dag_prior_yaml=dag_prior_yaml,
                           lineage_marker_yaml=lineage_yaml,
                           prior_inconsistent_penalty=args.prior_inconsistent_penalty,
                           prior_fusion_weight=args.prior_fusion_weight,
                           prior_fusion_hops=args.prior_fusion_hops,
                           ambiguity_margin=args.ambiguity_margin)
    G = adata.uns["dcta"]["G_data"]
    chosen = adata.uns["dcta"]["direction_method_chosen"]
    acc = adata.uns["dcta"]["direction_methods"].get("accuracy", {})
    print(f"[DCTA] G_data: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"{len(adata.uns['dcta']['cycle_removal_log'])} edges removed, "
          f"direction chosen={chosen}, accuracy={acc}", flush=True)

    # Phase 2.3: cluster-level annotation via absorbing chain on G_prior
    if dag_prior_yaml is not None:
        print("[DCTA] annotation (absorbing markov on G_prior)", flush=True)
        adata = dcta.markov.annotate_clusters(adata, dag_prior_yaml=dag_prior_yaml)
        annot_counts = adata.obs["dcta_annotation"].astype(object).value_counts(dropna=False)
        print(f"[DCTA] annotation distribution: "
              f"{dict(annot_counts.head(10))} ...", flush=True)

        if args.label_column and args.label_column in adata.obs.columns:
            eval_out = dcta.annotation_eval.evaluate_annotation(
                adata, gt_column=args.label_column,
                pred_column="dcta_annotation",
                lineage_marker_yaml=lineage_yaml,
            )
            print(f"[DCTA] annotation eval (DCTA full): "
                  f"strict_acc={eval_out['strict_accuracy']:.3f} "
                  f"lineage_acc={eval_out['lineage_accuracy']:.3f} "
                  f"(n_evaluable={eval_out['n_evaluable']})", flush=True)
            adata.uns["dcta"]["annotation_eval"] = eval_out

            # ---- Marker-only baseline (control): broadcast cluster_marker_only_label ----
            mol = adata.uns["dcta"].get("cluster_marker_only_label", {})
            if mol:
                leiden_str = adata.obs["dcta_leiden"].astype(str).to_numpy()
                adata.obs["dcta_marker_only_annotation"] = [
                    mol.get(c) for c in leiden_str
                ]
                eval_mo = dcta.annotation_eval.evaluate_annotation(
                    adata, gt_column=args.label_column,
                    pred_column="dcta_marker_only_annotation",
                    lineage_marker_yaml=lineage_yaml,
                )
                print(f"[DCTA] annotation eval (marker-only baseline): "
                      f"strict_acc={eval_mo['strict_accuracy']:.3f} "
                      f"lineage_acc={eval_mo['lineage_accuracy']:.3f}",
                      flush=True)
                adata.uns["dcta"]["annotation_eval_marker_only"] = eval_mo
                cl_annot = adata.uns["dcta"]["cluster_annotation"]
                diffs = [
                    (c, mol.get(c), cl_annot.get(c))
                    for c in mol if mol.get(c) != cl_annot.get(c)
                ]
                print(f"[DCTA] clusters changed by DCTA vs marker-only: "
                      f"{len(diffs)}/{len(mol)}", flush=True)
                for c, b, d in diffs[:10]:
                    print(f"         {c} : {b} -> {d}", flush=True)
                adata.uns["dcta"]["cluster_label_diff_vs_marker_only"] = [
                    {"cluster": c, "marker_only": b, "dcta": d}
                    for c, b, d in diffs
                ]

    out_h5ad = os.path.join(args.outdir, "dcta_phase1.h5ad")
    out_gml = os.path.join(args.outdir, "G_data.gml")
    # Pop non-h5ad-serializable items before writing; restore after so
    # downstream viz still works on adata.
    G_obj = adata.uns["dcta"].pop("G_data", None)
    bb_obj = adata.uns["dcta"].pop("backbone_validation", None)
    dm_obj = adata.uns["dcta"].pop("direction_methods", None)
    pe_obj = adata.uns["dcta"].pop("paper_embedding_source", None)
    cls_obj = adata.uns["dcta"].pop("cluster_lineage_scores", None)
    cbl_obj = adata.uns["dcta"].pop("cluster_best_lineage", None)
    pc_obj = adata.uns["dcta"].pop("prior_consistency", None)
    cpa_obj = adata.uns["dcta"].pop("cluster_prior_alternatives", None)
    am_obj = adata.uns["dcta"].pop("absorption_matrix", None)
    ae_obj_full = adata.uns["dcta"].pop("annotation_eval", None)
    ae_obj = ae_obj_full  # alias for summary
    mol_obj = adata.uns["dcta"].pop("cluster_marker_only_label", None)
    moa_obj = adata.uns["dcta"].pop("cluster_marker_only_alternatives", None)
    emo_obj = adata.uns["dcta"].pop("annotation_eval_marker_only", None)
    cld_obj = adata.uns["dcta"].pop("cluster_label_diff_vs_marker_only", None)
    dr_obj = adata.uns["dcta"].pop("cluster_disambig_reason", None)
    # absorption_matrix has nested lists+strs, JSON-safe but huge; keep popped
    print(f"[DCTA] writing {out_h5ad}", flush=True)
    adata.write_h5ad(out_h5ad)
    if G_obj is not None:
        adata.uns["dcta"]["G_data"] = G_obj
    if bb_obj is not None:
        adata.uns["dcta"]["backbone_validation"] = bb_obj
    if dm_obj is not None:
        adata.uns["dcta"]["direction_methods"] = dm_obj
    if pe_obj is not None:
        adata.uns["dcta"]["paper_embedding_source"] = pe_obj
    if cls_obj is not None:
        adata.uns["dcta"]["cluster_lineage_scores"] = cls_obj
    if cbl_obj is not None:
        adata.uns["dcta"]["cluster_best_lineage"] = cbl_obj
    if pc_obj is not None:
        adata.uns["dcta"]["prior_consistency"] = pc_obj
    if cpa_obj is not None:
        adata.uns["dcta"]["cluster_prior_alternatives"] = cpa_obj
    if am_obj is not None:
        adata.uns["dcta"]["absorption_matrix"] = am_obj
    if ae_obj is not None:
        adata.uns["dcta"]["annotation_eval"] = ae_obj
    if mol_obj is not None:
        adata.uns["dcta"]["cluster_marker_only_label"] = mol_obj
    if moa_obj is not None:
        adata.uns["dcta"]["cluster_marker_only_alternatives"] = moa_obj
    if emo_obj is not None:
        adata.uns["dcta"]["annotation_eval_marker_only"] = emo_obj
    if cld_obj is not None:
        adata.uns["dcta"]["cluster_label_diff_vs_marker_only"] = cld_obj
    if dr_obj is not None:
        adata.uns["dcta"]["cluster_disambig_reason"] = dr_obj
    print(f"[DCTA] writing {out_gml}", flush=True)
    # nx.write_gml requires str-keyed attrs and no None values
    G_clean = nx.DiGraph()
    for n_, data in G.nodes(data=True):
        G_clean.add_node(str(n_), **{k: float(v) if v is not None else 0.0
                                     for k, v in data.items()})
    for u, v, data in G.edges(data=True):
        G_clean.add_edge(str(u), str(v), **{k: float(val) for k, val in data.items()})
    nx.write_gml(G_clean, out_gml)

    print("[DCTA] exporting figures", flush=True)
    dcta.viz.export_all(adata, outdir=args.outdir)

    methods = adata.uns["dcta"].get("direction_methods", {})
    pc = pc_obj or adata.uns["dcta"].get("prior_consistency", {})
    summary = {
        "n_cells_after_qc": int(adata.n_obs),
        "n_genes_after_qc": int(adata.n_vars),
        "n_rep_dims": int(adata.uns["dcta"].get("n_rep_dims", 0)),
        "use_rep": adata.uns["dcta"].get("use_rep", ""),
        "n_clusters": int(n_clust),
        "n_components": int(adata.uns["dcta"]["component_count"]),
        "root_cluster": adata.uns["dcta"]["root_cluster"],
        "n_edges_initial": int(G.number_of_edges() + len(adata.uns["dcta"]["cycle_removal_log"])),
        "n_edges_after_fas": int(G.number_of_edges()),
        "n_edges_removed": int(len(adata.uns["dcta"]["cycle_removal_log"])),
        "direction": {
            "method_chosen": chosen,
            "accuracy": acc,
            "agreement_rate": methods.get("agreement_rate"),
            "n_pairs_compared": methods.get("n_pairs_compared"),
        },
        "backbone": {
            k: v for k, v in adata.uns["dcta"].get("backbone_validation", {}).items()
            if k not in ("backbone_seniority",)  # keep summary JSON-clean
        },
        "prior": {
            k: pc.get(k) for k in (
                "consistency_rate", "consistent_count", "scored_count",
                "unlabeled_count", "penalty", "g_prior_path",
            )
        } if pc else None,
        "annotation": ({
            "strict_accuracy": ae_obj["strict_accuracy"],
            "lineage_accuracy": ae_obj["lineage_accuracy"],
            "n_evaluable": ae_obj["n_evaluable"],
            "n_unmapped_gt_labels": len(ae_obj.get("unmapped_gt_labels", [])),
            "n_unmapped_pred_labels": len(ae_obj.get("unmapped_pred_labels", [])),
        } if ae_obj else None),
        "marker_only_baseline": ({
            "strict_accuracy": emo_obj["strict_accuracy"],
            "lineage_accuracy": emo_obj["lineage_accuracy"],
            "n_clusters_changed_by_dcta": len(cld_obj) if cld_obj else 0,
            "fusion_weight": args.prior_fusion_weight,
            "fusion_hops": args.prior_fusion_hops,
        } if emo_obj else None),
    }
    # backbone_seniority contains tuples — JSON-serialize as list of [c, s]
    bb_full = adata.uns["dcta"].get("backbone_validation", {})
    if bb_full.get("backbone_seniority"):
        summary["backbone"]["backbone_seniority"] = [
            [c, s] for (c, s) in bb_full["backbone_seniority"]
        ]
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[DCTA] DONE. summary: {summary}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
