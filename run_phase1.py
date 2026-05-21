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
    p.add_argument("--n-comps-pca", type=int, default=D.N_COMPS_PCA)
    p.add_argument("--n-comps-diffmap", type=int, default=D.N_COMPS_DIFFMAP)
    p.add_argument("--knn-diffmap-k", type=int, default=D.KNN_DIFFMAP_K)
    p.add_argument("--resolution", type=float, default=D.LEIDEN_RESOLUTION)
    p.add_argument("--edge-threshold", type=float, default=D.EDGE_THRESHOLD)
    p.add_argument("--robustness-k", type=int, default=D.ROBUSTNESS_K)
    p.add_argument("--regress-cell-cycle", action="store_true", default=True)
    p.add_argument("--no-regress-cell-cycle", action="store_false",
                   dest="regress_cell_cycle")
    p.add_argument("--transitive-reduction", action="store_true", default=False)
    p.add_argument("--seniority-yaml", default="data/priors/lineage_seniority.yaml",
                   help="YAML mapping cell-type labels to seniority rank (for direction-accuracy)")
    p.add_argument("--label-column", default="cell_type",
                   help="adata.obs column used to look up seniority")
    p.add_argument("--direction-strategy", default="auto_select",
                   choices=["auto_select", "marker", "dpt"])
    p.add_argument("--dag-prior-yaml", default="data/priors/dag_prior_bm.yaml",
                   help="Optional prior hematopoietic DAG. When provided "
                        "together with --lineage-marker-yaml, edges of G_data "
                        "are scored for prior consistency and backbone "
                        "penalizes anti-prior edges. Pass '' to disable.")
    p.add_argument("--prior-inconsistent-penalty", type=float, default=0.5,
                   help="Multiplicative penalty (0..1) applied to edge weight "
                        "during backbone calculation when edge is prior-"
                        "inconsistent. 0=no penalty, 1=drop entirely.")
    p.add_argument("--lineage-marker-yaml", default="data/priors/lineage_mature_markers.yaml",
                   help="Optional cell-type marker panel YAML; when provided, "
                        "marker direction uses (mature_panel - hsc_score) instead "
                        "of plain HSC immaturity. Pass '' to disable.")
    p.add_argument("--paper-embedding-h5ad", default=None,
                   help="Optional source h5ad (e.g. cellxgene download) "
                        "to graft paper-original X_umap/X_pca onto adata.obsm. "
                        "Matched on observation_joinid.")
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

    print("[DCTA] preprocessing", flush=True)
    adata = dcta.preprocessing.run(adata, regress_cell_cycle=args.regress_cell_cycle)
    print(f"[DCTA] after QC: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    print("[DCTA] diffusion", flush=True)
    adata = dcta.diffusion.run(adata, n_comps_pca=args.n_comps_pca,
                               n_comps=args.n_comps_diffmap)
    print(f"[DCTA] n_dcs selected: {adata.uns['dcta']['n_dcs']}", flush=True)

    print("[DCTA] connectivity", flush=True)
    adata = dcta.connectivity.run(adata, knn_k=args.knn_diffmap_k,
                                  resolution=args.resolution)
    n_clust = adata.obs["dcta_leiden"].nunique()
    print(f"[DCTA] leiden clusters: {n_clust}", flush=True)

    print("[DCTA] pseudotime", flush=True)
    lineage_yaml = args.lineage_marker_yaml or None
    adata = dcta.pseudotime.run(adata, hsc_marker_yaml=args.hsc_yaml,
                                lineage_marker_yaml=lineage_yaml,
                                robustness_k=args.robustness_k)
    print(f"[DCTA] root cluster: {adata.uns['dcta']['root_cluster']}, "
          f"components: {adata.uns['dcta']['component_count']}", flush=True)

    print("[DCTA] dag", flush=True)
    dag_prior_yaml = args.dag_prior_yaml or None
    adata = dcta.dag.build(adata, edge_threshold=args.edge_threshold,
                           transitive_reduction=args.transitive_reduction,
                           direction_strategy=args.direction_strategy,
                           seniority_yaml=args.seniority_yaml,
                           label_column=args.label_column,
                           dag_prior_yaml=dag_prior_yaml,
                           lineage_marker_yaml=lineage_yaml,
                           prior_inconsistent_penalty=args.prior_inconsistent_penalty)
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
            print(f"[DCTA] annotation eval: "
                  f"strict_acc={eval_out['strict_accuracy']:.3f} "
                  f"lineage_acc={eval_out['lineage_accuracy']:.3f} "
                  f"(n_evaluable={eval_out['n_evaluable']})", flush=True)
            adata.uns["dcta"]["annotation_eval"] = eval_out

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
        "n_dcs": int(adata.uns["dcta"]["n_dcs"]),
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
