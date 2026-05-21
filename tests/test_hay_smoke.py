"""End-to-end smoke test on a 1000-cell stratified subset of the Hay BM dataset.

Marked slow — runs only when explicitly selected.
"""
import os
import pytest
import numpy as np
import anndata as ad
import networkx as nx
from pathlib import Path

import dcta

HAY_PATH = "/Users/edemhw/Honglab/CS/04_Claude/10_bm_autoannotation/test_datasets/cxg/hay_hca_balanced_bm.h5ad"
HSC_YAML = str(Path(__file__).resolve().parent.parent / "data" / "priors" / "hsc_markers.yaml")
SENIORITY_YAML = str(Path(__file__).resolve().parent.parent / "data" / "priors" / "lineage_seniority.yaml")


pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def hay_subset():
    if not os.path.exists(HAY_PATH):
        pytest.skip(f"Hay dataset not found at {HAY_PATH}")
    a = ad.read_h5ad(HAY_PATH)
    # CXG datasets often store integer var_names with gene symbols in 'feature_name'.
    # Remap to symbols (dedup) so HSC marker lookup works.
    if "feature_name" in a.var.columns and not any(
        g in a.var_names for g in ("CRHBP", "AVP", "HLF")
    ):
        symbols = a.var["feature_name"].astype(str).to_numpy()
        # dedup: keep first occurrence
        seen = {}
        keep = np.zeros(len(symbols), dtype=bool)
        for i, s in enumerate(symbols):
            if s not in seen:
                seen[s] = i
                keep[i] = True
        a = a[:, keep].copy()
        a.var_names = a.var["feature_name"].astype(str).to_numpy()
        a.var_names_make_unique()
    rng = np.random.default_rng(0)
    target = 1000
    if "cell_type" in a.obs.columns:
        groups = a.obs.groupby("cell_type", observed=True).indices
        per_group = max(1, target // max(1, len(groups)))
        picks = []
        for _, idx in groups.items():
            n = min(per_group, len(idx))
            picks.extend(rng.choice(idx, size=n, replace=False).tolist())
        picks = picks[:target]
    else:
        picks = rng.choice(a.n_obs, size=min(target, a.n_obs), replace=False).tolist()
    sub = a[picks, :].copy()
    return sub


def test_hay_phase1_end_to_end(hay_subset, tmp_path):
    a = hay_subset
    a = dcta.preprocessing.run(a, regress_cell_cycle=True)
    a = dcta.diffusion.run(a, n_comps=15)
    a = dcta.connectivity.run(a, resolution=1.0)
    a = dcta.pseudotime.run(a, hsc_marker_yaml=HSC_YAML, robustness_k=3)
    a = dcta.dag.build(a, edge_threshold=0.05, seniority_yaml=SENIORITY_YAML,
                       label_column="cell_type")

    G = a.uns["dcta"]["G_data"]
    assert isinstance(G, nx.DiGraph)
    assert nx.is_directed_acyclic_graph(G)
    assert G.number_of_nodes() >= 3, f"only {G.number_of_nodes()} clusters survived QC"

    root = a.uns["dcta"]["root_cluster"]
    assert root in G.nodes
    leaves = [n for n in G.nodes if G.out_degree(n) == 0]
    reach = [n for n in leaves if nx.has_path(G, root, n)]
    assert len(reach) >= 1, "root cluster cannot reach any leaf"

    dcta.viz.export_all(a, outdir=str(tmp_path))
    for f in ["diffmap_eigvals.png", "paga_connectivities.png",
              "dpt_distribution.png", "dpt_root_std_hist.png",
              "G_data_diagram.png", "paga_vs_dag_overlay.png"]:
        assert (tmp_path / f).exists()
    assert a.uns["dcta"]["direction_method_chosen"] in ("marker", "dpt")
