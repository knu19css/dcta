"""Per-stage unit tests on synthetic 3-lineage star AnnData."""
import pytest
import anndata as ad
import numpy as np


def test_synthetic_fixture_shape(synthetic_adata):
    assert isinstance(synthetic_adata, ad.AnnData)
    assert synthetic_adata.n_obs == 500
    assert synthetic_adata.n_vars >= 50
    # 4 ground-truth lineage labels: HSC, Mye, Ery, Lym
    assert set(synthetic_adata.obs["true_lineage"].unique()) == {"HSC", "Mye", "Ery", "Lym"}
    # HSC marker genes from priors must be in var_names
    for g in ["CRHBP", "AVP", "HLF"]:
        assert g in synthetic_adata.var_names


import os
import dcta


def test_io_load_h5ad_returns_anndata(tmp_path, synthetic_adata):
    p = tmp_path / "tiny.h5ad"
    synthetic_adata.write_h5ad(p)
    loaded = dcta.io.load_h5ad(str(p))
    assert isinstance(loaded, ad.AnnData)
    assert loaded.n_obs == synthetic_adata.n_obs
    assert loaded.n_vars == synthetic_adata.n_vars


def test_io_load_h5ad_missing_file_raises():
    with pytest.raises(dcta.DCTADataError):
        dcta.io.load_h5ad("/nonexistent/path/no.h5ad")


def test_preprocessing_run_adds_qc_obs(synthetic_adata):
    out = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False, min_genes=10, n_top_genes=50)
    for col in ["n_counts", "pct_mito"]:
        assert col in out.obs.columns
    assert "highly_variable" in out.var.columns
    assert out.var["highly_variable"].sum() > 0
    assert out.raw is not None


def test_preprocessing_cell_cycle_columns(synthetic_adata):
    out = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=True,
                                 min_genes=10, n_top_genes=50)
    assert "S_score" in out.obs.columns
    assert "G2M_score" in out.obs.columns


def test_preprocessing_rejects_negative_X(synthetic_adata):
    a = synthetic_adata.copy()
    a.X = a.X.astype(float)
    a.X[0, 0] = -1.0
    with pytest.raises(dcta.DCTADataError):
        dcta.preprocessing.run(a, min_genes=10, n_top_genes=50)


def test_preprocessing_raises_when_min_genes_filters_all(synthetic_adata):
    # synthetic has ~100 detectable genes per cell; min_genes=999 wipes all
    with pytest.raises(dcta.DCTADataError, match="filtered out every cell"):
        dcta.preprocessing.run(synthetic_adata, min_genes=999)


def test_diffusion_adds_pca_and_diffmap(synthetic_adata):
    a = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False, min_genes=10, n_top_genes=50)
    out = dcta.diffusion.run(a, n_comps=10, n_comps_pca=20)
    assert "X_pca" in out.obsm
    assert "X_diffmap" in out.obsm
    assert "dcta" in out.uns
    assert "diffmap_eigvals" in out.uns["dcta"]
    assert "n_dcs" in out.uns["dcta"]
    assert 1 <= out.uns["dcta"]["n_dcs"] <= 10
    eigs = np.asarray(out.uns["dcta"]["diffmap_eigvals"])
    assert (eigs[:-1] >= eigs[1:] - 1e-9).all()


def test_diffusion_n_dcs_auto_gap_selection(synthetic_adata):
    a = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False, min_genes=10, n_top_genes=50)
    out = dcta.diffusion.run(a, n_comps=15, n_comps_pca=20)
    assert out.uns["dcta"]["n_dcs"] < 15


def test_connectivity_runs_leiden_and_paga(synthetic_adata):
    a = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False, min_genes=10, n_top_genes=50)
    a = dcta.diffusion.run(a, n_comps=10, n_comps_pca=20)
    out = dcta.connectivity.run(a, resolution=1.0)
    assert "dcta_leiden" in out.obs.columns
    assert out.obs["dcta_leiden"].nunique() >= 2
    paga = out.uns["dcta"]["paga"]
    assert "connectivities" in paga
    assert "cluster_sizes" in paga
    conn = np.asarray(paga["connectivities"])
    k = out.obs["dcta_leiden"].nunique()
    assert conn.shape == (k, k)


def test_connectivity_resolution_affects_cluster_count(synthetic_adata):
    a = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False, min_genes=10, n_top_genes=50)
    a = dcta.diffusion.run(a, n_comps=10, n_comps_pca=20)
    low = dcta.connectivity.run(a.copy(), resolution=0.3)
    high = dcta.connectivity.run(a.copy(), resolution=2.0)
    assert low.obs["dcta_leiden"].nunique() <= high.obs["dcta_leiden"].nunique()


import yaml
from pathlib import Path

HSC_YAML = str(Path(__file__).resolve().parent.parent / "data" / "priors" / "hsc_markers.yaml")


def _prepped(adata):
    a = dcta.preprocessing.run(adata.copy(), regress_cell_cycle=False,
                               min_genes=10, n_top_genes=50)
    a = dcta.diffusion.run(a, n_comps=10, n_comps_pca=20)
    a = dcta.connectivity.run(a, resolution=1.0)
    return a


def test_pseudotime_run_adds_dpt(synthetic_adata):
    a = _prepped(synthetic_adata)
    out = dcta.pseudotime.run(a, hsc_marker_yaml=HSC_YAML)
    assert "dcta_dpt" in out.obs.columns
    assert "dcta_component" in out.obs.columns
    assert np.isfinite(out.obs["dcta_dpt"]).sum() > 0
    uns = out.uns["dcta"]
    for k in ["root_strategy", "root_cells", "root_cluster",
              "component_count", "dpt_root_std",
              "cluster_hsc_score", "cluster_maturation_score"]:
        assert k in uns, f"missing key: {k}"
    assert uns["component_count"] >= 1
    assert len(uns["root_cells"]) >= 1
    # cluster_hsc_score keys must match the leiden cluster labels
    leiden_set = set(out.obs["dcta_leiden"].astype(str).unique())
    assert set(uns["cluster_hsc_score"].keys()) == leiden_set


def test_pseudotime_raises_when_no_hsc_markers_present(synthetic_adata):
    a = _prepped(synthetic_adata)
    # Remove every HSC marker gene
    keep = [g for g in a.var_names if g not in
            ["CRHBP", "AVP", "HLF", "MEIS1", "MLLT3", "HOPX", "HEMGN", "PROM1", "CD34"]]
    a2 = a[:, keep].copy()
    with pytest.raises(dcta.DCTARootError):
        dcta.pseudotime.run(a2, hsc_marker_yaml=HSC_YAML)


def test_pseudotime_virtual_hsc_hook_not_implemented(synthetic_adata):
    a = _prepped(synthetic_adata)
    with pytest.raises(NotImplementedError):
        dcta.pseudotime.run(a, root_strategy="virtual_hsc_ref", ref_adata=a)


def test_pseudotime_multi_component(synthetic_adata_disconnected):
    a = _prepped(synthetic_adata_disconnected)
    out = dcta.pseudotime.run(a, hsc_marker_yaml=HSC_YAML)
    # Disconnected manifold should produce at least 1 component
    assert out.uns["dcta"]["component_count"] >= 1
