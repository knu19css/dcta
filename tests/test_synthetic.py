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
    out = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=False)
    for col in ["n_counts", "pct_mito"]:
        assert col in out.obs.columns
    assert "highly_variable" in out.var.columns
    assert out.var["highly_variable"].sum() > 0
    assert out.raw is not None


def test_preprocessing_cell_cycle_columns(synthetic_adata):
    out = dcta.preprocessing.run(synthetic_adata, regress_cell_cycle=True)
    assert "S_score" in out.obs.columns
    assert "G2M_score" in out.obs.columns


def test_preprocessing_rejects_negative_X(synthetic_adata):
    a = synthetic_adata.copy()
    a.X = a.X.astype(float)
    a.X[0, 0] = -1.0
    with pytest.raises(dcta.DCTADataError):
        dcta.preprocessing.run(a)
