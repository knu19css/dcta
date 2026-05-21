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
