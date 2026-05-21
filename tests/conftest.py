"""Synthetic AnnData fixtures for DCTA tests."""
import numpy as np
import pandas as pd
import anndata as ad
import pytest


@pytest.fixture
def synthetic_adata():
    """3-lineage star graph: HSC -> {Mye, Ery, Lym}, 500 cells, 100 genes.

    Cells per group:
      HSC: 125, Mye: 125, Ery: 125, Lym: 125
    Genes:
      - 8 HSC markers (CRHBP, AVP, HLF, MEIS1, MLLT3, HOPX, HEMGN, PROM1)
      - 10 per-lineage markers (MYE_1..10, ERY_1..10, LYM_1..10)
      - 62 noise genes (NOISE_1..62)
    Counts are negative-binomial-like ints. Lineage-marker expression is
    enriched in the corresponding lineage; HSC markers are enriched in HSC.
    Returns a raw-count AnnData (no normalization done).
    """
    rng = np.random.default_rng(seed=0)
    n_per_group = 125
    n_obs = 4 * n_per_group  # 500
    lineages = (["HSC"] * n_per_group + ["Mye"] * n_per_group +
                ["Ery"] * n_per_group + ["Lym"] * n_per_group)

    hsc_markers = ["CRHBP", "AVP", "HLF", "MEIS1", "MLLT3", "HOPX", "HEMGN", "PROM1"]
    mye_markers = [f"MYE_{i}" for i in range(1, 11)]
    ery_markers = [f"ERY_{i}" for i in range(1, 11)]
    lym_markers = [f"LYM_{i}" for i in range(1, 11)]
    noise_genes = [f"NOISE_{i}" for i in range(1, 63)]
    var_names = hsc_markers + mye_markers + ery_markers + lym_markers + noise_genes
    n_vars = len(var_names)

    # Baseline low expression for all
    base = rng.negative_binomial(n=2, p=0.4, size=(n_obs, n_vars)).astype(np.float32)

    # Enrichment helper
    def enrich(group_label, gene_list, scale=20):
        gidx = [var_names.index(g) for g in gene_list]
        cmask = np.array([l == group_label for l in lineages])
        base[np.ix_(cmask, gidx)] += rng.negative_binomial(n=8, p=0.3,
                                                           size=(cmask.sum(), len(gidx))) * (scale / 10)

    enrich("HSC", hsc_markers, scale=25)
    enrich("Mye", mye_markers, scale=25)
    enrich("Ery", ery_markers, scale=25)
    enrich("Lym", lym_markers, scale=25)

    # Add a couple of mitochondrial-prefix genes (small, low expression)
    mt_names = ["MT-CO1", "MT-CO2", "MT-ND1"]
    mt_counts = rng.negative_binomial(n=1, p=0.5, size=(n_obs, len(mt_names))).astype(np.float32)
    X = np.concatenate([base, mt_counts], axis=1)
    var_names = var_names + mt_names

    obs = pd.DataFrame({"true_lineage": lineages,
                        "cell_id": [f"cell_{i:04d}" for i in range(n_obs)]})
    obs = obs.set_index("cell_id")
    var = pd.DataFrame(index=pd.Index(var_names, name="gene"))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


@pytest.fixture
def synthetic_adata_disconnected(synthetic_adata):
    """Same as synthetic_adata but with two cells groups artificially split.

    Used to test multi-root / multi-component handling in pseudotime.
    Splits cells into 2 disconnected components by zeroing out shared genes
    between the two halves.
    """
    a = synthetic_adata.copy()
    half = a.n_obs // 2
    # Shared HSC markers — zero them in second half to disconnect manifold
    hsc_idx = [a.var_names.get_loc(g) for g in ["CRHBP", "AVP", "HLF"]]
    a.X[half:, hsc_idx] = 0
    return a
