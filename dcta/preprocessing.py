"""DCTA preprocessing: QC, normalization, HVG, optional cell-cycle regress-out."""
from __future__ import annotations
import numpy as np
import scanpy as sc
import anndata as ad
from dcta import DCTADataError
from dcta import defaults as D


# Tirosh et al. 2016 cell cycle gene lists (subset, human symbols).
# Source: regev_lab_cell_cycle_genes.txt (commonly distributed with scanpy).
_S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
    "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
]
_G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


def _validate(adata: ad.AnnData) -> None:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise DCTADataError("AnnData is empty (0 cells or 0 genes)")
    X = adata.X
    Xmin = X.min() if not hasattr(X, "toarray") else X.min()
    if Xmin < 0:
        raise DCTADataError(f"AnnData.X contains negative values (min={Xmin})")


def run(
    adata: ad.AnnData,
    *,
    copy: bool = False,
    regress_cell_cycle: bool = True,
    min_genes: int = D.MIN_GENES_PER_CELL,
    min_cells: int = D.MIN_CELLS_PER_GENE,
    pct_mito_max: float = D.PCT_MITO_MAX,
    n_top_genes: int = D.N_TOP_GENES,
) -> ad.AnnData:
    """Preprocess raw count AnnData for DCTA.

    Steps:
      1. Validate (no empty matrix, no negative values).
      2. QC filter (min_genes, min_cells, pct_mito).
      3. Save raw counts to .raw.
      4. Normalize total + log1p.
      5. Highly-variable genes.
      6. Optional cell-cycle scoring + regress-out.

    Returns the same AnnData (mutated) when copy=False.
    """
    a = adata.copy() if copy else adata
    _validate(a)

    # Mitochondrial percentage
    a.var["mt"] = a.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    a.obs["n_counts"] = a.obs["total_counts"]
    a.obs["pct_mito"] = a.obs["pct_counts_mt"]

    # QC filter
    if min_genes > 0:
        sc.pp.filter_cells(a, min_genes=min_genes)
    if a.n_obs == 0:
        raise DCTADataError(
            f"min_genes={min_genes} filtered out every cell. "
            f"Pass a smaller min_genes in tests or check input data."
        )
    sc.pp.filter_genes(a, min_cells=min_cells)
    a = a[a.obs["pct_mito"] <= pct_mito_max, :].copy()

    # Save raw counts
    a.raw = a

    # Normalize + log
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # HVG
    sc.pp.highly_variable_genes(a, n_top_genes=n_top_genes, flavor="seurat")

    # Cell cycle
    s_present = [g for g in _S_GENES if g in a.var_names]
    g2m_present = [g for g in _G2M_GENES if g in a.var_names]
    if s_present and g2m_present:
        sc.tl.score_genes_cell_cycle(a, s_genes=s_present, g2m_genes=g2m_present)
    else:
        # No cell-cycle genes available — write NaN columns so downstream code
        # can rely on the columns existing
        a.obs["S_score"] = np.nan
        a.obs["G2M_score"] = np.nan
        a.obs["phase"] = "unknown"

    if regress_cell_cycle and s_present and g2m_present:
        sc.pp.regress_out(a, ["S_score", "G2M_score"])
        sc.pp.scale(a, max_value=10)

    # Mirror result into original adata if copy=False and slicing happened
    return a
