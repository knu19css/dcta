"""Prepare external BM datasets (Triana 2021, Strati 2023) for DCTA Phase-1.

Operations:
  - Read CXG-format h5ad (Ensembl IDs in var.index)
  - Swap var.index to feature_name (gene symbol); drop duplicates by averaging
  - For Strati: subset to disease_setting == 'Normal' (healthy donors)
  - Ensure 'paper_AuthorCellType_Broad' obs column exists
  - Ensure X_harmony or X_pca obsm key exists for DCTA's `use_rep`

This is a one-shot data-prep utility OUTSIDE the dcta package. It does not
modify any DCTA algorithm code.

Usage:
  python scripts/prep_external_bm.py --which triana
  python scripts/prep_external_bm.py --which strati
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external"
OUT = ROOT / "data" / "prepared"
OUT.mkdir(parents=True, exist_ok=True)


# Mapping from Triana 2021 Cell_label values → canonical DAG node names.
# Best-effort: ambiguous / rare populations left as None (treated as
# unmapped by DCTA eval).
TRIANA_LABEL_MAP = {
    "HSCs & MPPs": "HSC",
    "Lymphomyeloid prog": "MPP",
    "Erythro-myeloid progenitors": "CMP",
    "Megakaryocyte progenitors": "Megakaryocyte",
    "Eosinophil-basophil-mast cell progenitors": "BMCP",
    "Early erythroid progenitor": "BFU_E",
    "Late erythroid progenitor": "CFU_E",
    "Aberrant erythroid": "Erythroblast",
    "Pre-pro-B cells": "CLP",
    "Pro-B cells": "Pro_B_cell",
    "Pre-B cells": "Pre_B_cell",
    "Small pre-B cell": "Pre_B_cell",
    "Immature B cells": "B_cell",
    "Mature naive B cells": "B_cell",
    "Nonswitched memory B cells": "B_cell",
    "Class switched memory B cells": "B_cell",
    "CD11c+ memory B cells": "B_cell",
    "Plasma cells": "Plasma_cell",
    "GammaDelta T cells": "T_gd",
    "NK T cells": "NKT_cell",
    "NK cell progenitors": "CLP",
    "CD56dimCD16+ NK cells": "NK_cell",
    "CD56brightCD16- NK cells": "NK_cell",
    "CD4+ naive T cells": "Naive_CD4_T_cell",
    "CD8+ naive T cells": "Naive_CD8_T_cell",
    "CD8+ central memory T cells": "CD8_central_memory",
    "CD8+ effector memory T cells": "CD8_effector_memory",
    "CD8+CD103+ tissue resident memory T cells": "CD8_effector_memory",
    "Classical Monocytes": "Monocyte",
    "Non-classical monocytes": "Monocyte",
    "Early promyelocytes": "GMP",
    "Late promyelocytes": "GMP",
    "Myelocytes": "Neutrophil",
    "Conventional dendritic cell 1": "cDC1",
    "Conventional dendritic cell 2": "cDC2",
    "Plasmacytoid dendritic cells": "pDC",
    "Plasmacytoid dendritic cell progenitors": "pDC",
    # Unmapped (best left blank for eval):
    # "CD4+ memory T cells" — generic memory CD4 not in DAG
    # "CD4+ cytotoxic T cells" — rare subset
    # "CD69+PD-1+ memory CD4+ T cells" — exhaustion subset
    # "Mesenchymal cells_1" — not hematopoietic
}


# Mapping from Strati 2023 cluster column ("C00:CD8 T", ...) → canonical
# DAG node names.
STRATI_LABEL_MAP = {
    "C00:CD8 T": "Naive_CD8_T_cell",
    "C01:CD8 T": "Naive_CD8_T_cell",
    "C02:CD8 T": "Naive_CD8_T_cell",
    "C03:CD8 T": "Naive_CD8_T_cell",
    "C04:CD8 T": "Naive_CD8_T_cell",
    "C05:CD8 T": "Naive_CD8_T_cell",
    "C06:MAIT": "Naive_CD8_T_cell",   # MAIT not in DAG; CD8 fallback
    "C07:CD4 T": "Naive_CD4_T_cell",
    "C08:Treg": "Treg",
    "C09:CD16 NK": "NK_cell",
    "C10:NK": "NK_cell",
    "C11:Naive B": "B_cell",
    "C12:Memory B": "B_cell",
    "C13:Plasma": "Plasma_cell",
    "C14:Small pre-B": "Pre_B_cell",
    "C15:pro-B/large pre-B": "Pro_B_cell",
    "C16:CD14 Mono": "Monocyte",
    "C17:CD16 Mono": "Monocyte",
    "C18:DC": "DC",
    "C19:pDC": "pDC",
    "C20:Macrophage": "Macrophage",
    "C21:Megakaryocyte": "Megakaryocyte",
    "C23:Myeloid progenitor": "CMP",
    "C24:Erythroid progenitor": "BFU_E",
    "C25:HSC": "HSC",
    # "C22:Cycling lymphocyte" left unmapped (proliferating lymphoid pool)
}


def _dedup_var_to_symbol(a: ad.AnnData) -> ad.AnnData:
    """Replace var.index (Ensembl) with feature_name, averaging duplicates."""
    if "feature_name" not in a.var.columns:
        raise KeyError("var.feature_name missing")
    sym = a.var["feature_name"].astype(str).values
    # Mask out empty / NA symbol
    keep = np.array([s and s != "nan" and s != "None" for s in sym])
    a = a[:, keep].copy()
    sym = a.var["feature_name"].astype(str).values
    # If duplicates, aggregate by symbol via groupby-mean on a DataFrame
    if pd.Series(sym).duplicated().any():
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        df = pd.DataFrame(X.T, index=sym)
        df = df.groupby(df.index).mean()  # average duplicates
        new_var = a.var.copy()
        new_var = new_var.drop_duplicates(subset=["feature_name"])
        new_var = new_var.set_index("feature_name")
        new_var = new_var.loc[df.index]
        new_X = df.T.values  # cells x genes
        a2 = ad.AnnData(X=new_X, obs=a.obs.copy(), var=new_var, obsm=dict(a.obsm),
                        uns=dict(a.uns))
        return a2
    a.var = a.var.set_index("feature_name")
    a.var.index.name = None
    return a


def prep_triana() -> Path:
    src = EXT / "triana_young_old_bm.h5ad"
    print(f"[triana] reading {src}")
    a = ad.read_h5ad(src)
    print(f"[triana] before: {a.shape}, obsm={list(a.obsm.keys())}")

    a = _dedup_var_to_symbol(a)
    # Paper cluster column = Cell_label (verbose). Use the verbose value
    # as the cluster source so cluster resolution is preserved; add a
    # canonical (DAG-node-name) column for eval lookups.
    if "Cell_label" not in a.obs.columns:
        raise KeyError("Cell_label missing in Triana obs")
    a.obs["paper_AuthorCellType_Broad"] = a.obs["Cell_label"].astype("category")
    canon = a.obs["Cell_label"].astype(str).map(TRIANA_LABEL_MAP)
    a.obs["paper_AuthorCellType_canonical"] = canon.astype("category")
    n_mapped = int(canon.notna().sum())
    n_total = int(len(canon))
    uniq_unmapped = sorted(set(a.obs.loc[canon.isna(), "Cell_label"].astype(str)))
    print(f"[triana] label canonicalization: {n_mapped}/{n_total} mapped, "
          f"{len(uniq_unmapped)} unique unmapped labels: {uniq_unmapped[:8]}")

    # Confirm an embedding for KNN (prefer X_pca; X_mofa second)
    if "X_pca" not in a.obsm and "X_mofa" in a.obsm:
        a.obsm["X_pca"] = a.obsm["X_mofa"]
    print(f"[triana] after: {a.shape}, paper col uniq="
          f"{a.obs['paper_AuthorCellType_Broad'].nunique()}")

    out = OUT / "triana_bm_prepared.h5ad"
    a.write_h5ad(out)
    print(f"[triana] wrote {out}")
    return out


def prep_strati() -> Path:
    src = EXT / "strati2023_bmmc_healthy.h5ad"
    print(f"[strati] reading {src}")
    a = ad.read_h5ad(src)
    print(f"[strati] before: {a.shape}, obsm={list(a.obsm.keys())}")

    # Subset to Healthy donors
    if "disease_setting" in a.obs.columns:
        mask = a.obs["disease_setting"].astype(str) == "Normal"
        a = a[mask].copy()
        print(f"[strati] after Normal-subset: {a.shape}")

    a = _dedup_var_to_symbol(a)
    if "cluster" not in a.obs.columns:
        raise KeyError("cluster column missing in Strati obs")
    a.obs["paper_AuthorCellType_Broad"] = a.obs["cluster"].astype("category")
    canon = a.obs["cluster"].astype(str).map(STRATI_LABEL_MAP)
    a.obs["paper_AuthorCellType_canonical"] = canon.astype("category")
    n_mapped = int(canon.notna().sum())
    n_total = int(len(canon))
    uniq_unmapped = sorted(set(a.obs.loc[canon.isna(), "cluster"].astype(str)))
    print(f"[strati] label canonicalization: {n_mapped}/{n_total} mapped, "
          f"{len(uniq_unmapped)} unique unmapped labels: {uniq_unmapped}")

    # NOTE: Strati CXG download ships only X_umap (no X_pca / X_harmony).
    # We DO NOT recompute PCA — the paper UMAP is the canonical embedding,
    # and DCTA's KNN runs on whatever `--use-rep` points at. For Strati
    # callers should pass `--use-rep X_umap`.

    print(f"[strati] after: {a.shape}, paper col uniq="
          f"{a.obs['paper_AuthorCellType_Broad'].nunique()}")

    out = OUT / "strati_bm_prepared.h5ad"
    a.write_h5ad(out)
    print(f"[strati] wrote {out}")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--which", required=True, choices=["triana", "strati", "both"])
    args = p.parse_args()
    if args.which in ("triana", "both"):
        prep_triana()
    if args.which in ("strati", "both"):
        prep_strati()


if __name__ == "__main__":
    main()
