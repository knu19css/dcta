"""Update only the `paper_AuthorCellType_canonical` obs column in prepared
h5ads to use the new 32-node canonical names (from CSV-rebuilt priors).

Avoids re-running expensive PCA on Strati. Reads/writes:
  data/prepared/triana_bm_prepared.h5ad
  data/prepared/strati_bm_prepared.h5ad
"""
from __future__ import annotations
from pathlib import Path
import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "data" / "prepared"

# New 32-node DAG canonical mapping for Triana Cell_label values.
TRIANA = {
    "HSCs & MPPs": "HSC",
    "Lymphomyeloid prog": "MPP",
    "Erythro-myeloid progenitors": "CMP",
    "Megakaryocyte progenitors": "MEP",
    "Eosinophil-basophil-mast cell progenitors": "BMCP",
    "Early erythroid progenitor": "Erythroid_precursor",
    "Late erythroid progenitor": "Erythroid_precursor",
    "Aberrant erythroid": "Erythroblast",
    "Pre-pro-B cells": "CLP",
    "NK cell progenitors": "CLP",
    "Pro-B cells": "Pro_B_cell",
    "Pre-B cells": "Pre_B_cell",
    "Small pre-B cell": "Pre_B_cell",
    "Immature B cells": "B_cell",
    "Mature naive B cells": "B_cell",
    "Nonswitched memory B cells": "Memory_B",
    "Class switched memory B cells": "Memory_B",
    "CD11c+ memory B cells": "Memory_B",
    "Plasma cells": "Plasma_cell",
    "GammaDelta T cells": "Cytotoxic_CD8",
    "NK T cells": "NKT_cell",
    "CD56dimCD16+ NK cells": "NK_cell",
    "CD56brightCD16- NK cells": "NK_cell",
    "CD4+ naive T cells": "Naive_T",
    "CD8+ naive T cells": "Naive_T",
    "CD8+ central memory T cells": "Cytotoxic_CD8",
    "CD8+ effector memory T cells": "Cytotoxic_CD8",
    "CD8+CD103+ tissue resident memory T cells": "Cytotoxic_CD8",
    "CD4+ memory T cells": "CD4_T",
    "CD4+ cytotoxic T cells": "CD4_T",
    "CD69+PD-1+ memory CD4+ T cells": "CD4_T",
    "Classical Monocytes": "Monocyte",
    "Non-classical monocytes": "Monocyte",
    "Early promyelocytes": "GMP",
    "Late promyelocytes": "GMP",
    "Myelocytes": "Neutrophil",
    "Conventional dendritic cell 1": "DC",
    "Conventional dendritic cell 2": "DC",
    "Plasmacytoid dendritic cells": "pDC",
    "Plasmacytoid dendritic cell progenitors": "pDC",
    # left unmapped:
    # "Mesenchymal cells_1", "Mesenchymal cells_2", "Monocyte-like blasts"
}

# Strati cluster column → new canonical names.
STRATI = {
    "C00:CD8 T": "Cytotoxic_CD8",
    "C01:CD8 T": "Cytotoxic_CD8",
    "C02:CD8 T": "Cytotoxic_CD8",
    "C03:CD8 T": "Cytotoxic_CD8",
    "C04:CD8 T": "Cytotoxic_CD8",
    "C05:CD8 T": "Cytotoxic_CD8",
    "C06:MAIT": "Cytotoxic_CD8",
    "C07:CD4 T": "CD4_T",
    "C08:Treg": "Treg",
    "C09:CD16 NK": "NK_cell",
    "C10:NK": "NK_cell",
    "C11:Naive B": "B_cell",
    "C12:Memory B": "Memory_B",
    "C13:Plasma": "Plasma_cell",
    "C14:Small pre-B": "Pre_B_cell",
    "C15:pro-B/large pre-B": "Pro_B_cell",
    "C16:CD14 Mono": "Monocyte",
    "C17:CD16 Mono": "Monocyte",
    "C18:DC": "DC",
    "C19:pDC": "pDC",
    "C20:Macrophage": "Macrophage",
    "C21:Megakaryocyte": "MEP",
    "C23:Myeloid progenitor": "CMP",
    "C24:Erythroid progenitor": "Erythroid_precursor",
    "C25:HSC": "HSC",
    # C22:Cycling lymphocyte left unmapped
}


def relabel(h5_path: Path, source_col: str, mapping: dict[str, str], ds_name: str):
    a = ad.read_h5ad(h5_path)
    if source_col not in a.obs.columns:
        print(f"[{ds_name}] source column {source_col!r} missing; skip")
        return
    src = a.obs[source_col].astype(str)
    canon = src.map(mapping)
    a.obs["paper_AuthorCellType_canonical"] = pd.Categorical(canon)
    n_mapped = int(canon.notna().sum())
    n_total = int(len(canon))
    uniq_unmapped = sorted(set(src[canon.isna()]))
    print(f"[{ds_name}] {n_mapped}/{n_total} mapped, "
          f"unmapped={uniq_unmapped[:8]}")
    a.write_h5ad(h5_path)
    print(f"[{ds_name}] wrote {h5_path}")


def main():
    relabel(PREP / "triana_bm_prepared.h5ad",
            source_col="Cell_label", mapping=TRIANA, ds_name="triana")
    relabel(PREP / "strati_bm_prepared.h5ad",
            source_col="cluster", mapping=STRATI, ds_name="strati")


if __name__ == "__main__":
    main()
