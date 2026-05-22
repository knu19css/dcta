"""Rebuild DCTA prior YAMLs from data/priors/final_unique_gene.csv.

CSV layout: 1 header row of cell-type names; 12 data rows of marker genes
(empty cells allowed).

Outputs (back up originals to *.bak.csv-rebuild):
  data/priors/dag_prior_bm.yaml             — 32 nodes + topology edges
  data/priors/lineage_mature_markers.yaml   — markers from CSV (normalized)

T-maturation YAML is intentionally NOT regenerated — CSV has only one
broad thymic node ("Thymocytes"), so the per-stage DN/DP/SP subtree is
dropped for this run.
"""
from __future__ import annotations
import csv
from pathlib import Path
import shutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "priors" / "final_unique_gene.csv"
PRIOR_DIR = ROOT / "data" / "priors"


# CSV column header → canonical DAG node name.
COL_TO_NODE: dict[str, str] = {
    "HSC": "HSC",
    "MPP": "MPP",
    "CMP": "CMP",
    "CLP": "CLP",
    "GMP": "GMP",
    "MEP": "MEP",
    "MDP": "MDP",
    "BMCP": "BMCP",
    "Erythroblast": "Erythroblast",
    "Macrophages": "Macrophage",
    "Mono": "Monocyte",
    "Neutrophils": "Neutrophil",
    "DC": "DC",
    "pDC": "pDC",
    "Erythroid-like and erythroid precursor cells": "Erythroid_precursor",
    "Memory B": "Memory_B",
    "Plasma B": "Plasma_cell",
    "Na_ve B cell": "B_cell",
    "Pre B cell": "Pre_B_cell",
    "Pro B cell": "Pro_B_cell",
    "Basophils": "Basophil",
    "Eosinophil": "Eosinophil",
    "NK": "NK_cell",
    "NK_T": "NKT_cell",
    "Mast": "Mast_cell",
    "Na_ve T cells": "Naive_T",
    "CD4 T": "CD4_T",
    "Cytotoxic CD8": "Cytotoxic_CD8",
    "T reg": "Treg",
    "Platelates": "Platelet",
    "Reticulocytes": "Reticulocyte",
    "Thymocytes": "Thymocytes",
}


# Per-node stage / lineage / synthetic flags (matches old YAML semantics).
NODE_STAGE: dict[str, str] = {
    "HSC": "stem",
    "MPP": "progenitor",
    "CMP": "progenitor",
    "CLP": "progenitor",
    "GMP": "progenitor",
    "MEP": "progenitor",
    "MDP": "progenitor",
    "BMCP": "progenitor",
    "Erythroid_precursor": "progenitor",
    "Erythroblast": "progenitor",
    "Reticulocyte": "progenitor",
    "Pro_B_cell": "progenitor",
    "Pre_B_cell": "progenitor",
    "Thymocytes": "progenitor",
    "Naive_T": "mature",
    "CD4_T": "mature",
    "Cytotoxic_CD8": "mature",
    "Treg": "mature",
    "NK_cell": "mature",
    "NKT_cell": "mature",
    "B_cell": "mature",
    "Memory_B": "mature",
    "Plasma_cell": "terminal",
    "Monocyte": "mature",
    "Macrophage": "terminal",
    "DC": "mature",
    "pDC": "mature",
    "Neutrophil": "terminal",
    "Eosinophil": "terminal",
    "Basophil": "terminal",
    "Mast_cell": "terminal",
    "Platelet": "terminal",
}
NODE_LINEAGE: dict[str, str] = {
    "HSC": "stem", "MPP": "stem",
    "CMP": "myeloid", "GMP": "myeloid", "MDP": "myeloid", "BMCP": "myeloid",
    "MEP": "erythroid", "Erythroid_precursor": "erythroid",
    "Erythroblast": "erythroid", "Reticulocyte": "erythroid",
    "Platelet": "megakaryocyte",
    "Monocyte": "myeloid", "Macrophage": "myeloid",
    "Neutrophil": "myeloid", "Eosinophil": "myeloid",
    "Basophil": "myeloid", "Mast_cell": "myeloid",
    "DC": "dendritic", "pDC": "dendritic",
    "CLP": "lymphoid",
    "Pro_B_cell": "lymphoid_B", "Pre_B_cell": "lymphoid_B",
    "B_cell": "lymphoid_B", "Memory_B": "lymphoid_B", "Plasma_cell": "lymphoid_B",
    "NK_cell": "lymphoid_NK", "NKT_cell": "lymphoid_NK",
    "Thymocytes": "lymphoid_T",
    "Naive_T": "lymphoid_T", "CD4_T": "lymphoid_T",
    "Cytotoxic_CD8": "lymphoid_T", "Treg": "lymphoid_T",
}


# Per-node aliases (paper-label canonicalization helpers).
NODE_ALIASES: dict[str, list[str]] = {
    "HSC": ["hematopoietic stem cell", "Stem", "HSCs & MPPs", "HSC MPP"],
    "MPP": ["Multipotent_progenitor", "Lymphomyeloid prog", "LMPP",
            "Cycling Progenitor", "Early Lymphoid"],
    "CMP": ["Common_myeloid_progenitor", "Erythro-myeloid progenitors",
            "Myeloid progenitor", "common myeloid progenitor"],
    "CLP": ["Common_lymphoid_progenitor", "Pre-pro-B cells",
            "Lymphoid progenitor", "NK cell progenitors"],
    "GMP": ["Granulocyte-monocyte progenitor", "Early promyelocytes",
            "Late promyelocytes", "Promyelocyte",
            "Early GMP", "Late GMP", "Pro-Monocyte"],
    "MEP": ["Megakaryocyte-erythroid progenitor",
            "Megakaryocyte progenitors", "Megakaryocyte Precursor"],
    "MDP": ["Monocyte-DC progenitor"],
    "BMCP": ["Basophil-mast progenitor",
             "Eosinophil-basophil-mast cell progenitors",
             "EoBasoMast Precursor"],
    "Erythroid_precursor": ["BFU_E", "CFU_E", "BFU-E", "CFU-E",
                            "Early erythroid progenitor",
                            "Late erythroid progenitor",
                            "Early Erythroid", "Late Erythroid",
                            "erythroid progenitor cell"],
    "Erythroblast": ["Aberrant erythroid", "erythroid lineage cell"],
    "Reticulocyte": [],
    "Pro_B_cell": ["pro-B cell", "Pro-B cells", "Pro-B"],
    "Pre_B_cell": ["precursor B cell", "Pre-B cells", "Small pre-B cell",
                   "Pre-B"],
    "B_cell": ["naive B cell", "Mature naive B cells", "Naive_B",
               "Immature B cells", "immature B cell", "B"],
    "Memory_B": ["memory B cell", "Nonswitched memory B cells",
                 "Class switched memory B cells", "CD11c+ memory B cells",
                 "unswitched memory B cell", "class switched memory B cell"],
    "Plasma_cell": ["Plasma cells", "plasma cell", "Plasma_Cell"],
    "NK_cell": ["natural killer cell", "NK",
                "CD56dimCD16+ NK cells", "CD56brightCD16- NK cells",
                "CD16-positive, CD56-dim natural killer cell, human",
                "CD16-negative, CD56-bright natural killer cell, human"],
    "NKT_cell": ["NK T cells", "natural killer T cell", "iNKT"],
    "Mast_cell": ["mast cell"],
    "Basophil": ["basophil"],
    "Eosinophil": ["eosinophil"],
    "Neutrophil": ["neutrophil", "Myelocytes", "myelocyte"],
    "Monocyte": ["Classical Monocytes", "Non-classical monocytes",
                 "classical monocyte", "non-classical monocyte"],
    "Macrophage": ["macrophage"],
    "DC": ["dendritic cell", "Conventional dendritic cell 1",
           "Conventional dendritic cell 2", "cDC1", "cDC2", "cDC",
           "conventional dendritic cell"],
    "pDC": ["plasmacytoid dendritic cell",
            "Plasmacytoid dendritic cells",
            "Plasmacytoid dendritic cell progenitors",
            "common dendritic progenitor", "CDP"],
    "Megakaryocyte": [],  # not in CSV directly; left in case alias hits
    "Platelet": ["Platelets", "platelet"],
    "Thymocytes": ["T_cell_progenitor", "Pro-T", "DN", "DP",
                   "CD4_SP", "CD8_SP", "Treg_thymic", "thymocyte"],
    "Naive_T": ["CD4+ naive T cells", "CD8+ naive T cells",
                "Naive_CD4_T_cell", "Naive_CD8_T_cell", "Naive T",
                "naive thymus-derived CD4-positive, alpha-beta T cell",
                "naive thymus-derived CD8-positive, alpha-beta T cell"],
    "CD4_T": ["CD4+ memory T cells", "CD4+ cytotoxic T cells",
              "CD69+PD-1+ memory CD4+ T cells", "CD4 Memory T",
              "CD4-positive, alpha-beta memory T cell",
              "Th1", "Th2", "Th17", "Tfh"],
    "Cytotoxic_CD8": ["CD8+ effector memory T cells",
                       "CD8+ central memory T cells",
                       "CD8+CD103+ tissue resident memory T cells",
                       "effector memory CD8-positive, alpha-beta T cell",
                       "central memory CD8-positive, alpha-beta T cell",
                       "CD8_effector", "CD8_effector_memory",
                       "CD8_central_memory", "CD8_Treg", "CD8 Memory T",
                       "GammaDelta T cells", "T_gd", "gamma-delta T cell"],
    "Treg": ["regulatory T cell", "T_reg"],
}


# Gene-symbol normalization (typos / historical CD-numbers → HGNC symbols).
GENE_ALIAS: dict[str, str] = {
    "VREB1": "VPREB1",
    "VREB2": "VPREB2",
    "VREB": "VPREB1",
    "ILR7": "IL7R",
    "IL7RA": "IL7R",
    "FCERIA": "FCER1A",
    "HBA": "HBA1",
    "HBBA": "HBB",
    "CD20": "MS4A1",
    "CD21": "CR2",
    "CD23": "FCER2",
    "CD11C": "ITGAX",
    "CD11B": "ITGAM",
    "CD11A": "ITGAL",
    "CD10": "MME",
    "CD123": "IL3RA",
    "CD303": "CLEC4C",
    "CD304": "NRP1",
    "CD141": "THBD",
    "PLZF": "ZBTB16",
    "CD115": "CSF1R",
    "CD150": "SLAMF1",
    "B220": "PTPRC",
}
# Genes to drop outright (textual noise / 1-letter / not a gene)
GENE_DROP = {"", "T", "T cell", "erythroblast", "B"}


# Edges for the 32-node DAG (best-effort biological topology).
EDGES: list[tuple[str, str]] = [
    # Stem → multipotent
    ("HSC", "MPP"),
    # MPP → myeloid / lymphoid bifurcation
    ("MPP", "CMP"),
    ("MPP", "CLP"),
    # CMP → progenitors of myeloid / erythroid / mono-DC / mast
    ("CMP", "GMP"),
    ("CMP", "MEP"),
    ("CMP", "MDP"),
    ("CMP", "BMCP"),
    # GMP → granulocyte/monocyte terminals
    ("GMP", "Neutrophil"),
    ("GMP", "Eosinophil"),
    ("GMP", "Monocyte"),
    # MEP → erythroid + platelet (no Megakaryocyte node in CSV)
    ("MEP", "Erythroid_precursor"),
    ("MEP", "Platelet"),
    ("Erythroid_precursor", "Erythroblast"),
    ("Erythroblast", "Reticulocyte"),
    # MDP → DC/pDC/Monocyte
    ("MDP", "DC"),
    ("MDP", "pDC"),
    ("MDP", "Monocyte"),
    # BMCP → basophil / mast
    ("BMCP", "Basophil"),
    ("BMCP", "Mast_cell"),
    # Monocyte → Macrophage
    ("Monocyte", "Macrophage"),
    # B-cell chain
    ("CLP", "Pro_B_cell"),
    ("Pro_B_cell", "Pre_B_cell"),
    ("Pre_B_cell", "B_cell"),
    ("B_cell", "Memory_B"),
    ("Memory_B", "Plasma_cell"),
    ("B_cell", "Plasma_cell"),  # short-circuit naive → plasma also seen
    # NK lineage
    ("CLP", "NK_cell"),
    # Thymic / T branch
    ("CLP", "Thymocytes"),
    ("Thymocytes", "Naive_T"),
    ("Thymocytes", "NKT_cell"),  # iNKT branches at DP-equivalent stage
    ("Naive_T", "CD4_T"),
    ("Naive_T", "Cytotoxic_CD8"),
    ("Naive_T", "Treg"),
]


def normalize_gene(g: str) -> str | None:
    g = (g or "").strip()
    if g in GENE_DROP:
        return None
    return GENE_ALIAS.get(g.upper(), g)


def main():
    if not CSV_PATH.exists():
        raise SystemExit(f"CSV not found: {CSV_PATH}")

    # Parse CSV
    with open(CSV_PATH) as f:
        reader = csv.reader(f)
        rows = list(reader)
    header = rows[0]
    data_rows = rows[1:]

    print(f"[csv] {CSV_PATH.name}: {len(data_rows)} rows × {len(header)} cols")
    unmapped_cols = [h for h in header if h not in COL_TO_NODE]
    if unmapped_cols:
        raise SystemExit(f"unmapped CSV columns: {unmapped_cols}")

    # Build node -> deduped marker list (preserving CSV order)
    node_markers: dict[str, list[str]] = {}
    for col_idx, col_name in enumerate(header):
        node = COL_TO_NODE[col_name]
        seen: set[str] = set(node_markers.get(node, []))
        ordered = list(node_markers.get(node, []))
        for r in data_rows:
            if col_idx >= len(r):
                continue
            g = normalize_gene(r[col_idx])
            if g is None or g in seen:
                continue
            seen.add(g)
            ordered.append(g)
            # HBA / HBA1 special: also include HBA2 paralog when HBA collapsed
            if r[col_idx].strip().upper() == "HBA" and "HBA2" not in seen:
                seen.add("HBA2")
                ordered.append("HBA2")
        node_markers[node] = ordered

    # Build node set: every node we declared in COL_TO_NODE
    all_nodes = sorted(set(COL_TO_NODE.values()),
                       key=lambda x: (NODE_STAGE.get(x, "z"), x))

    # ---- Write dag_prior_bm.yaml ----
    bm_yaml_path = PRIOR_DIR / "dag_prior_bm.yaml"
    if bm_yaml_path.exists():
        bak = bm_yaml_path.with_suffix(".yaml.bak.csv-rebuild")
        shutil.copy2(bm_yaml_path, bak)
        print(f"[backup] {bm_yaml_path} → {bak.name}")

    bm_doc = {
        "nodes": all_nodes,
        "edges": [list(e) for e in EDGES],
    }
    # Pretty preamble comment
    header_comment = (
        "# Rebuilt from data/priors/final_unique_gene.csv (32 cell types).\n"
        "# Edges = best-effort hematopoietic precedence covering the 32-node\n"
        "# set. T-maturation subtree (DN/DP/SP/etc.) collapsed into a single\n"
        "# Thymocytes node since the CSV groups them broadly.\n\n"
    )
    with open(bm_yaml_path, "w") as f:
        f.write(header_comment)
        yaml.safe_dump(bm_doc, f, sort_keys=False, default_flow_style=False)
    print(f"[wrote] {bm_yaml_path}  ({len(all_nodes)} nodes, {len(EDGES)} edges)")

    # ---- Write lineage_mature_markers.yaml ----
    mk_yaml_path = PRIOR_DIR / "lineage_mature_markers.yaml"
    if mk_yaml_path.exists():
        bak = mk_yaml_path.with_suffix(".yaml.bak.csv-rebuild")
        shutil.copy2(mk_yaml_path, bak)
        print(f"[backup] {mk_yaml_path} → {bak.name}")

    mk_doc: dict[str, dict] = {}
    for node in all_nodes:
        markers = node_markers.get(node, [])
        mk_doc[node] = {
            "markers": markers,
            "aliases": NODE_ALIASES.get(node, []),
            "stage": NODE_STAGE.get(node, "mature"),
            "lineage": NODE_LINEAGE.get(node, "unspecified"),
            "synthetic": False,
        }
    mk_header = (
        "# Rebuilt from data/priors/final_unique_gene.csv.\n"
        "# `markers` are taken verbatim from the CSV (column = cell type,\n"
        "# rows = genes), with light symbol normalization (VREB→VPREB,\n"
        "# CD20→MS4A1, CD11C→ITGAX, etc.). Aliases are kept for paper-label\n"
        "# canonicalization across Zeng/Triana/Strati datasets.\n\n"
    )
    with open(mk_yaml_path, "w") as f:
        f.write(mk_header)
        yaml.safe_dump(mk_doc, f, sort_keys=False, default_flow_style=False)
    print(f"[wrote] {mk_yaml_path}  ({len(mk_doc)} nodes)")

    # ---- Summary table ----
    print("\nnode marker-panel sizes:")
    for n in all_nodes:
        print(f"  {n:25s} n={len(node_markers.get(n, [])):2d}  "
              f"{node_markers.get(n, [])[:6]}")


if __name__ == "__main__":
    main()
