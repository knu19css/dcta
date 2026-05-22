"""Export BM hematopoiesis DAG for Cytoscape (DEFAULT/canonical export).

Node set       = slide 2 of data/DAG_editable.pptx (33 nodes, broad-level).
Marker source  = data/priors/final_unique_gene.csv (one column per cell type).
Edges          = canonical hematopoiesis topology covering the 33-node set,
                 verified against slide 2.

Outputs (docs/cytoscape/):
  bm_dag_broad.graphml
  bm_dag_broad.sif
  bm_dag_broad_nodes.tsv
  bm_dag_broad_edges.tsv

Nodes whose label does not map to any CSV column will have an empty
`markers` field and `has_marker_panel = "False"` -> in Cytoscape, map
Border Paint via Discrete Mapping (False -> red) to highlight them.
"""
from __future__ import annotations
from pathlib import Path
import csv
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "priors" / "final_unique_gene.csv"
OUT_DIR = ROOT / "docs" / "cytoscape"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Slide 2 node list (extracted from DAG_editable.pptx, slide 2).
# Order roughly top-down for readability; topology defined separately below.
# ---------------------------------------------------------------------------
SLIDE2_NODES = [
    "HSC", "MPP",
    "CMP", "CLP",
    "GMP", "MEP", "MDP",
    "BMCP", "CDP", "BFU",
    "Neutrophil", "Eosinophil", "Basophil",
    "Monocyte", "Mast_cell", "Megakaryocyte",
    "Erythroblast", "CFU",
    "DC", "pDC",
    "Pro_B_cell", "Pre_B_cell", "NK_cell", "T_cell_progenitor",
    "Macrophage", "Platelet", "Erythrocyte",
    "cDC1", "cDC2",
    "CD14_Monocyte", "CD16_Monocyte",
    "B_cell", "Plasma_cell",
]

# ---------------------------------------------------------------------------
# Mapping: slide-2 node label  ->  CSV column header.
# Nodes not listed here will get empty markers (red border in Cytoscape).
# ---------------------------------------------------------------------------
SLIDE2_TO_CSV = {
    "HSC": "HSC",
    "MPP": "MPP",
    "CMP": "CMP",
    "CLP": "CLP",
    "GMP": "GMP",
    "MEP": "MEP",
    "MDP": "MDP",
    "BMCP": "BMCP",
    "Pro_B_cell": "Pro B cell",
    "NK_cell": "NK",
    "Neutrophil": "Neutrophils",
    "Eosinophil": "Eosinophil",
    "Basophil": "Basophils",
    "Monocyte": "Mono",
    "Mast_cell": "Mast",
    "Erythroblast": "Erythroblast",
    "Pre_B_cell": "Pre B cell",
    "Macrophage": "Macrophages",
    "Platelet": "Platelates",
    "DC": "DC",
    "pDC": "pDC",
    "B_cell": "Na_ve B cell",   # slide has Plasma_cell separately -> map to naive B
    "Erythrocyte": "Reticulocytes",
    "Plasma_cell": "Plasma B",
}

# Nodes intentionally without a CSV column (highlight with red border).
EXPECTED_EMPTY = {
    "BFU", "CFU", "CDP", "Megakaryocyte",
    "cDC1", "cDC2",
    "CD14_Monocyte", "CD16_Monocyte",
    "T_cell_progenitor",
}

# ---------------------------------------------------------------------------
# Edges: canonical hematopoiesis (DAG) covering the 33-node set.
# ---------------------------------------------------------------------------
EDGES = [
    ("HSC", "MPP"),
    ("MPP", "CMP"), ("MPP", "CLP"),
    ("CMP", "GMP"), ("CMP", "MEP"), ("CMP", "MDP"), ("CMP", "BMCP"),
    ("GMP", "Neutrophil"), ("GMP", "Eosinophil"), ("GMP", "Monocyte"),
    ("BMCP", "Basophil"), ("BMCP", "Mast_cell"),
    ("MDP", "Monocyte"), ("MDP", "CDP"),
    ("CDP", "DC"),
    ("CLP", "pDC"), ("DC", "pDC"),
    ("DC", "cDC1"), ("DC", "cDC2"),
    ("Monocyte", "CD14_Monocyte"), ("Monocyte", "CD16_Monocyte"),
    ("Monocyte", "Macrophage"),
    ("MEP", "Megakaryocyte"), ("MEP", "BFU"),
    ("Megakaryocyte", "Platelet"),
    ("BFU", "CFU"),
    ("CFU", "Erythroblast"),
    ("Erythroblast", "Erythrocyte"),
    ("CLP", "Pro_B_cell"), ("CLP", "NK_cell"), ("CLP", "T_cell_progenitor"),
    ("Pro_B_cell", "Pre_B_cell"),
    ("Pre_B_cell", "B_cell"),
    ("B_cell", "Plasma_cell"),
]


def load_csv_columns(path: Path) -> dict[str, list[str]]:
    """Read final_unique_gene.csv: first row = column headers (cell types),
    subsequent rows = markers (one per row, may be blank)."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    headers = rows[0]
    out: dict[str, list[str]] = {h: [] for h in headers}
    for row in rows[1:]:
        for i, val in enumerate(row):
            v = (val or "").strip()
            if not v or i >= len(headers):
                continue
            out[headers[i]].append(v)
    # dedupe while preserving order
    for k in out:
        seen = set()
        deduped = []
        for g in out[k]:
            if g not in seen:
                deduped.append(g)
                seen.add(g)
        out[k] = deduped
    return out


def lineage_for(node: str) -> str:
    """Coarse lineage label for styling (purely cosmetic)."""
    if node in {"HSC", "MPP"}:
        return "stem"
    if node in {"CMP", "CLP", "GMP", "MEP", "MDP", "BMCP", "CDP", "BFU", "CFU"}:
        return "progenitor"
    if node in {"Neutrophil", "Eosinophil", "Basophil", "Monocyte",
                "Mast_cell", "Macrophage", "DC", "pDC", "cDC1", "cDC2",
                "CD14_Monocyte", "CD16_Monocyte"}:
        return "myeloid"
    if node in {"MEP", "Megakaryocyte", "Platelet",
                "Erythroblast", "Erythrocyte"}:
        return "erythroid_meg"
    if node in {"Pro_B_cell", "Pre_B_cell", "B_cell", "Plasma_cell"}:
        return "lymphoid_B"
    if node in {"NK_cell", "T_cell_progenitor"}:
        return "lymphoid_T_NK"
    return "other"


def main() -> None:
    csv_cols = load_csv_columns(CSV_PATH)
    print(f"CSV columns: {len(csv_cols)}")

    G: nx.DiGraph = nx.DiGraph()
    red_border_nodes = []

    for node in SLIDE2_NODES:
        csv_col = SLIDE2_TO_CSV.get(node)
        markers = csv_cols.get(csv_col, []) if csv_col else []
        attrs = {
            "name": node,
            "stage": (
                "stem" if node in {"HSC", "MPP"}
                else "progenitor" if lineage_for(node) == "progenitor"
                else "mature"
            ),
            "lineage": lineage_for(node),
            "csv_column": csv_col or "",
            "n_markers": len(markers),
            "markers": "|".join(markers),         # pipe-delim for Cytoscape List import
            "markers_display": ", ".join(markers),  # for labels/tooltips
            "has_marker_panel": str(bool(markers)),
        }
        if not markers:
            red_border_nodes.append(node)
        G.add_node(node, **attrs)

    for s, d in EDGES:
        G.add_edge(s, d, relation="differentiates_to")

    # --- write graphml -----------------------------------------------------
    gpath = OUT_DIR / "bm_dag_broad.graphml"
    nx.write_graphml(G, gpath)
    print(f"wrote {gpath}  ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")

    sif = OUT_DIR / "bm_dag_broad.sif"
    with open(sif, "w") as f:
        for s, d in G.edges():
            f.write(f"{s}\tdifferentiates_to\t{d}\n")
    print(f"wrote {sif}")

    ntsv = OUT_DIR / "bm_dag_broad_nodes.tsv"
    fields = ["name", "stage", "lineage", "csv_column",
              "n_markers", "has_marker_panel", "markers", "markers_display"]
    with open(ntsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(fields)
        for n, d in G.nodes(data=True):
            w.writerow([d.get(k, "") for k in fields])
    print(f"wrote {ntsv}")

    etsv = OUT_DIR / "bm_dag_broad_edges.tsv"
    with open(etsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source", "target", "relation"])
        for s, d, dd in G.edges(data=True):
            w.writerow([s, d, dd.get("relation", "")])
    print(f"wrote {etsv}")

    # --- topology sanity ---------------------------------------------------
    roots = [n for n, deg in G.in_degree() if deg == 0]
    leaves = [n for n, deg in G.out_degree() if deg == 0]
    print(f"\nroots: {roots}")
    print(f"leaves: {leaves}")
    print(f"is_dag: {nx.is_directed_acyclic_graph(G)}")

    print(f"\nnodes with NO markers (red border in Cytoscape) "
          f"[{len(red_border_nodes)}]: {red_border_nodes}")
    unexpected_empty = set(red_border_nodes) - EXPECTED_EMPTY
    unexpected_filled = EXPECTED_EMPTY - set(red_border_nodes)
    if unexpected_empty:
        print(f"  ! unexpectedly empty (check mapping): {sorted(unexpected_empty)}")
    if unexpected_filled:
        print(f"  ! expected empty but got markers: {sorted(unexpected_filled)}")


if __name__ == "__main__":
    main()
