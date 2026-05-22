"""Aggregate summary.json from 3 BM DCTA runs into one comparison table.

Reads:
  results/phase1_hay_paperclust/summary.json   (Zeng 2025 BMM 30k subset)
  results/phase1_triana/summary.json           (Triana 2021 Nat Immunol)
  results/phase1_strati/summary.json           (Strati 2023 Cell Rep Med — Normal donors)

Writes:
  docs/cmp_3datasets.tsv   — TSV (also printed to stdout)
"""
from __future__ import annotations
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]


def load_summary(p: Path) -> dict | None:
    if not p.exists():
        return None
    with open(p) as f:
        return json.load(f)


def get(d, *keys, default=None):
    cur = d
    for k in keys:
        if cur is None or not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def main():
    runs = [
        ("Zeng 2025 BMM 30k", ROOT / "results/phase1_hay_paperclust/summary.json"),
        ("Triana 2021",        ROOT / "results/phase1_triana/summary.json"),
        ("Strati 2023 Healthy", ROOT / "results/phase1_strati/summary.json"),
    ]

    rows = []
    fields = [
        ("dataset", lambda s, name: name),
        ("n_cells", lambda s, name: s.get("n_cells_after_qc")),
        ("n_genes", lambda s, name: s.get("n_genes_after_qc")),
        ("n_clusters", lambda s, name: s.get("n_clusters")),
        ("n_components", lambda s, name: s.get("n_components")),
        ("use_rep", lambda s, name: s.get("use_rep")),
        ("n_edges (init)", lambda s, name: s.get("n_edges_initial")),
        ("n_edges (after FAS)", lambda s, name: s.get("n_edges_after_fas")),
        ("edges removed", lambda s, name: s.get("n_edges_removed")),
        ("direction method", lambda s, name: get(s, "direction", "method_chosen")),
        ("direction agreement", lambda s, name: get(s, "direction", "agreement_rate")),
        ("prior consistency", lambda s, name: get(s, "prior", "consistency_rate")),
        ("backbone len", lambda s, name: get(s, "backbone", "backbone_length")),
        ("backbone kendall_tau", lambda s, name: get(s, "backbone", "kendall_tau")),
        ("backbone monotonic", lambda s, name: get(s, "backbone", "seniority_monotonic")),
        ("DCTA strict_acc", lambda s, name: get(s, "annotation", "strict_accuracy")),
        ("DCTA lineage_acc", lambda s, name: get(s, "annotation", "lineage_accuracy")),
        ("DCTA n_evaluable", lambda s, name: get(s, "annotation", "n_evaluable")),
        ("DCTA n_unmapped_gt", lambda s, name: get(s, "annotation", "n_unmapped_gt_labels")),
        ("marker-only strict_acc", lambda s, name: get(s, "marker_only_baseline", "strict_accuracy")),
        ("marker-only lineage_acc", lambda s, name: get(s, "marker_only_baseline", "lineage_accuracy")),
        ("DCTA - marker-only (strict, pp)",
         lambda s, name: _gap(get(s, "annotation", "strict_accuracy"),
                              get(s, "marker_only_baseline", "strict_accuracy"))),
        ("clusters changed by DCTA",
         lambda s, name: get(s, "marker_only_baseline", "n_clusters_changed_by_dcta")),
    ]

    summaries = []
    for name, path in runs:
        s = load_summary(path)
        if s is None:
            print(f"[warn] {path} missing — run not yet completed", file=sys.stderr)
            summaries.append((name, {}))
        else:
            summaries.append((name, s))

    # Build columns: field × dataset
    header = ["metric"] + [name for name, _ in summaries]
    table_rows = [header]
    for label, fn in fields:
        row = [label]
        for name, s in summaries:
            v = fn(s, name) if s else None
            row.append(_fmt(v))
        table_rows.append(row)

    out = ROOT / "docs" / "cmp_3datasets.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        for r in table_rows:
            f.write("\t".join(r) + "\n")

    # Pretty print
    widths = [max(len(r[i]) for r in table_rows) for i in range(len(header))]
    for r in table_rows:
        print("  ".join(c.ljust(w) for c, w in zip(r, widths)))
    print(f"\nwrote {out}")


def _gap(a, b):
    if a is None or b is None:
        return None
    return f"{(a - b) * 100:+.2f}"


def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


if __name__ == "__main__":
    main()
