"""Evaluate dcta_annotation against a ground-truth label column.

Phase 2.3. Cell-level accuracy + per-class metrics, plus a forgiving
hierarchical match using the lineage_mature_markers.yaml alias/lineage
field — e.g. a 'CD14-positive monocyte' GT label counts as correct if
the prediction is 'Monocyte' (its canonical) OR 'Macrophage' (same
'myeloid' lineage; configurable).
"""
from __future__ import annotations
import numpy as np
import yaml
import anndata as ad


def _norm(s):
    if s is None:
        return None
    s = str(s).strip()
    return s if s and s.lower() != "nan" else None


def _build_alias_index(lineage_marker_yaml: str) -> tuple[dict, dict]:
    """Return (alias_to_canonical, canonical_to_lineage)."""
    with open(lineage_marker_yaml) as f:
        types = yaml.safe_load(f) or {}
    alias_to_canon: dict[str, str] = {}
    canon_to_lineage: dict[str, str] = {}
    for canon, td in types.items():
        td = td or {}
        alias_to_canon[canon] = canon
        for a in (td.get("aliases") or []):
            alias_to_canon[str(a)] = canon
        canon_to_lineage[canon] = str(td.get("lineage") or "unknown")
    return alias_to_canon, canon_to_lineage


def evaluate_annotation(
    adata: ad.AnnData,
    *,
    gt_column: str = "cell_type",
    pred_column: str = "dcta_annotation",
    lineage_marker_yaml: str | None = None,
) -> dict:
    """Compute strict & lenient accuracy of pred vs GT.

    Returns:
      {
        "n_cells": int,
        "n_evaluable": int,        # both GT and pred non-null and mapped
        "strict_accuracy": float,  # canonical(pred) == canonical(gt)
        "lineage_accuracy": float, # same `lineage` field (myeloid, etc.)
        "by_class": {gt_canonical: {n, correct_strict, correct_lineage}},
        "confusion_top": list of (gt, pred, count) sorted desc,
        "unmapped_gt_labels": list[str],
        "unmapped_pred_labels": list[str],
      }
    """
    if pred_column not in adata.obs.columns:
        raise ValueError(f"prediction column {pred_column!r} not in adata.obs")
    if gt_column not in adata.obs.columns:
        raise ValueError(f"GT column {gt_column!r} not in adata.obs")

    alias_to_canon: dict = {}
    canon_to_lineage: dict = {}
    if lineage_marker_yaml:
        alias_to_canon, canon_to_lineage = _build_alias_index(lineage_marker_yaml)

    def to_canon(label):
        if label is None:
            return None
        return alias_to_canon.get(str(label), str(label))

    gt = adata.obs[gt_column].astype(object).map(_norm).to_numpy()
    pred = adata.obs[pred_column].astype(object).map(_norm).to_numpy()

    by_class: dict = {}
    confusion: dict = {}
    unmapped_gt: set = set()
    unmapped_pred: set = set()
    n_eval = 0
    n_strict = 0
    n_lineage = 0
    for g, p in zip(gt, pred):
        if g is None or p is None:
            continue
        gc = to_canon(g)
        pc = to_canon(p)
        if alias_to_canon:
            if gc not in canon_to_lineage:
                unmapped_gt.add(str(g))
            if pc not in canon_to_lineage:
                unmapped_pred.add(str(p))
        n_eval += 1
        strict_ok = gc == pc
        lin_ok = (
            strict_ok
            or (canon_to_lineage.get(gc) is not None
                and canon_to_lineage.get(gc) == canon_to_lineage.get(pc))
        )
        if strict_ok:
            n_strict += 1
        if lin_ok:
            n_lineage += 1
        rec = by_class.setdefault(gc, {"n": 0, "correct_strict": 0, "correct_lineage": 0})
        rec["n"] += 1
        rec["correct_strict"] += int(strict_ok)
        rec["correct_lineage"] += int(lin_ok)
        confusion[(gc, pc)] = confusion.get((gc, pc), 0) + 1

    confusion_top = sorted(confusion.items(), key=lambda kv: -kv[1])[:30]
    return {
        "n_cells": int(adata.n_obs),
        "n_evaluable": n_eval,
        "strict_accuracy": (n_strict / n_eval) if n_eval else None,
        "lineage_accuracy": (n_lineage / n_eval) if n_eval else None,
        "by_class": {
            k: {**v,
                "strict_rate": (v["correct_strict"] / v["n"]) if v["n"] else None,
                "lineage_rate": (v["correct_lineage"] / v["n"]) if v["n"] else None}
            for k, v in by_class.items()
        },
        "confusion_top": [[g, p, c] for (g, p), c in confusion_top],
        "unmapped_gt_labels": sorted(unmapped_gt),
        "unmapped_pred_labels": sorted(unmapped_pred),
    }
