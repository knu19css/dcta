# DCTA — DAG-Consistent Trajectory Annotation

Phase 1: data-driven DAG (G_data) on a single BM dataset (Hay).

See `docs/specs/2026-05-21-dcta-phase1-design.md` for design.
See `docs/plans/2026-05-21-dcta-phase1.md` for implementation plan.
See `STATE.md` for current progress and parameter registry.

## Quick start (after implementation)

```bash
conda activate sc
pytest -v
python run_phase1.py --h5ad ../10_bm_autoannotation/test_datasets/cxg/hay_hca_balanced_bm.h5ad
```
