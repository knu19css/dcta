# DCTA — DAG-Consistent Trajectory Annotation

## 이 branch의 작업: RNA 중심 DAG 검토안

`codex/rna-dag-review-2026-10-16`에서 HSC·MPP를 별도 후보로 유지하는 노드·연결표와 검토용 그림을 개발합니다. **사용자가 그림과 변경표를 확인하고 승인하기 전에는 병합하지 않습니다.**

- [진행 상태와 검토 안내](README_rna_dag_review_ko.md)
- [구현 범위와 완료 기준](docs/IMPLEMENTATION_RNA_DAG_REVIEW_2026-10-16.md)

## 기존 Phase 1 안내

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
