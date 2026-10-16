# DCTA — DAG-Consistent Trajectory Annotation

## 이 branch의 작업: 단일 경로–분기 경로 비교

`codex/branch-candidates-2026-10-16`에서 고정된 경로 후보의 RNA 발현을 한 곡선과 여러 갈래의 곡선으로 설명하고, 학습에 쓰지 않은 자료의 예측을 비교하는 R 기능을 개발합니다. **사용자의 결과 검토 전에는 병합하지 않습니다.**

- [진행 상태와 사용 범위](README_branch_candidates_ko.md)
- [모형·공여자 수별 평가·완료 기준](docs/IMPLEMENTATION_BRANCH_COMPARISON_2026-10-16.md)

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
