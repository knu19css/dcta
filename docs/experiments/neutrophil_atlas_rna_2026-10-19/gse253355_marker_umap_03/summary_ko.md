# GSE253355 실제 R UMAP 에 고정 6 단계 RNA 패널 대응 (표시 · 집계 전용) — 요약

- 그림 3 장은 모두 R (R version 4.5.2 (2025-10-31), base graphics). 좌표 = gse253355_umap_03 의 R uwot UMAP 그대로, 점수 = gse253355_02 군집 점수 그대로 (재계산 없음).
- 입력 sha256: 좌표 1f91bc928f06… · lossless c423ecfe33be… · 객체 83fd17576c48… · calc_manifest 5842fd04cf6e… (고정값 일치, 실행 뒤에도 같음).
- 군집 26 개의 1 위 단계 (winner): Myeloblast 11 · PM 4 · MY 0 · MM 2 · BN 2 · SN 7 · 동률 0 · 0 점 0.
- 단계별 전체 군집 최고 (peak, 저자 최다 라벨): Myeloblast 8 (HSC) · PM 6 (Neutrophil) · MY 12 (Early Myeloid Progenitor) · MM 15 (Late Myeloid) · BN 25 (Ba/Eo/Ma) · SN 1 (Adipo-MSC).
- peak 이지만 그 군집 1 위는 다른 단계: 군집 12 (peak MY, 1 위 MM) · 군집 25 (peak BN, 1 위 Myeloblast).
- 주요 mismatch 군집의 1 위 단계 기여 상위 유전자 (단계 점수 대비 비율):
  - 6 (Neutrophil 93%) 1 위 PM 0.95: G0S2 23% · FPR1 23% · AQP9 19%
  - 12 (Early Myeloid Progenitor 59%) 1 위 MM 0.51: LTF 40% · LCN2 14% · CEACAM1 12%
  - 15 (Late Myeloid 100%) 1 위 MM 1.77: LTF 27% · LCN2 20% · CAMP 18%
  - 1 (Adipo-MSC 55%) 1 위 SN 0.48: ZFP36L1 51% · ALPL 23% · MME 12%
  - 25 (Ba/Eo/Ma 98%) 1 위 Myeloblast 0.36: GATA2 46% · FCER1A 23% · RGS1 16%
- 부분집합 18,892 세포 (저자 라벨 8 개): 번호 군집 6 · 7 · 8 · 12 · 15 · 16 · 17; 번호 없는 군집 14 개 89 세포는 색만.
- 유전자 기여 = 군집 평균 log1p(CP10k) / 10, 기여 합 = 군집 점수 (최대 차 1.3e-15, 156 군집 x 단계). 없는 유전자 CCL3L3 은 발현 NA · 기여 0.
- 두 단계에 겹치는 패널 유전자: MMP9 (PM · BN), CCL4L2 (PM · SN) — 두 단계 점수에 모두 들어간다.
- 저자 라벨은 계산 고정 뒤 표시 · 사후 대조에만 썼다. 후보는 확률 · 확정 세포 유형 · 분화 판정이 아니며 정확도를 주장하지 않는다.
- 표: cluster_marker_summary · cluster_stage_gene_contributions · focus_cluster_top_genes · panel_gene_diagnostics · contribution_check. 그림: figures/fig1 · fig2 · fig3 (PNG).
- 검사 13 / 13 통과: inputs_pinned, frozen_scores_equal_lossless, cells_clusters_coordinates_unchanged, gene_values_reproduce_frozen_scores, contribution_sum_equals_cluster_score, missing_gene_na_and_zero, author_composition_matches_existing, author_colours_match_r73, figure_colours_match_winner, subset_same_coordinates, existing_inputs_unchanged, tables_roundtrip_ok, public_copy_aggregate_only.
