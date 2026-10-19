# GSE253355 새 UMAP (표시 전용) — 요약

- 입력: gse253355_02 lossless.rds (sha256 c423ecfe33be…) 의 저장 PCA 30 · kNN 20, 로컬 객체 (sha256 83fd17576c48…) 의 저자 주석 cluster_anno_l2.
- 계산: uwot::umap 0.2.5, n_neighbors=20, min_dist=0.3, spread=1, n_components=2, n_epochs=200, metric=euclidean, scale=FALSE, init=spectral, n_threads=2, n_sgd_threads=1, fast_sgd=FALSE, seed=42, set.seed(42). 이웃 거리는 저장 PCA 로 정확 계산(자기 0). 재척도화 · HVG · 그래프 · 점수 재계산 없음.
- 그림 1: 전체 82,742 세포, 저자 라벨 35 개 범례. 그림 2: HSC · MPP · Cycling HSPC · GMP · Early Myeloid Progenitor · Late Myeloid · Neutrophil · Monocyte (18,892 세포), 같은 좌표의 확대.
- 계통 라벨 세포 수: HSC 2,323 · MPP 1,727 · Cycling HSPC 1,108 · GMP 3,599 · Early Myeloid Progenitor 1,638 · Late Myeloid 3,294 · Neutrophil 3,857 · Monocyte 1,346.
- 이 좌표는 새로 계산한 UMAP 이다: 저자 원 UMAP 도, gse253355_02 의 DrL 레이아웃도 아니다. 연결선 · 방향 · 성숙 단계 가정을 넣지 않았다.
- 색은 저자 주석(정답 라벨 아님)이며 DCTA 단계 후보와 다르다. 세포 좌표 · 군집 번호는 results 의 umap_coordinates*.{rds,csv.gz} 에만 둔다(공개 · git 제외).
- 검사: pca_finite_dims=TRUE, cell_ids_exact_alignment=TRUE, knn_self_first=TRUE, knn_dist_exact=TRUE, labels_complete=TRUE, umap_finite_dims=TRUE, all_cells_kept=TRUE, umap_args_recorded=TRUE, subset_coordinates_identical=TRUE, existing_inputs_results_unchanged=TRUE, tables_roundtrip_ok=TRUE.
