# =============================================================================
# test_gse253355_application.R — GSE253355 (성숙 호중구 포함 골수) 에 고정 atlas 6x10 패널을 군집 단위로 다시 적용 (R68–72) 테스트
#   계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_gse253355_instruction.md
#   Seurat 객체(SeuratObject 미설치)를 슬롯 그대로 읽기 · 이중 gzip, 정규화 분모 = 모든 유전자 library, 0 · 동률 · 결측 유전자(고정 분모 vs 관측 유전자),
#   donor 동일 가중(>= 5 세포, 대체 없음), 그래프에 없는 edge 를 만들지 않음 · 같은 군집 최고 = 자기 전이 없음, 세포 ID 정렬, 실행기 순서 · 덮어쓰기 거절
# =============================================================================
library(testthat)
gsx_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- gsx_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R",
            "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R",
            "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R", "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R",
            "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R", "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R",
            "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

# 실제 객체와 같은 슬롯 이름의 가짜 Seurat 클래스 (코드는 클래스 정의가 아니라 슬롯 속성으로 읽는다)
suppressMessages({
  setClass("LogMap", contains = "matrix")
  setClass("Graph", contains = "dgCMatrix", representation(assay.used = "character"))
  setClass("Assay5", representation(layers = "list", cells = "LogMap", features = "LogMap", key = "character"))
  setClass("Assay", representation(counts = "dgCMatrix", data = "dgCMatrix", key = "character"))
  setClass("DimReduc", representation(cell.embeddings = "matrix", key = "character"))
  setClass("Seurat", representation(assays = "list", meta.data = "data.frame", active.assay = "character", graphs = "list", reductions = "list", project.name = "character"))
})
gsx_fake <- function(counts, meta, graph = NULL, umap = NULL, v5 = TRUE, split = NULL) {
  genes <- rownames(counts); cells <- colnames(counts); unm <- function(m) { dimnames(m) <- list(NULL, NULL); m }
  a <- if (v5) {
    grp <- if (is.null(split)) rep("counts", length(cells)) else paste0("counts.", split); ln <- unique(grp)
    lay <- stats::setNames(lapply(ln, function(l) unm(counts[, grp == l, drop = FALSE])), ln)
    fm <- matrix(TRUE, length(genes), length(ln), dimnames = list(genes, ln)); cm <- vapply(ln, function(l) grp == l, logical(length(cells))); dim(cm) <- c(length(cells), length(ln)); dimnames(cm) <- list(cells, ln)
    new("Assay5", layers = lay, cells = new("LogMap", cm), features = new("LogMap", fm), key = "rna_")
  } else new("Assay", counts = counts, data = counts, key = "rna_")
  new("Seurat", assays = list(RNA = a), meta.data = meta, active.assay = "RNA", graphs = if (is.null(graph)) list() else list(RNA_snn = new("Graph", graph, assay.used = "RNA")),
      reductions = if (is.null(umap)) list() else list(UMAP_dim30 = new("DimReduc", cell.embeddings = umap, key = "UMAPdim30_")), project.name = "t")
}
spm <- function(m) methods::as(methods::as(Matrix::Matrix(m, sparse = TRUE), "generalMatrix"), "CsparseMatrix")   # 대칭 자동 변환(dsC) 없이 dgC

test_that("RDS 읽기: gzip · 이중 gzip · 무압축을 모두 같은 객체로 읽고 형식을 기록한다", {
  d <- tempfile("gsx_rds"); dir.create(d); obj <- list(a = 1:3, b = "x")
  p1 <- file.path(d, "o.rds.gz"); saveRDS(obj, p1); p0 <- file.path(d, "o0.rds"); saveRDS(obj, p0, compress = FALSE)
  p2 <- file.path(d, "o2.rds.gz"); con <- gzfile(p2, "wb"); writeBin(readBin(p1, "raw", file.size(p1)), con); close(con)
  r1 <- gsx_open_rds(p1); r2 <- gsx_open_rds(p2); r0 <- gsx_open_rds(p0)
  expect_identical(r1$object, obj); expect_identical(r2$object, obj); expect_identical(r0$object, obj)
  expect_equal(c(r1$format, r2$format, r0$format), c("gzip", "double_gzip", "none"))
  expect_error(gsx_open_rds(paste0(p1, ".part")), "part")
})

test_that("raw count 추출: Seurat v5 (단일 · 표본별 분할 layer) 와 v4 가 같은 유전자 x 세포 행렬, 세포 순서 = meta.data 순서; 분할 layer 의 유전자가 다르면 멈춘다", {
  cnt <- spm(matrix(c(1, 0, 2, 0, 3, 0, 0, 4, 5, 1, 0, 2), 3, 4, dimnames = list(c("G1", "G2", "G3"), paste0("c", 1:4))))
  meta <- data.frame(donor = c("a", "b", "a", "b"), row.names = paste0("c", 1:4))
  for (ob in list(gsx_fake(cnt, meta), gsx_fake(cnt, meta, split = meta$donor), gsx_fake(cnt, meta, v5 = FALSE))) {
    x <- gsx_counts(ob, "RNA"); expect_equal(as.matrix(x$counts), as.matrix(cnt)); expect_equal(x$layers, x$layers) }
  expect_equal(gsx_counts(gsx_fake(cnt, meta, split = meta$donor), "RNA")$layers, c("counts.a", "counts.b"))
  bad <- gsx_fake(cnt, meta, split = meta$donor); f <- attr(attr(bad, "assays")$RNA, "features"); f[3, 2] <- FALSE
  a <- attr(bad, "assays")$RNA; attr(a, "features") <- f; L <- attr(a, "layers"); L$counts.b <- L$counts.b[1:2, ]; attr(a, "layers") <- L
  attr(bad, "assays") <- list(RNA = a); expect_error(gsx_counts(bad, "RNA"), "유전자")
})

test_that("그래프 · UMAP · 메타데이터: 슬롯 그대로 꺼내고 대칭 · 음수 · 대각 규칙을 기록하며, 세포 ID 가 어긋나면 멈춘다", {
  cnt <- spm(diag(3) + 1); dimnames(cnt) <- list(c("G1", "G2", "G3"), c("c1", "c2", "c3")); meta <- data.frame(cl = c("A", "A", "B"), row.names = c("c1", "c2", "c3"))
  g <- spm(matrix(c(1, 0.5, 0, 0.5, 1, 0.2, 0, 0.2, 1), 3)); dimnames(g) <- list(c("c1", "c2", "c3"), c("c1", "c2", "c3"))
  U <- cbind(UMAPdim30_1 = 1:3, UMAPdim30_2 = 3:1); rownames(U) <- c("c1", "c2", "c3"); ob <- gsx_fake(cnt, meta, graph = g, umap = U)
  G <- gsx_graph(ob, "RNA_snn"); expect_true(methods::is(G, "dgCMatrix")); expect_equal(dimnames(G), dimnames(g))
  ck <- gsx_graph_checks(G); expect_true(ck$symmetric); expect_true(ck$nonnegative); expect_equal(ck$n_diagonal, 3L); expect_equal(ck$n_upper_edges, 2L)
  expect_equal(unname(gsx_umap(ob, "UMAP_dim30")[, 1]), 1:3)
  expect_true(gsx_check_ids(c("c1", "c2", "c3"), meta = meta, umap = U, graph = G))
  U2 <- U[c(2, 1, 3), ]; expect_error(gsx_check_ids(c("c1", "c2", "c3"), meta = meta, umap = U2, graph = G), "세포 ID")
  expect_error(gsx_graph(ob, "RNA_nn"), "RNA_nn")
})

test_that("패널 유전자 대응: 정확한 symbol 만 쓴다; 없는 유전자는 missing 으로 남기고 HGNC 별칭 · 접미사 변형은 진단으로만 적는다", {
  hg <- data.frame(symbol = c("CCL3L3", "ID1", "AVP"), prev_symbol = c("", "", ""), alias_symbol = c("CCL3L1", "ID", "MMP9"), stringsAsFactors = FALSE)   # AVP 는 정확 대응 → 별칭 진단 안 함
  pm <- gsx_panel_map(c("ID1", "AVP", "CCL3L1", "MMP9", "MMP9.1"), data.frame(stage = c("Myeloblast", "Myeloblast", "SN", "PM"), rank = c(1L, 2L, 1L, 1L),
                                                                              symbol = c("ID1", "AVP", "CCL3L3", "MMP9"), stringsAsFactors = FALSE), hg)
  expect_equal(pm$status, c("exact", "exact", "missing", "exact")); expect_equal(pm$feature_index[3], NA_integer_)
  expect_equal(pm$alias_candidates_in_data[3], "CCL3L1"); expect_equal(pm$suffix_variants_in_data[4], "MMP9.1"); expect_equal(pm$feature_index[c(1, 2, 4)], c(1L, 2L, 4L))
  expect_equal(pm$alias_candidates_in_data[2], "")                                                        # 정확히 측정된 유전자에는 다른 유전자 이름을 별칭으로 붙이지 않는다
})

test_that("세포 점수: library 는 모든 유전자 합, 고정 분모(패널 유전자 수) vs 관측 유전자 분모를 따로, 패널 유전자가 하나도 없으면 결측", {
  cnt <- spm(matrix(c(5, 0, 95, 5, 0, 195, 0, 0, 100), 3, 3, dimnames = list(c("g1", "g2", "X"), c("c1", "c2", "c3"))))
  pan <- data.frame(stage = c("Myeloblast", "Myeloblast", "PM"), rank = c(1L, 2L, 1L), symbol = c("g1", "gMiss", "gMiss2"), stringsAsFactors = FALSE)
  pm <- gsx_panel_map(rownames(cnt), pan, NULL); sc <- gsx_cell_scores(cnt, pm, stages = c("Myeloblast", "PM"))
  expect_equal(sc$library, c(100, 200, 100))                                                              # X(패널 밖) 도 library 에 들어간다
  expect_equal(unname(sc$fixed[, "Myeloblast"]), c(log1p(500) / 2, log1p(250) / 2, 0)); expect_equal(unname(sc$available[, "Myeloblast"]), c(log1p(500), log1p(250), 0))
  expect_true(all(is.na(sc$fixed[, "PM"]))); expect_true(all(is.na(sc$available[, "PM"]))); expect_equal(unname(sc$n_detected[, "Myeloblast"]), c(1, 1, 0))
  expect_equal(sc$stage_genes$n_listed, c(2L, 1L)); expect_equal(sc$stage_genes$n_measured, c(1L, 0L))
  expect_error(gsx_cell_scores(spm(matrix(c(0, 0, 0, 1), 2, 2, dimnames = list(c("g1", "X"), c("a", "b")))), gsx_panel_map(c("g1", "X"), pan[1, ], NULL), "Myeloblast"), "library")
})

test_that("군집 점수: 0 인 세포를 포함한 모든 세포 평균, donor 동일 가중은 donor x 군집 >= 5 세포만 (없으면 결측, 대체 없음), 검출 비율", {
  v <- c(1, 2, 0, 3, 3, 3, 3, 3, 9, 9); cl <- c("A", "A", "A", "B", "B", "B", "B", "B", "C", "C"); dn <- c("d1", "d1", "d1", "d1", "d1", "d1", "d1", "d1", "d2", "d2")
  M <- cbind(Myeloblast = v, PM = 0); N <- cbind(Myeloblast = c(1, 1, 0, 2, 2, 2, 2, 2, 1, 1), PM = 0)
  cs <- gsx_cluster_scores(M, N, cl, dn, min_cells = 5L)
  a <- cs[cs$cluster == "A" & cs$stage == "Myeloblast", ]; expect_equal(a$mean_all_cells, 1); expect_equal(a$frac_cells_ge1, 2 / 3); expect_true(is.na(a$mean_donor_equal))
  b <- cs[cs$cluster == "B" & cs$stage == "Myeloblast", ]; expect_equal(c(b$mean_donor_equal, b$n_donors, b$n_donors_min_cells), c(3, 1, 1)); expect_equal(b$frac_cells_ge2, 1)
  expect_equal(cs$n_cells[cs$cluster == "C" & cs$stage == "PM"], 2L); expect_equal(cs$mean_all_cells[cs$cluster == "C" & cs$stage == "PM"], 0)
})

test_that("그래프 감사 (Codex 검토): 실제 edge 만 세고, 군집 요약 경로와 실제 세포 경로를 따로 적는다 — B 가 끊기면 cluster_path_only, 이어지면 cell_path", {
  cells <- paste0("c", 1:6); cl <- c("A", "A", "B", "B", "C", "C"); dn <- c("d1", "d1", "d1", "d2", "d2", "d2")
  mkG <- function(i, j) { G <- Matrix::sparseMatrix(i = i, j = j, x = 0.5, dims = c(6, 6), symmetric = TRUE); G <- methods::as(methods::as(G, "generalMatrix"), "CsparseMatrix"); dimnames(G) <- list(cells, cells); G }
  e <- gsx_upper_edges(mkG(c(1, 2, 4), c(2, 3, 5))); expect_equal(nrow(e), 3L)
  ce <- gsx_cluster_edges(e, cl, dn); expect_equal(sum(ce$n_cell_edges), 2L); expect_equal(ce$n_cell_edges[ce$cluster_a == "A" & ce$cluster_b == "B"], 1L)
  expect_false(any(ce$cluster_a == "A" & ce$cluster_b == "C"))                                             # 실제 edge 가 없는 쌍은 표에 없다
  peaks <- data.frame(stage = c("Myeloblast", "PM", "MY"), peak_clusters = c("A", "A", "C"), stringsAsFactors = FALSE)
  pa <- gsx_peak_paths(e, cl, dn, peaks, length(cells))
  s1 <- pa[pa$from_stage == "Myeloblast" & pa$to_stage == "PM", ]; expect_equal(s1$relation, "same_cluster"); expect_true(is.na(s1$quotient_hops)); expect_true(is.na(s1$cells_connected))
  r <- pa[pa$from_stage == "PM" & pa$to_stage == "MY", ]                                                   # B 안 c3–c4 가 끊겼다: 군집 요약에는 A>B>C, 실제 세포 경로는 없다
  expect_equal(c(r$relation, r$direct_cell_edges, r$quotient_status, r$quotient_hops, r$quotient_path), c("different_clusters", "0", "quotient_path", "2", "A>B>C"))
  expect_false(r$cells_connected); expect_true(is.na(r$cell_hops)); expect_equal(r$path_agreement, "cluster_path_only")
  expect_equal(c(r$from_cells_in_graph, r$from_cluster_cells, r$to_cells_in_graph, r$to_cluster_cells), c(2L, 2L, 2L, 2L))
  e2 <- gsx_upper_edges(mkG(c(1, 2, 3, 4), c(2, 3, 4, 5))); r2 <- gsx_peak_paths(e2, cl, dn, peaks, 6L); r2 <- r2[r2$from_stage == "PM" & r2$to_stage == "MY", ]
  expect_true(r2$cells_connected); expect_equal(r2$cell_hops, 3L); expect_equal(r2$cell_path_clusters, "A>B>C"); expect_equal(r2$path_agreement, "cell_path")
  expect_equal(attr(gsx_peak_paths(e2, cl, dn, peaks, 6L), "cell_paths")[["PM>MY"]], 2:5)
  pg <- gsx_peak_paths(e2, cl, dn, peaks, 6L, keep = cl != "B"); g1 <- pg[pg$from_stage == "PM" & pg$to_stage == "MY", ]   # B 를 뺀 유도 그래프: 둘 다 없다
  expect_equal(c(g1$quotient_status, g1$path_agreement), c("no_quotient_path", "no_path")); expect_false(g1$cells_connected)
  k2 <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE); g2 <- gsx_peak_paths(e2, cl, dn, peaks, 6L, keep = k2); g2 <- g2[g2$from_stage == "PM" & g2$to_stage == "MY", ]
  expect_equal(c(g2$to_cells_in_graph, g2$to_cluster_cells), c(1L, 2L))                                  # 최고 군집 C 의 2 세포 중 1 개만 유도 그래프에 있다
  expect_equal(gsx_components(e, length(cells))$n_components, 3L); expect_equal(gsx_components(e2, length(cells))$n_components, 2L)   # {1,2,3} {4,5} {6} / {1..5} {6}
})

test_that("새 표현 그래프 (원 그래프 · 군집이 없을 때): vst 형 HVG → 척도화 PCA → 자기 포함 kNN → SNN → Louvain(seed) → DrL 표시 좌표; 패널 · 라벨을 받지 않는다", {
  set.seed(1); G <- 50; n <- 90; grp <- rep(1:3, each = 30); lam <- matrix(1, G, n); for (k in 1:3) lam[(k - 1) * 5 + 1:5, grp == k] <- 20
  lam[16:35, ] <- (20 * 30 + 60) / 90                                                                     # 같은 평균의 하우스키핑 유전자: vst 추세가 Poisson 분산을 잡게
  X <- spm(matrix(stats::rpois(G * n, lam), G, n, dimnames = list(sprintf("g%02d", 1:G), sprintf("c%02d", 1:n))))
  hv <- gsx_hvg(X, n = 15L); expect_equal(length(hv$genes), 15L); expect_true(all(sprintf("g%02d", 1:15) %in% hv$genes))
  expect_error(suppressWarnings(gsx_hvg(X[c(1, 6, 11, 16), ], n = 2L)), "loess")                                             # 유전자 4 개로는 추세를 못 만든다: 조용히 고르지 않는다
  pc <- gsx_pca(X, hv$genes, n_pcs = 5L); expect_equal(dim(pc$scores), c(90L, 5L)); expect_identical(pc, gsx_pca(X, hv$genes, n_pcs = 5L))
  nn <- gsx_knn(pc$scores, k = 10L, chunk = 37L); expect_equal(nn[, 1], 1:90)
  d <- as.matrix(stats::dist(pc$scores)); expect_equal(sort(nn[7, -1]), sort(order(d[7, ])[2:10]))
  S <- gsx_snn(nn, prune = 1 / 15); expect_true(Matrix::isSymmetric(S)); expect_true(all(Matrix::diag(S) == 1))
  cl <- gsx_louvain(S, resolution = 0.5, seed = 42L); expect_equal(length(unique(cl)), 3L); expect_equal(sum(table(cl, grp) > 0), 3L)
  expect_identical(cl, gsx_louvain(S, resolution = 0.5, seed = 42L))
  lay <- gsx_layout(S, seed = 42L); expect_equal(dim(lay), c(90L, 2L)); expect_true(all(is.finite(lay)))
  for (f in list(gsx_hvg, gsx_pca, gsx_knn, gsx_snn, gsx_louvain, gsx_layout)) expect_false(any(grepl("panel|label|stage|author", names(formals(f)))))
})

test_that("SNN 손계산: 자기 포함 이웃 집합의 Jaccard", {
  nn <- rbind(c(1L, 2L), c(2L, 1L), c(3L, 2L)); S <- as.matrix(gsx_snn(nn, prune = 1 / 15))
  expect_equal(unname(S), rbind(c(1, 1, 1 / 3), c(1, 1, 1 / 3), c(1 / 3, 1 / 3, 1)))
})

test_that("순위 · 단계별 최고 군집은 R66 순수 함수를 그대로 쓴다 (0 · 동률 · 일부 결측)", {
  S <- rbind(A = c(0, 0, 0, 0, 0, 0), B = c(0.2, 0.2, 0, 0, 0, 0), C = c(0.1, NA, 0, 0, 0, 0.3)); colnames(S) <- NSO_STAGES
  r <- atlas_cl_rank(S, ATLAS_CL_TIE_TOL); expect_equal(r$score_top_status, c("unassigned_zero", "tied", "score_top_candidate"))
  expect_equal(r$comparison_scope[3], "incomplete_panel_comparison")
})


# ── 소형 합성 Seurat 객체로 실행기 전체 ────────────────────────────────────────
gsx_world <- function() {
  w <- tempfile("gsx_world"); dir.create(w); out <- file.path(w, "out", "gse253355_01"); pub <- file.path(w, "pub", "gse253355_01")
  sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  # 패널: 6 단계 x 2 유전자. SN 둘째 유전자 S6b 는 자료에 없다 (불완전 비교)
  pan <- data.frame(stage = rep(NSO_STAGES, each = 2), rank = rep(1:2, 6), ensg = sprintf("ENSG%011d", 1:12), symbol = c("M1", "M2", "P1", "P2", "Y1", "Y2", "N1", "N2", "B1", "B2", "S1", "S6b"),
                    min_adjacent_effect = 1, effect_vs_prev = 1, effect_vs_next = -1, stringsAsFactors = FALSE)
  att <- file.path(w, "attempt"); dir.create(file.path(att, "source"), recursive = TRUE); atlas_write_csv(pan, file.path(att, "source", "source_panels.csv"))
  writeLines(jsonlite::toJSON(list(phase = "test", files = list(`source/source_panels.csv` = sha(file.path(att, "source", "source_panels.csv")))), auto_unbox = TRUE), file.path(att, "source_freeze.json"))
  pf <- file.path(w, "atlas_stage_markers.csv"); file.copy(file.path(att, "source", "source_panels.csv"), pf)
  genes <- c(setdiff(pan$symbol, "S6b"), "OTHER"); n <- 24; cells <- sprintf("H%02d_cell%02d", rep(1:2, each = 12), 1:24)
  cl <- rep(c("0", "1", "2", "3"), each = 6); dn <- rep(c("H01", "H02"), times = 12); au <- rep(c("GMP", "Neutrophil", "Monocyte", "CD4+ T-Cell"), each = 6)
  X <- matrix(0, length(genes), n, dimnames = list(genes, cells))
  X["M1", cl == "0"] <- 3; X["P1", cl == "0"] <- 1; X["S1", cl == "1"] <- 4; X["B1", cl == "1"] <- 2; X["S1", cl == "2"] <- 5; X["OTHER", ] <- 50
  meta <- data.frame(orig.ident = dn, seurat_clusters = factor(cl), cluster_anno_l2 = au, row.names = cells, stringsAsFactors = FALSE)
  ii <- c(1:5, 6, 7:11, 12, 13:17, 19:23); jj <- c(2:6, 7, 8:12, 13, 14:18, 20:24)                       # 군집 안 사슬 + 0–1, 1–2 사이 edge 하나씩, 3 은 따로
  G <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n), symmetric = TRUE); G <- methods::as(methods::as(G, "generalMatrix"), "CsparseMatrix"); dimnames(G) <- list(cells, cells)
  U <- cbind(UMAPdim30_1 = rep(1:4, each = 6) * 5 + seq_len(n) %% 3, UMAPdim30_2 = seq_len(n) %% 4); rownames(U) <- cells
  ob <- gsx_fake(spm(X), meta, graph = G, umap = U, split = dn); inp <- file.path(w, "GSE253355_test.rds.gz"); saveRDS(ob, inp)
  hg <- file.path(w, "hgnc.txt"); utils::write.table(data.frame(symbol = "S6b", prev_symbol = "", alias_symbol = "S6", status = "Approved"), hg, sep = "\t", row.names = FALSE, quote = FALSE)
  pins <- list(input = sha(inp), source_freeze = sha(file.path(att, "source_freeze.json")))
  argv <- c("--input", inp, "--panels", pf, "--panel-attempt", att, "--hgnc", hg, "--out-dir", out, "--public-dir", pub, "--assay", "RNA", "--cluster-col", "seurat_clusters",
            "--donor-col", "orig.ident", "--author-col", "cluster_anno_l2", "--umap", "UMAP_dim30", "--graph", "RNA_snn", "--granulocyte-labels", "GMP;Neutrophil",
            "--neutrophil-labels", "Neutrophil")
  list(argv = argv, out = out, pub = pub, pins = pins, inp = inp, cells = cells)
}
gsx_run <- function(w, pins = w$pins, argv = w$argv) { code <- NA_integer_; utils::capture.output(code <- gsx_main(argv, root = ROOT, pins = pins)); code }
gsx_rd <- function(d, f) utils::read.csv(file.path(d, f), colClasses = "character", na.strings = character(0), check.names = FALSE)

test_that("실행기 (소형 합성 Seurat): 객체 점검 → 고정 패널 점수 → 계산 고정 → (그 뒤) 저자 라벨 · 과립구 그래프 · 그림 · 공개 사본", {
  w <- gsx_world(); expect_equal(gsx_run(w), 0L); out <- w$out
  need <- c("object_inspection.json", "object_metadata_columns.csv", "panel_gene_map.csv", "cluster_scores.csv", "cluster_candidates.csv", "stage_peak_clusters.csv",
            "cluster_edges.csv", "peak_paths.csv", "graph_checks.csv", "calc_manifest.json", "posthoc_author_composition.csv", "posthoc_donor_composition.csv",
            "posthoc_peak_paths_granulocyte.csv", "posthoc_mature_neutrophil.csv", "checks.csv", "input_sha256.csv", "environment.csv", "source_snapshot.csv", "manifest.json",
            "summary_ko.md", file.path("figures", c("fig1_umap_author_vs_candidate.png", "fig1_umap_author_vs_candidate.pdf", "fig2_granulocyte_zoom_heatmap.png",
                                                    "fig2_granulocyte_zoom_heatmap.pdf")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  pm <- gsx_rd(out, "panel_gene_map.csv"); expect_equal(pm$status[pm$symbol == "S6b"], "missing"); expect_equal(pm$alias_candidates_in_data[pm$symbol == "S6b"], "")
  cd <- gsx_rd(out, "cluster_candidates.csv"); expect_equal(nrow(cd), 4L); expect_equal(cd$score_top_candidate[cd$cluster == "0"], "Myeloblast")
  expect_equal(cd$comparison_scope[cd$cluster == "0"], "incomplete_panel_comparison"); expect_equal(cd$score_top_status[cd$cluster == "3"], "unassigned_zero")
  sc <- gsx_rd(out, "cluster_scores.csv"); s1 <- sc[sc$cluster == "1" & sc$stage == "SN", ]
  expect_equal(as.numeric(s1$mean_all_cells), log1p(1e4 * 4 / 56) / 2, tolerance = 1e-12); expect_equal(as.numeric(s1$mean_available_genes), log1p(1e4 * 4 / 56), tolerance = 1e-12)
  pk <- gsx_rd(out, "stage_peak_clusters.csv"); expect_equal(pk$peak_clusters[pk$stage == "SN"], "2"); expect_equal(pk$peak_clusters[pk$stage == "Myeloblast"], "0")
  pp <- gsx_rd(out, "peak_paths.csv"); r <- pp[pp$from_stage == "Myeloblast" & pp$to_stage == "PM", ]; expect_equal(r$relation, "same_cluster")
  expect_true(all(pp$quotient_hops[pp$quotient_status == "quotient_path"] != "NA")); expect_true(all(pp$cells_connected[pp$relation == "different_clusters"] %in% c("TRUE", "FALSE")))
  b <- pp[pp$from_stage == "BN" & pp$to_stage == "SN", ]; expect_equal(c(b$path_agreement, b$cell_path_clusters), c("cell_path", "1>2"))
  ce <- gsx_rd(out, "cluster_edges.csv"); expect_equal(sum(as.integer(ce$n_cell_edges)), 2L)
  gp <- gsx_rd(out, "posthoc_peak_paths_granulocyte.csv"); expect_true(all(gp$relation %in% c("same_cluster", "different_clusters", "peak_not_in_subset", "unavailable")))
  expect_equal(gp$relation[gp$from_stage == "BN" & gp$to_stage == "SN"], "peak_not_in_subset")
  mn <- gsx_rd(out, "posthoc_mature_neutrophil.csv"); expect_equal(mn$n_cells[mn$author_label == "Neutrophil"], "6")
  ck <- gsx_rd(out, "checks.csv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  cmf <- names(jsonlite::fromJSON(file.path(out, "calc_manifest.json"))$files); expect_false(any(grepl("posthoc|figure", cmf)))
  man <- jsonlite::fromJSON(file.path(out, "manifest.json"))$files
  expect_true(all(vapply(names(man), function(f) identical(digest::digest(file.path(out, f), algo = "sha256", file = TRUE), man[[f]]), TRUE)))
  expect_lte(length(readLines(file.path(out, "summary_ko.md"))), 25L)
  pubf <- list.files(w$pub, recursive = TRUE); expect_true(all(c("summary_ko.md", "cluster_candidates.csv", "stage_peak_clusters.csv") %in% pubf))
  expect_false(any(grepl("cell_level|lossless", pubf)))                                                   # 공개 사본에는 세포 단위 · 원 자료 파일이 없다
  ft <- gsx_rd(out, "figure_text.csv"); expect_false(any(ft$text %in% c("NA", "NaN", "Inf")))
})

test_that("실행기: .part 입력 · hash 불일치 · 덮어쓰기 · 인자 누락은 멈춘다", {
  w <- gsx_world(); a2 <- w$argv; a2[2] <- paste0(w$inp, ".part"); file.copy(w$inp, a2[2]); expect_equal(gsx_run(w, argv = a2), 4L)
  w2 <- gsx_world(); bad <- w2$pins; bad$input <- strrep("0", 64); expect_equal(gsx_run(w2, pins = bad), 4L); expect_true(any(grepl("hash", readLines(file.path(w2$out, "run.log")))))
  w3 <- gsx_world(); expect_equal(gsx_run(w3), 0L); expect_equal(gsx_run(w3), 6L)
  w4 <- gsx_world(); expect_equal(gsx_run(w4, argv = w4$argv[1:4]), 2L)
})

test_that("실행기 새 그래프 모드 (실제 자료처럼 graph · reduction · 군집이 없음): 새 HVG/PCA/kNN/SNN/Louvain/DrL 을 점수 전에 만들고 '새'로 기록한다", {
  w <- gsx_world(); ob <- readRDS(w$inp); meta <- attr(ob, "meta.data"); meta$seurat_clusters <- NULL
  X <- gsx_counts(ob, "RNA")$counts; set.seed(3); bg <- matrix(stats::rpois(40 * ncol(X), rep(seq(0.5, 5, length.out = 40), ncol(X))), 40, ncol(X), dimnames = list(sprintf("BG%02d", 1:40), colnames(X)))
  X <- rbind(X, spm(bg)); ob2 <- gsx_fake(X, meta, v5 = FALSE); saveRDS(ob2, w$inp)                    # 배경 유전자 40 개: vst loess 추세가 성립하게
  pins <- w$pins; pins$input <- digest::digest(w$inp, algo = "sha256", file = TRUE); a <- w$argv
  for (k in c("--cluster-col", "--umap", "--graph")) a[match(k, a) + 1L] <- "NEW"
  a <- c(a, "--n-hvg", "8", "--n-pcs", "4", "--k", "5")
  expect_equal(gsx_run(w, pins = pins, argv = a), 0L); out <- w$out
  gp <- gsx_rd(out, "new_graph_parameters.csv"); expect_equal(gp$value[gp$name == "n_hvg"], "8"); expect_equal(gp$value[gp$name == "k"], "5")
  expect_equal(gp$value[gp$name == "seed"], "42"); expect_equal(gp$value[gp$name == "resolution"], "0.5")
  ck <- gsx_rd(out, "checks.csv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true("new_graph_built_before_scores" %in% ck$check)
  gc <- gsx_rd(out, "graph_checks.csv"); expect_equal(gc$value[gc$item == "graph_source"], "new_expression_snn (not original)")
  cd <- gsx_rd(out, "cluster_candidates.csv"); sc <- gsx_rd(out, "cluster_scores.csv"); expect_equal(sort(unique(sc$cluster)), sort(cd$cluster))
  expect_equal(sum(as.integer(cd$n_cells)), 24L); expect_true(any(grepl("새", readLines(file.path(out, "summary_ko.md")))))
  expect_true(all(file.exists(file.path(out, "figures", c("fig1_umap_author_vs_candidate.png", "fig2_granulocyte_zoom_heatmap.pdf")))))
})
