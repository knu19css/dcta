# =============================================================================
# test_gse253355_umap.R — GSE253355 새 UMAP 표시 (R73) 테스트: 저장 PCA · kNN 재사용, 이웃 거리 정확 계산(자기 0), uwot 인자, 입력 검사, 부분집합 = 같은 좌표
#   계약: codex 전달 (2026-10-19, 사용자 새 요청: 전체 UMAP 1 개 + HSC→Neutrophil/Monocyte 계통 UMAP 1 개). uwot 은 설치 승인 뒤에만 쓴다 — 테스트는 기록용 가짜 함수를 넣는다.
# =============================================================================
library(testthat)
gsu_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- gsu_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R",
            "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R",
            "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R", "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R",
            "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R", "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R",
            "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R", "73_gse253355_umap.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

fake_umap <- function(calls) function(X, n_neighbors, min_dist, spread, n_components, n_epochs, metric, scale, init, nn_method, n_threads, n_sgd_threads, fast_sgd, seed, verbose = FALSE) {
  calls$args <- list(X = X, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, n_components = n_components, n_epochs = n_epochs, metric = metric, scale = scale,
                     init = init, nn_method = nn_method, n_threads = n_threads, n_sgd_threads = n_sgd_threads, fast_sgd = fast_sgd, seed = seed)
  cbind(X[, 1], X[, 2] + 0.01 * stats::runif(nrow(X)))                                                 # 결정적 가짜 좌표 (set.seed 뒤 난수)
}

test_that("이웃 거리: 저장 kNN 번호마다 저장 PCA 로 정확한 유클리드 거리, 자기 자신은 첫 열 0", {
  set.seed(5); Y <- matrix(stats::rnorm(60 * 4), 60, 4); idx <- gsx_knn(Y, k = 6L)
  d <- gsu_knn_dist(Y, idx); D <- as.matrix(stats::dist(Y))
  expect_equal(dim(d), c(60L, 6L)); expect_true(all(d[, 1] == 0)); expect_equal(d, matrix(D[cbind(rep(1:60, 6), as.vector(idx))], 60, 6), tolerance = 1e-12)
  expect_true(all(apply(d, 1, function(z) all(diff(z) >= -1e-12))))
})

test_that("uwot 호출: 입력은 저장 PCA 와 미리 계산한 idx / dist 뿐이고 고정 파라미터를 그대로 넘긴다; 필수 인자가 API 에 없거나 uwot 이 없으면 멈춘다 (DrL 대체 없음)", {
  set.seed(5); Y <- matrix(stats::rnorm(40 * 3), 40, 3, dimnames = list(sprintf("c%02d", 1:40), NULL)); idx <- gsx_knn(Y, k = 5L); d <- gsu_knn_dist(Y, idx)
  calls <- new.env(); r <- gsu_umap_call(Y, idx, d, GSU_UMAP_ARGS, 42L, umap_fun = fake_umap(calls))
  a <- calls$args; expect_identical(a$X, Y); expect_identical(a$nn_method$idx, idx); expect_identical(a$nn_method$dist, d)
  expect_equal(c(a$n_neighbors, a$min_dist, a$spread, a$n_components, a$n_epochs, a$n_threads, a$n_sgd_threads, a$seed), c(20, 0.3, 1, 2, 200, 2, 1, 42))
  expect_equal(c(a$metric, a$init), c("euclidean", "spectral")); expect_false(a$scale); expect_false(a$fast_sgd)
  expect_equal(dim(r$U), c(40L, 2L)); expect_equal(rownames(r$U), rownames(Y)); expect_identical(r$U, gsu_umap_call(Y, idx, d, GSU_UMAP_ARGS, 42L, umap_fun = fake_umap(new.env()))$U)
  expect_error(gsu_umap_call(Y, idx, d, GSU_UMAP_ARGS, 42L, umap_fun = function(X, nn_method, n_neighbors) X[, 1:2]), "API")
  if (!nzchar(system.file(package = "uwot"))) expect_error(gsu_umap_call(Y, idx, d, GSU_UMAP_ARGS, 42L), "uwot")
})

test_that("색: 계통 8 라벨은 고정된 뚜렷한 색, 35 라벨 전체에 같은 색이 두 번 쓰이지 않고 두 그림이 같은 매핑을 쓴다", {
  lab <- c(GSU_SUBSET, sprintf("Other%02d", 1:27)); cols <- .gsu_cols(lab, GSU_SUBSET)
  expect_equal(unname(cols[GSU_SUBSET]), unname(GSU_LINEAGE_COLS[GSU_SUBSET])); expect_false(anyDuplicated(unname(cols)) > 0); expect_equal(length(cols), 35L)
  expect_identical(cols, .gsu_cols(rev(lab), GSU_SUBSET))                                                # 입력 순서와 무관한 같은 매핑
})

test_that("부분집합 라벨: 정확한 저자 라벨만 고르고 수를 센다; 없는 라벨은 멈춘다", {
  lab <- c(rep("HSC", 3), rep("Neutrophil", 2), "Plasma Cell", "Monocyte")
  s <- gsu_subset(lab, c("HSC", "Neutrophil", "Monocyte")); expect_equal(sum(s$mask), 6L); expect_equal(s$counts$n_cells, c(3L, 2L, 1L))
  expect_error(gsu_subset(lab, c("HSC", "Neutrophils")), "Neutrophils")
})


# ── 실행기 (가짜 lossless · 가짜 객체) ─────────────────────────────────────────
gsu_world <- function() {
  w <- tempfile("gsu_world"); dir.create(w); sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  set.seed(7); n <- 60; cells <- sprintf("H%02d_c%03d", rep(1:2, each = 30), 1:n); Y <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(cells, paste0("PC", 1:5)))
  ls <- list(cells = cells, clusters = as.character(rep(0:2, each = 20)), graph_build = list(pca_scores = Y, knn = gsx_knn(Y, k = 8L)))
  lp <- file.path(w, "lossless.rds"); saveRDS(ls, lp)
  lab <- rep(c("HSC", "GMP", "Neutrophil", "Monocyte", "Plasma Cell", "CD4+ T-Cell"), each = 10)
  setClass("SeuratU", representation(meta.data = "data.frame")); ob <- new("SeuratU", meta.data = data.frame(cluster_anno_l2 = factor(lab), row.names = cells))
  op <- file.path(w, "object.rds"); saveRDS(ob, op)
  list(w = w, argv = c("--lossless", lp, "--object", op, "--out-dir", file.path(w, "out", "gse253355_umap_01"), "--public-dir", file.path(w, "pub", "gse253355_umap_01"),
                       "--subset-labels", "HSC;GMP;Neutrophil;Monocyte"),
       pins = list(lossless = sha(lp), object = sha(op)), out = file.path(w, "out", "gse253355_umap_01"), pub = file.path(w, "pub", "gse253355_umap_01"), n = n)
}
gsu_rd <- function(d, f) utils::read.csv(file.path(d, f), colClasses = "character", na.strings = character(0), check.names = FALSE)

test_that("실행기: 사전 검사 모드는 UMAP 없이 입력만 확인하고, 계산 모드는 같은 좌표로 두 그림 · 좌표 표 · 검사를 남긴다", {
  w <- gsu_world(); pre <- c(w$argv, "--preflight", "TRUE"); pre[match("--out-dir", pre) + 1L] <- file.path(w$w, "pre")
  expect_equal(gsu_main(pre, root = ROOT, pins = w$pins), 0L)
  pc <- gsu_rd(file.path(w$w, "pre"), "preflight_checks.csv"); expect_true(all(pc$pass == "TRUE")); expect_false(file.exists(file.path(w$w, "pre", "figures")))
  calls <- new.env(); expect_equal(gsu_main(w$argv, root = ROOT, pins = w$pins, umap_fun = fake_umap(calls)), 0L); out <- w$out
  need <- c("umap_coordinates.rds", "umap_coordinates_local.csv.gz", "annotation_counts.csv", "checks.csv", "environment.csv", "manifest.json", "summary_ko.md", "figure_text.csv",
            file.path("figures", c("fig_umap_full_author.png", "fig_umap_hsc_neutrophil_monocyte.png")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  ck <- gsu_rd(out, "checks.csv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true(all(c("all_cells_kept", "subset_coordinates_identical", "umap_args_recorded") %in% ck$check))
  ac <- gsu_rd(out, "annotation_counts.csv"); expect_equal(ac$n_cells[ac$author_label == "Neutrophil"], "10"); expect_equal(sum(ac$in_subset == "TRUE"), 4L)
  co <- utils::read.csv(gzfile(file.path(out, "umap_coordinates_local.csv.gz")), stringsAsFactors = FALSE); expect_equal(nrow(co), w$n); expect_true(all(c("cell_id", "umap_1", "umap_2", "new_cluster", "author_label") %in% names(co)))
  pubf <- list.files(w$pub, recursive = TRUE); expect_true(all(c("summary_ko.md", "fig_umap_full_author.png", "fig_umap_hsc_neutrophil_monocyte.png") %in% basename(pubf)))
  expect_false(any(grepl("coordinates", pubf)))                                                          # 공개 사본에 세포 좌표 없음
  expect_lte(length(readLines(file.path(out, "summary_ko.md"))), 25L); ft <- gsu_rd(out, "figure_text.csv"); expect_false(any(ft$text %in% c("NA", "NaN", "Inf")))
  expect_false(any(grepl("→|->", ft$text)))                                                              # 그림 글자에도 방향 화살표가 없다
  expect_equal(gsu_main(w$argv, root = ROOT, pins = w$pins, umap_fun = fake_umap(new.env())), 6L)
  w2 <- gsu_world(); bad <- w2$pins; bad$lossless <- strrep("0", 64); expect_equal(gsu_main(w2$argv, root = ROOT, pins = bad, umap_fun = fake_umap(new.env())), 4L)
})
