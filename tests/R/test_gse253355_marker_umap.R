# =============================================================================
# test_gse253355_marker_umap.R — 실제 R UMAP 에 고정 6 단계 RNA 패널 대응 (R74) 테스트: 유전자 기여 표(평균 / 목록 수, 없는 유전자 NA 와 0, 겹치는 유전자 진단),
#   기여 합 = 기존 군집 점수, 저자 구성(사후), 군집 안 1 위(winner) vs 단계별 최고(peak), mismatch 군집 기여 유전자, 그림 색 · 같은 좌표 · 번호 규칙, 실행기
#   계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_marker_umap_instruction.md (2026-10-19, 사용자 "R로 그려봐 · 우리 마커들이랑 대응시켜봐")
# =============================================================================
library(testthat)
gsm_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- gsm_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R",
            "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R",
            "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R", "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R",
            "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R", "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R",
            "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R", "73_gse253355_umap.R", "74_gse253355_marker_umap.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

spm <- function(m) methods::as(methods::as(Matrix::Matrix(m, sparse = TRUE), "generalMatrix"), "CsparseMatrix")
gsm_rd <- function(d, f) utils::read.csv(file.path(d, f), colClasses = "character", na.strings = character(0), check.names = FALSE)

# 소형 세계: 두 단계 x 목록 3 개. G1 은 두 단계에 겹치고, GX 는 자료에 없다. 군집 2 는 패널 유전자가 모두 0
gsm_toy <- function() {
  pan <- data.frame(stage = rep(c("Myeloblast", "PM"), each = 3), rank = rep(1:3, 2), symbol = c("G1", "G2", "GX", "G3", "G1", "G4"), stringsAsFactors = FALSE)
  genes <- c("G1", "G2", "G3", "G4", "OTHER"); cells <- sprintf("c%02d", 1:9); cl <- rep(c("0", "1", "2"), each = 3); dn <- rep(c("d1", "d2", "d3"), 3)
  X <- matrix(0, length(genes), 9, dimnames = list(genes, cells))
  X["G1", ] <- c(2, 0, 1, 5, 3, 0, 0, 0, 0); X["G2", ] <- c(1, 1, 0, 0, 0, 0, 0, 0, 0); X["G3", ] <- c(0, 0, 0, 2, 2, 2, 0, 0, 0); X["G4", ] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
  X["OTHER", ] <- c(40, 30, 50, 20, 60, 45, 70, 35, 55)
  pm <- gsx_panel_map(genes, pan); cs <- gsx_cell_scores(spm(X), pm, stages = c("Myeloblast", "PM"))
  clus <- gsx_cluster_scores(cs$fixed, cs$n_detected, cl, dn, 1L, A = cs$available)
  list(X = X, pm = pm, cs = cs, clus = clus, cl = cl)
}
gsm_toy_table <- function(t) gsm_gene_table(t$cs$gene_log, t$cs$gene_det, t$pm, t$cl, t$clus)

test_that("유전자 기여 표: 군집 x 단계 x 유전자 평균 log1p(CP10k) · 검출률 · 기여(평균 / 목록 수); 없는 유전자는 NA 와 0, 겹치는 유전자는 단계마다 한 행 + 진단", {
  t <- gsm_toy(); gt <- gsm_toy_table(t)
  expect_equal(nrow(gt), 3L * 6L); expect_equal(unique(gt$cluster), c("0", "1", "2")); expect_true(all(gt$n_listed_in_stage == 3L))
  lib <- colSums(t$X); w <- t$cl == "1"; m1 <- mean(log1p(1e4 * t$X["G1", w] / lib[w]))
  r <- gt[gt$cluster == "1" & gt$stage == "PM" & gt$gene == "G1", ]
  expect_equal(r$mean_log1p_cp10k, m1, tolerance = 1e-14); expect_equal(r$frac_cells_detected, 2 / 3); expect_equal(r$contribution, m1 / 3, tolerance = 1e-14)
  expect_equal(unique(gt$listed_in_stages[gt$gene == "G1"]), "Myeloblast|PM"); expect_equal(unique(gt$listed_in_stages[gt$gene == "G3"]), "PM")
  x <- gt[gt$gene == "GX", ]; expect_true(all(x$gene_status == "missing")); expect_true(all(is.na(x$mean_log1p_cp10k)) && all(is.na(x$frac_cells_detected)))
  expect_identical(x$contribution, c(0, 0, 0)); expect_true(all(gt$gene_status[gt$gene != "GX"] == "measured")); expect_false(anyNA(gt$mean_log1p_cp10k[gt$gene != "GX"]))
  s <- tapply(gt$contribution, paste(gt$cluster, gt$stage), sum); ref <- stats::setNames(t$clus$mean_all_cells, paste(t$clus$cluster, t$clus$stage))
  expect_equal(as.numeric(s[names(ref)]), unname(ref), tolerance = 1e-12)
  z <- gt[gt$cluster == "2", ]; expect_true(all(z$stage_score == 0)); expect_true(all(is.na(z$share_of_stage_score)))
  y <- gt[gt$cluster == "0" & gt$stage == "Myeloblast", ]; expect_equal(sum(y$share_of_stage_score), 1, tolerance = 1e-12)
  expect_equal(y$rank_in_stage[match(c("G1", "G2", "GX"), y$gene)], 1:3)                                 # 기여 큰 순서, 0 인 없는 유전자는 마지막
})

test_that("기여 합 검사: 유전자 기여 합 = 기존 군집 점수 (허용 1e-10); 어긋나면 그 행만 FALSE", {
  t <- gsm_toy(); gt <- gsm_toy_table(t); cc <- gsm_contribution_check(gt, t$clus, 1e-10)
  expect_equal(nrow(cc), 6L); expect_true(all(cc$pass)); expect_true(all(cc$abs_diff <= 1e-12))
  bad <- t$clus; bad$mean_all_cells[2] <- bad$mean_all_cells[2] + 1e-9; cb <- gsm_contribution_check(gt, bad, 1e-10)
  expect_equal(sum(!cb$pass), 1L); expect_equal(cb$abs_diff[!cb$pass], 1e-9, tolerance = 1e-5)
})

test_that("저자 구성 (사후): 군집마다 많은 순서와 비율 (R72 posthoc 표와 같은 정렬)", {
  d <- gsm_author_composition(c("0", "0", "0", "0", "1", "1", "1"), c("B", "A", "A", "A", "C", "B", "C"), c("0", "1"))
  expect_equal(d$cluster, c("0", "0", "1", "1")); expect_equal(d$author_label, c("A", "B", "C", "B")); expect_equal(d$n_cells, c(3L, 1L, 2L, 1L))
  expect_equal(d$frac_of_cluster, c(0.75, 0.25, 2 / 3, 1 / 3))
})

# 세 군집: 1 은 MY 의 단계별 최고(peak)지만 군집 안 1 위(winner)는 MM (실제 12 군집과 같은 꼴), 2 는 여러 단계의 peak
gsm_S <- function() matrix(c(0.5, 0.1, 0, 0, 0.2, 0, 0.1, 0.3, 0.42, 0.5, 0, 0, 0, 0.6, 0.1, 0.9, 0.3, 0.05), 3, byrow = TRUE, dimnames = list(c("0", "1", "2"), NSO_STAGES))
gsm_cand <- function(S, n = c(10L, 20L, 30L)) {
  rk <- atlas_cl_rank(S, 1e-12); d <- data.frame(cluster = rownames(S), n_cells = n, n_donors = c(2L, 3L, 3L), stringsAsFactors = FALSE)
  for (s in NSO_STAGES) d[[paste0("score_", s)]] <- unname(S[, s])
  d <- cbind(d, rk[, setdiff(names(rk), "cluster")]); d$donor_equal_top_candidate <- c("Myeloblast", "MY", "MM"); d$donor_equal_agreement <- c("same_top", "different_top", "same_top"); d
}

test_that("군집 요약: 군집 안 1 위(winner) · 2 위 · 차이 · donor 동일 가중 일치와 단계별 전체 군집 최고(peak)를 따로 적는다", {
  S <- gsm_S(); cand <- gsm_cand(S); peaks <- atlas_cl_stage_peaks(S, atlas_cl_rank(S, 1e-12), 1e-12, stats::setNames(cand$n_cells, cand$cluster))
  comp <- gsm_author_composition(rep(c("0", "1", "2"), c(10, 20, 30)), c(rep("HSC", 7), rep("GMP", 3), rep("Late Myeloid", 12), rep("Early Myeloid Progenitor", 8), rep("Neutrophil", 30)), c("0", "1", "2"))
  sm <- gsm_cluster_summary(cand, peaks, comp, c(`0` = 10L, `1` = 0L, `2` = 30L))
  expect_equal(sm$cluster, c("0", "1", "2")); expect_equal(sm$winner_stage, c("Myeloblast", "MM", "MM")); expect_equal(sm$second_stage, c("BN", "MY", "PM"))
  expect_equal(sm$winner_minus_second, c(0.3, 0.08, 0.3), tolerance = 1e-12)
  expect_equal(sm$stage_peak_of, c("Myeloblast", "MY", "PM|MM|BN|SN")); expect_equal(sm$winner_is_stage_peak, c(TRUE, FALSE, TRUE))
  expect_equal(sm$author_top_label, c("HSC", "Late Myeloid", "Neutrophil")); expect_equal(sm$author_top_frac, c(0.7, 0.6, 1))
  expect_equal(sm$author_second_label, c("GMP", "Early Myeloid Progenitor", "")); expect_true(is.na(sm$author_second_frac[3]))
  expect_equal(sm$author_composition[1], "HSC 7 (70%) | GMP 3 (30%)"); expect_equal(sm$n_subset_cells, c(10L, 0L, 30L))
  expect_equal(sm$donor_equal_agreement, c("same_top", "different_top", "same_top")); expect_equal(sm$score_MM, unname(S[, "MM"]))
  expect_false(anyNA(sm[vapply(sm, is.character, TRUE)]))
})

test_that("mismatch 군집 기여 유전자: 1 위 · 2 위 · peak 단계마다 기여 > 0 상위 n 개와 누적 비율; 없는 군집은 멈춘다", {
  t <- gsm_toy(); gt <- gsm_toy_table(t)
  sm <- data.frame(cluster = c("0", "1", "2"), winner_stage = c("Myeloblast", "PM", ""), second_stage = c("PM", "Myeloblast", ""), stage_peak_of = c("Myeloblast", "PM", ""),
                   author_top_label = c("HSC", "Neutrophil", "T"), author_top_frac = c(1, 1, 1), stringsAsFactors = FALSE)
  fg <- gsm_focus_genes(gt, sm, c("0", "1"), n_top = 2L)
  r0 <- fg[fg$cluster == "0", ]; expect_equal(unique(r0$stage), c("Myeloblast", "PM"))
  expect_equal(unique(r0$role[r0$stage == "Myeloblast"]), "winner+peak"); expect_equal(unique(r0$role[r0$stage == "PM"]), "second")
  expect_equal(r0$gene[r0$stage == "Myeloblast"], c("G1", "G2")); expect_equal(r0$gene[r0$stage == "PM"], "G1")    # 0 인 유전자(GX · G3 · G4)는 적지 않는다
  expect_equal(r0$cumulative_share[r0$stage == "Myeloblast"][2], 1, tolerance = 1e-12); expect_equal(unique(r0$listed_in_stages[r0$gene == "G1"]), "Myeloblast|PM")
  expect_true(all(fg$contribution > 0)); expect_lte(max(table(paste(fg$cluster, fg$stage))), 2L)
  expect_equal(unique(fg$role[fg$cluster == "1" & fg$stage == "PM"]), "winner+peak")
  expect_error(gsm_focus_genes(gt, sm, c("0", "9")), "9")
})

test_that("그림: 오른쪽 색 = CSV 1 위 단계 색 (군집마다 한 색), 부분집합 = 같은 좌표의 모든 선택 세포, 번호 = 부분집합 세포 >= 비율 기준 군집, heatmap 의 winner · peak 표식 구분", {
  set.seed(3); n <- 90; cl <- rep(c("0", "1", "2"), c(40, 40, 10)); U <- cbind(stats::rnorm(n) + rep(c(0, 5, 10), c(40, 40, 10)), stats::rnorm(n)); rownames(U) <- sprintf("c%03d", 1:n)
  au <- rep(c("HSC", "Neutrophil", "Plasma Cell"), c(40, 40, 10)); au[81] <- "HSC"                         # 군집 2 에 부분집합 세포 1 개
  S <- gsm_S(); cand <- gsm_cand(S, c(40L, 40L, 10L)); peaks <- atlas_cl_stage_peaks(S, atlas_cl_rank(S, 1e-12), 1e-12, stats::setNames(cand$n_cells, cand$cluster))
  sl <- c("HSC", "Neutrophil"); acols <- .gsu_cols(au, sl); mask <- au %in% sl
  p1 <- tempfile(fileext = ".png"); r1 <- gsm_draw_pair(p1, U, au, acols, cl, cand, ko = FALSE, fam = "")
  expect_true(file.exists(p1)); expect_identical(r1$coords, U); expect_equal(r1$numbered, c("0", "1", "2"))
  p2 <- tempfile(fileext = ".png"); r2 <- gsm_draw_pair(p2, U, au, acols, cl, cand, mask = mask, subset_labels = sl, label_min_frac = 0.05, ko = FALSE, fam = "")
  expect_true(file.exists(p2)); expect_identical(r2$coords, U[mask, , drop = FALSE]); expect_equal(nrow(r2$coords), 81L)
  expect_equal(r2$numbered, c("0", "1")); expect_equal(r2$unnumbered$cluster, "2"); expect_equal(r2$unnumbered$n_cells, 1L)
  for (r in list(r1, r2)) { expect_true(all(r$colours$n_colours == 1L)); expect_true(all(r$colours$colour == ATLAS_CL_STAGE_COL[r$colours$winner]))
    expect_false(any(r$text$text %in% c("NA", "NaN", "Inf"))); expect_false(any(grepl("→|->", r$text$text))) }
  summ <- data.frame(cluster = c("0", "1", "2"), author_top_label = c("HSC", "Neutrophil", "Plasma Cell"), author_top_frac = c(1, 1, 0.9), n_cells = c(40L, 40L, 10L),
                     donor_equal_agreement = c("same_top", "different_top", "same_top"), stringsAsFactors = FALSE)
  h <- gsm_draw_heatmap(tempfile(fileext = ".png"), S, cand, peaks, summ, ko = FALSE, fam = "")
  w <- h$marks[h$marks$mark == "winner", ]; expect_equal(paste(w$cluster, w$stage), c("0 Myeloblast", "1 MM", "2 MM"))
  pk <- h$marks[h$marks$mark == "peak", ]; expect_equal(sort(paste(pk$cluster, pk$stage)), sort(c("0 Myeloblast", "2 PM", "1 MY", "2 MM", "2 BN", "2 SN")))
  expect_false(any(h$text$text %in% c("NA", "NaN", "Inf")))
})


# ── 실행기: gse253355 형 소형 결과 (R72 실행기를 가짜 Seurat 에 돌린 것) + 가짜 UMAP 좌표 RDS ────────────────
#   가짜 Seurat 클래스 · 세계는 test_gse253355_application.R 의 것과 같다 (테스트 도구, 실행 코드 아님)
suppressMessages({
  setClass("LogMap", contains = "matrix")
  setClass("Assay5", representation(layers = "list", cells = "LogMap", features = "LogMap", key = "character"))
  setClass("Graph", contains = "dgCMatrix", representation(assay.used = "character"))
  setClass("DimReduc", representation(cell.embeddings = "matrix", key = "character"))
  setClass("Seurat", representation(assays = "list", meta.data = "data.frame", active.assay = "character", graphs = "list", reductions = "list", project.name = "character"))
})
gsm_fake_seurat <- function(counts, meta, graph, umap, split) {
  genes <- rownames(counts); cells <- colnames(counts); unm <- function(m) { dimnames(m) <- list(NULL, NULL); m }
  grp <- paste0("counts.", split); ln <- unique(grp); lay <- stats::setNames(lapply(ln, function(l) unm(counts[, grp == l, drop = FALSE])), ln)
  fm <- matrix(TRUE, length(genes), length(ln), dimnames = list(genes, ln)); cm <- vapply(ln, function(l) grp == l, logical(length(cells))); dim(cm) <- c(length(cells), length(ln)); dimnames(cm) <- list(cells, ln)
  a <- new("Assay5", layers = lay, cells = new("LogMap", cm), features = new("LogMap", fm), key = "rna_")
  new("Seurat", assays = list(RNA = a), meta.data = meta, active.assay = "RNA", graphs = list(RNA_snn = new("Graph", graph, assay.used = "RNA")),
      reductions = list(UMAP_dim30 = new("DimReduc", cell.embeddings = umap, key = "UMAPdim30_")), project.name = "t")
}
gsm_world <- function() {
  w <- tempfile("gsm_world"); dir.create(w); sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE); sd <- file.path(w, "out", "gse253355_02")
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
  G <- Matrix::sparseMatrix(i = c(1:5, 6, 7:11, 12, 13:17, 19:23), j = c(2:6, 7, 8:12, 13, 14:18, 20:24), x = 1, dims = c(n, n), symmetric = TRUE)
  G <- methods::as(methods::as(G, "generalMatrix"), "CsparseMatrix"); dimnames(G) <- list(cells, cells)
  U0 <- cbind(UMAPdim30_1 = rep(1:4, each = 6) * 5 + seq_len(n) %% 3, UMAPdim30_2 = seq_len(n) %% 4); rownames(U0) <- cells
  inp <- file.path(w, "GSE253355_test.rds.gz"); saveRDS(gsm_fake_seurat(spm(X), meta, G, U0, dn), inp)
  hg <- file.path(w, "hgnc.txt"); utils::write.table(data.frame(symbol = "S6b", prev_symbol = "", alias_symbol = "S6", status = "Approved"), hg, sep = "\t", row.names = FALSE, quote = FALSE)
  a72 <- c("--input", inp, "--panels", pf, "--panel-attempt", att, "--hgnc", hg, "--out-dir", sd, "--public-dir", file.path(w, "pub72"), "--assay", "RNA", "--cluster-col", "seurat_clusters",
           "--donor-col", "orig.ident", "--author-col", "cluster_anno_l2", "--umap", "UMAP_dim30", "--graph", "RNA_snn", "--granulocyte-labels", "GMP;Neutrophil", "--neutrophil-labels", "Neutrophil")
  code <- NA_integer_; utils::capture.output(code <- gsx_main(a72, root = ROOT, pins = list(input = sha(inp), source_freeze = sha(file.path(att, "source_freeze.json")))))
  if (!identical(code, 0L)) stop("테스트 세계: R72 실행 실패")
  ls <- readRDS(file.path(sd, "lossless.rds")); sl <- c("GMP", "Neutrophil", "Monocyte"); U <- cbind(umap_1 = seq_len(n) %% 5 + rep(0:3, each = 6) * 4, umap_2 = seq_len(n) %/% 5)
  rownames(U) <- ls$cells; up <- file.path(w, "umap_coordinates.rds")
  saveRDS(list(cells = ls$cells, U = U, new_cluster = ls$clusters, author_label = au, subset_labels = sl, subset_mask = au %in% sl), up)
  ac <- .gsu_cols(au, sl); atlas_write_csv(data.frame(author_label = names(ac), n_cells = as.integer(table(au)[names(ac)]), in_subset = names(ac) %in% sl, colour = unname(ac),
                                                      stringsAsFactors = FALSE), file.path(w, "annotation_counts.csv"))   # R73 색 표 (좌표 RDS 옆)
  out <- file.path(w, "out", "gse253355_marker_umap_01"); pub <- file.path(w, "pub", "gse253355_marker_umap_01")
  argv <- c("--umap", up, "--scores-dir", sd, "--object", inp, "--panels", pf, "--hgnc", hg, "--out-dir", out, "--public-dir", pub, "--focus-clusters", "0;1;2")
  list(argv = argv, out = out, pub = pub, sd = sd, n = n,
       pins = list(umap = sha(up), lossless = sha(file.path(sd, "lossless.rds")), object = sha(inp), calc_manifest = sha(file.path(sd, "calc_manifest.json"))))
}
gsm_run <- function(w, pins = w$pins, argv = w$argv) { code <- NA_integer_; utils::capture.output(code <- gsm_main(argv, root = ROOT, pins = pins)); code }

test_that("실행기: 고정 입력 대조 → 유전자 기여 표 · 합 검사 (계산 고정) → 사후 저자 구성 · 군집 요약 · mismatch 유전자 · 그림 3 장 · 공개 사본 (세포 단위 없음)", {
  w <- gsm_world(); expect_equal(gsm_run(w), 0L); out <- w$out
  need <- c("cluster_stage_gene_contributions.csv", "contribution_check.csv", "panel_gene_diagnostics.csv", "calc_manifest.json", "posthoc_author_composition.csv",
            "cluster_marker_summary.csv", "focus_cluster_top_genes.csv", "figure_text.csv", "checks.csv", "environment.csv", "input_sha256.csv", "source_snapshot.csv",
            "manifest.json", "summary_ko.md", file.path("figures", c("fig1_umap_full_author_vs_marker.png", "fig2_umap_lineage_author_vs_marker.png", "fig3_cluster_stage_heatmap.png")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  ck <- gsm_rd(out, "checks.csv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true(all(c("inputs_pinned", "cells_clusters_coordinates_unchanged", "frozen_scores_equal_lossless", "gene_values_reproduce_frozen_scores", "contribution_sum_equals_cluster_score",
                    "missing_gene_na_and_zero", "author_composition_matches_existing", "author_colours_match_r73", "figure_colours_match_winner", "subset_same_coordinates", "existing_inputs_unchanged",
                    "tables_roundtrip_ok", "public_copy_aggregate_only") %in% ck$check))
  gt <- gsm_rd(out, "cluster_stage_gene_contributions.csv"); r <- gt[gt$cluster == "1" & gt$stage == "SN" & gt$gene == "S1", ]
  expect_equal(as.numeric(r$mean_log1p_cp10k), log1p(1e4 * 4 / 56), tolerance = 1e-12); expect_equal(as.numeric(r$contribution), log1p(1e4 * 4 / 56) / 2, tolerance = 1e-12)
  m <- gt[gt$gene == "S6b", ]; expect_equal(nrow(m), 4L); expect_true(all(m$mean_log1p_cp10k == "NA")); expect_true(all(as.numeric(m$contribution) == 0))
  sm <- gsm_rd(out, "cluster_marker_summary.csv"); expect_equal(sm$cluster, c("0", "1", "2", "3")); expect_equal(sm$winner_stage[sm$cluster == "0"], "Myeloblast")
  expect_equal(sm$stage_peak_of[sm$cluster == "2"], "SN"); expect_equal(sm$stage_peak_of[sm$cluster == "1"], "BN"); expect_equal(sm$winner_status[sm$cluster == "3"], "unassigned_zero")
  expect_equal(sm$author_top_label, c("GMP", "Neutrophil", "Monocyte", "CD4+ T-Cell")); expect_equal(sm$n_subset_cells, c("6", "6", "6", "0"))
  fg <- gsm_rd(out, "focus_cluster_top_genes.csv"); expect_true(all(fg$cluster %in% c("0", "1", "2"))); expect_equal(fg$gene[fg$cluster == "2"], "S1")
  pubf <- list.files(w$pub, recursive = TRUE)
  expect_true(all(c("summary_ko.md", "cluster_marker_summary.csv", "cluster_stage_gene_contributions.csv", "focus_cluster_top_genes.csv", "fig3_cluster_stage_heatmap.png") %in% basename(pubf)))
  expect_false(any(grepl("lossless|coordinates|\\.rds$|\\.gz$", pubf)))                                    # 공개 사본: 세포 단위 · 좌표 없음
  expect_lte(length(readLines(file.path(out, "summary_ko.md"))), 25L)
  ft <- gsm_rd(out, "figure_text.csv"); expect_false(any(ft$text %in% c("NA", "NaN", "Inf"))); expect_false(any(grepl("→|->", ft$text)))
  cmf <- names(jsonlite::fromJSON(file.path(out, "calc_manifest.json"))$files); expect_true("cluster_stage_gene_contributions.csv" %in% cmf); expect_false(any(grepl("posthoc|figure|summary", cmf)))
  man <- jsonlite::fromJSON(file.path(out, "manifest.json"))$files
  expect_true(all(vapply(names(man), function(f) identical(digest::digest(file.path(out, f), algo = "sha256", file = TRUE), man[[f]]), TRUE)))
  expect_equal(gsm_run(w), 6L)                                                                             # 덮어쓰지 않는다
  w2 <- gsm_world(); bad <- w2$pins; bad$umap <- strrep("0", 64); expect_equal(gsm_run(w2, pins = bad), 4L); expect_true(any(grepl("hash", readLines(file.path(w2$out, "run.log")))))
})
