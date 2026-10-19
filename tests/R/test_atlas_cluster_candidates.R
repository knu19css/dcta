# =============================================================================
# test_atlas_cluster_candidates.R — atlas 6 단계 패널의 군집 평균 점수 · 군집 내부 1 위 후보 · 단계별 최고 군집 (R66–67) 테스트
#   계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_cluster_instruction.md
#   세포 수 가중(pooled) vs 동일 block 가중 차이(손계산), 모두 0 · 결측 · 동률, 작은 군집 보존, 단계 최고 ≠ 군집 내 최고, 세포 ID 뒤섞임 중지,
#   실제 attempt_02 cluster_programs 와 숫자 일치, 덮어쓰기 거절 · 고정 hash 불일치 중지, 계산 고정 전 UMAP 읽기 금지
# =============================================================================
library(testthat)
acl_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- acl_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R",
            "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R",
            "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R", "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R",
            "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R", "67_atlas_cluster_view.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

ST <- NSO_STAGES

test_that("순위: 정상 1 위 · 2 위 · 차이, 모두 0 = unassigned_zero, 모두 결측 = unavailable, 허용오차 안 동률 = tied (후보 모두), 일부 결측 = incomplete_panel_comparison", {
  S <- rbind(c1 = c(0.5, 0.2, 0.1, 0, 0, 0), c2 = rep(0, 6), c3 = rep(NA_real_, 6), c4 = c(0.4, 0.4, NA, 0, 0, 0),
             c5 = c(0.1, NA, NA, NA, NA, 0.3), c6 = c(0.3, 0.3 + 1e-14, 0, 0, 0, 0), c7 = c(0.3, 0.3 + 1e-9, 0, 0, 0, 0))
  colnames(S) <- ST
  r <- atlas_cl_rank(S, tol = 1e-12)
  expect_equal(r$score_top_status, c("score_top_candidate", "unassigned_zero", "unavailable", "tied", "score_top_candidate", "tied", "score_top_candidate"))
  expect_equal(r$score_top_candidate, c("Myeloblast", "", "", "Myeloblast|PM", "SN", "Myeloblast|PM", "PM"))
  expect_equal(r$score_second_candidate[1], "PM"); expect_equal(r$score_top_minus_second[1], 0.3); expect_equal(r$n_top_tied[c(1, 4, 6)], c(1L, 2L, 2L))
  expect_equal(r$score_top_minus_second[4], 0); expect_equal(r$score_top_value[5], 0.3); expect_equal(r$score_second_value[5], 0.1)
  expect_equal(r$comparison_scope, c("complete", "complete", "none", "incomplete_panel_comparison", "incomplete_panel_comparison", "complete", "complete"))
  expect_true(is.na(r$score_top_value[3])); expect_true(is.na(r$score_top_minus_second[3]))
  one <- atlas_cl_rank(rbind(a = c(0.2, NA, NA, NA, NA, NA)), tol = 1e-12)                       # 비교할 2 위가 없으면 2 위 · 차이는 결측
  expect_equal(one$score_top_candidate, "Myeloblast"); expect_true(is.na(one$score_second_value)); expect_true(is.na(one$score_top_minus_second))
})

test_that("0 은 실제 0 이다 (Codex 검토): 1e-14 같은 작은 양수는 unassigned_zero / all_zero 가 아니고, 점수 0 단계는 작은 양수와 동률 후보로 묶지 않는다", {
  S <- rbind(z = c(1e-14, 0, 0, 0, 0, 0), y = c(0, 0, 0, 0, 0, 0)); colnames(S) <- ST
  r <- atlas_cl_rank(S, tol = 1e-12)
  expect_equal(r$score_top_status, c("score_top_candidate", "unassigned_zero")); expect_equal(r$score_top_candidate, c("Myeloblast", "")); expect_equal(r$n_top_tied, c(1L, 0L))
  pk <- atlas_cl_stage_peaks(S, r, tol = 1e-12)
  expect_equal(pk$peak_status[pk$stage == "Myeloblast"], "peak"); expect_equal(pk$peak_clusters[pk$stage == "Myeloblast"], "z"); expect_equal(pk$peak_status[pk$stage == "PM"], "all_zero")
})

test_that("단계별 최고 군집은 군집 내부 1 위와 다른 비교다: 동률 모두 유지, 어떤 군집의 1 위도 아닌 단계에도 최고 군집이 있다", {
  S <- rbind(A = c(0.9, 0.5, 0, 0, 0, 0.05), B = c(0.2, 0.4, 0, 0, 0, 0.05), C = c(0, 0, 0.3, 0.3, 0, 0)); colnames(S) <- ST
  rk <- atlas_cl_rank(S, tol = 1e-12); pk <- atlas_cl_stage_peaks(S, rk, tol = 1e-12)
  expect_equal(pk$stage, ST)
  pm <- pk[pk$stage == "PM", ]; expect_equal(pm$peak_clusters, "A"); expect_equal(pm$peak_cluster_top_candidate, "Myeloblast"); expect_false(pm$peak_is_cluster_top)
  expect_equal(pm$second_cluster, "B"); expect_equal(pm$peak_minus_second, 0.1)
  sn <- pk[pk$stage == "SN", ]; expect_equal(sn$peak_clusters, "A|B"); expect_equal(sn$n_tied, 2L); expect_equal(sn$n_clusters_with_stage_top, 0L)
  expect_equal(pk$peak_status[pk$stage == "BN"], "all_zero"); expect_equal(pk$peak_clusters[pk$stage == "BN"], "")
  expect_equal(pk$n_clusters_with_stage_top[pk$stage %in% c("MY", "MM")], c(1L, 1L))                 # C 는 MY · MM 동률 → 두 단계 모두 1 위 후보
})

test_that("군집 번호 배치: 자리가 있으면 겹치지 않게 놓고, 겹치지 않는 번호는 대표 위치 중심에 둔다", {
  lp <- atlas_cl_place_labels(x = c(0, 0.05, 5), y = c(0, 0, 5), w = rep(1, 3), h = rep(0.5, 3), xlim = c(-10, 10), ylim = c(-10, 10))
  ov <- function(i, j) !(lp$x1[i] <= lp$x0[j] || lp$x1[j] <= lp$x0[i] || lp$y1[i] <= lp$y0[j] || lp$y1[j] <= lp$y0[i])
  expect_false(ov(1, 2)); expect_false(ov(1, 3)); expect_false(ov(2, 3)); expect_true(all(lp$free))
  expect_equal(c(lp$lx[1], lp$ly[1], lp$lx[3], lp$ly[3]), c(0, 0, 5, 5)); expect_equal(lp$moved, c(FALSE, TRUE, FALSE))
})

test_that("세포 대표 위치: 군집 세포 중 UMAP 중앙값에 가장 가까운 실제 세포 (표시 전용)", {
  U <- cbind(c(0, 1, 2, 10, 11, 30), c(0, 0, 0, 5, 5, 0)); rc <- atlas_cl_rep_cells(U, c("a", "a", "a", "b", "b", "c"))
  expect_equal(rc$cell[rc$cluster == "a"], 2L); expect_true(rc$cell[rc$cluster == "b"] %in% c(4L, 5L)); expect_equal(rc$cell[rc$cluster == "c"], 6L)
})


# ── 소형 합성 attempt (R59 를 그대로 거친 lossless · cluster_programs) ────────────
acl_world <- function(shuffle_ids = FALSE) {
  w <- tempfile("acl_world"); dir.create(w); att <- file.path(w, "attempt"); dir.create(att); out <- file.path(w, "out", "cluster_view_01")
  sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  # 군집: A 40 (block b1 30 = S1 5 count, b2 10 = S2 7 count), B 3 (S1, 작은 군집), C 6 (패널 0), D 6 (S3 3 개 · S4 3 개, 동률), E 6 (S5 5 · 둘째 세포만 S2 1 · 첫 세포만 S6 1)
  cl <- c(rep("A", 40), rep("B", 3), rep("C", 6), rep("D", 6), rep("E", 6)); n <- length(cl)
  blk <- c(rep("b1", 30), rep("b2", 10), rep("b1", 3), rep("b1", 6), rep("b2", 6), rep("b1", 6))
  g <- c("S1", "S2", "S3", "S4", "S5", "S6"); X <- matrix(0, n, 7, dimnames = list(NULL, c(g, "OTHER")))
  X[1:30, "S1"] <- 5; X[31:40, "S2"] <- 7; X[41:43, "S1"] <- 5; X[50:52, "S3"] <- 5; X[53:55, "S4"] <- 5; X[56:61, "S5"] <- 5; X[57, "S2"] <- 1; X[56, "S6"] <- 1
  X[, "OTHER"] <- 100 - rowSums(X[, g]); X <- Matrix::Matrix(X, sparse = TRUE)
  ids <- sprintf("cell%02d", 1:n); md <- data.frame(cell_id = ids, cluster92 = cl, cluster51 = ifelse(cl %in% c("A", "B"), "X", "Y"), donor = blk, study = "s1", block = blk,
                                                     stringsAsFactors = FALSE)
  ref <- list(panels = data.frame(stage = ST, rank = 1L, ensg = sprintf("ENSG%011d", 1:6), symbol = g, stringsAsFactors = FALSE),
              transitions = data.frame(pair = character(0), direction = character(0), rank = integer(0), ensg = character(0), symbol = character(0), stringsAsFactors = FALSE))
  nref <- atlas_nso_reference(ref)
  groups <- list(p92 = data.frame(cluster = md$cluster92, block = md$block, stringsAsFactors = FALSE), p51 = data.frame(cluster = md$cluster51, block = md$block, stringsAsFactors = FALSE))
  ge <- nso_gene_evidence(nrv_chunks_matrix(X, 25L), nrv_gene_map(colnames(X), g), groups); pe <- nso_program_evidence(ge, nref)
  saveRDS(list(metadata = md, programs = pe, nso_reference = nref), file.path(att, "lossless.rds"))
  mbvt_write_tsv(pe$clusters, file.path(att, "cluster_programs.tsv"))
  cc <- do.call(rbind, lapply(c("p92", "p51"), function(p) { k <- if (p == "p92") md$cluster92 else md$cluster51; u <- fdt_cluster_order_safe(k)
    do.call(rbind, lapply(ST, function(s) { v <- pe$cells$n_detected[, paste0("marker_", s)]
      data.frame(partition = p, cluster = u, stage = s, atlas_frac_ge1 = as.numeric(tapply(v >= 1, factor(k, levels = u), mean)),
                 atlas_frac_ge2 = as.numeric(tapply(v >= 2, factor(k, levels = u), mean)), stringsAsFactors = FALSE) })) }))
  mbvt_write_tsv(cc, file.path(att, "cluster_comparison.tsv"))
  au <- c(rep("GMP-Neut", 30), rep("CD14 Mono", 10), rep("HSC", 3), rep(NA, 6), rep("CD4 Naive", 6), rep("NK", 6))
  mbvt_write_tsv(nso_posthoc_composition(md, au), file.path(att, "posthoc_author_cluster_composition.tsv"))
  h5 <- file.path(w, "fake.h5ad"); writeLines("h5-bytes", h5)
  mbvt_write_tsv(data.frame(input = "h5ad", path = h5, sha256 = sha(h5), stringsAsFactors = FALSE), file.path(att, "input_sha256.tsv"))
  fl <- c("lossless.rds", "cluster_programs.tsv", "cluster_comparison.tsv", "posthoc_author_cluster_composition.tsv", "input_sha256.tsv")
  writeLines(jsonlite::toJSON(list(files = as.list(stats::setNames(vapply(file.path(att, fl), sha, ""), fl))), auto_unbox = TRUE, pretty = TRUE), file.path(att, "manifest.json"))
  U <- cbind(seq_len(n) %% 7 + as.integer(factor(cl)) * 10, seq_len(n) %% 5)
  g_calc <- function(what) if (!file.exists(file.path(out, "calc_manifest.json"))) stop(what, " 를 계산 고정 전에 읽었다")
  oid <- if (shuffle_ids) rev(ids) else ids
  reader <- list(umap = function() { g_calc("UMAP"); U }, obs = function(name) { g_calc("obs"); if (name != "cell_id") stop("obs: ", name); oid },
                 input_files = function() c(h5ad = h5))
  pins <- list(attempt = sha(file.path(att, "manifest.json")), h5ad = sha(h5))
  list(argv = c("--attempt-dir", att, "--data-root", w, "--out-dir", out), out = out, reader = reader, pins = pins, pe = pe, md = md)
}
acl_run <- function(w, pins = w$pins, argv = w$argv) { code <- NA_integer_; utils::capture.output(code <- atlas_cl_main(argv, root = ROOT, reader = w$reader, pins = pins)); code }
acl_rd <- function(out, f) utils::read.csv(file.path(out, f), colClasses = "character", na.strings = character(0), check.names = FALSE)

test_that("실행기 (소형 합성): 손계산 가중 차이 · 작은 군집 보존 · 0 · 동률 · 단계 최고 ≠ 군집 내 최고가 표에 그대로 남고, 계산 고정 뒤에만 UMAP 을 읽는다", {
  w <- acl_world(); expect_equal(acl_run(w), 0L); out <- w$out
  need <- c("cluster_scores.csv", "cluster_candidates.csv", "stage_peak_clusters.csv", "lossless.rds", "calc_manifest.json", "figure_text.csv",
            "posthoc_cluster_author_context.csv", "checks.csv", "input_sha256.csv", "environment.csv", "source_snapshot.csv", "manifest.json", "README_ko.md",
            file.path("figures", c("figA_cluster_top_candidate_umap.png", "figB_stage_cluster_scores_p92.png", "figB_stage_cluster_scores_p51.png")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  a <- log1p(1e4 * 5 / 100); b <- log1p(1e4 * 7 / 100)
  cd <- acl_rd(out, "cluster_candidates.csv"); A <- cd[cd$partition == "p92" & cd$cluster == "A", ]
  expect_equal(as.numeric(c(A$score_Myeloblast, A$score_PM, A$eb_Myeloblast, A$eb_PM)), c(30 * a / 40, 10 * b / 40, a / 2, b / 2), tolerance = 1e-12)
  expect_equal(c(A$score_top_candidate, A$eb_top_candidate, A$equal_block_agreement), c("Myeloblast", "PM", "different_top"))
  B <- cd[cd$partition == "p92" & cd$cluster == "B", ]
  expect_equal(c(B$n_cells, B$small_cluster, B$score_top_candidate, B$eb_top_status, B$equal_block_agreement), c("3", "TRUE", "Myeloblast", "unavailable", "equal_block_unavailable"))
  expect_true(nzchar(B$warning)); expect_equal(B$eb_Myeloblast, "NA")                                # 가용 block 이 없으면 결측 (pooled 로 바꾸지 않는다)
  C <- cd[cd$partition == "p92" & cd$cluster == "C", ]; expect_equal(c(C$score_top_status, C$score_top_candidate), c("unassigned_zero", ""))
  D <- cd[cd$partition == "p92" & cd$cluster == "D", ]; expect_equal(c(D$score_top_status, D$score_top_candidate, D$score_top_minus_second), c("tied", "MY|MM", "0"))
  expect_equal(sort(cd$cluster[cd$partition == "p92"]), c("A", "B", "C", "D", "E")); expect_equal(sort(cd$cluster[cd$partition == "p51"]), c("X", "Y"))
  pk <- acl_rd(out, "stage_peak_clusters.csv"); pm <- pk[pk$partition == "p92" & pk$stage == "PM", ]
  expect_equal(c(pm$peak_clusters, pm$peak_cluster_top_candidate, pm$peak_is_cluster_top), c("A", "Myeloblast", "FALSE"))
  sn <- pk[pk$partition == "p92" & pk$stage == "SN", ]; expect_equal(c(sn$peak_clusters, sn$n_clusters_with_stage_top), c("E", "0"))
  sc <- acl_rd(out, "cluster_scores.csv"); z <- sc[sc$partition == "p92" & sc$cluster == "E" & sc$stage == "SN", ]
  expect_equal(as.numeric(c(z$n_cells, z$frac_cells_ge1, z$frac_cells_ge2)), c(6, 1 / 6, 0)); expect_equal(nrow(sc), (5 + 2) * 6)
  ck <- acl_rd(out, "checks.csv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true(all(c("pooled_scores_recomputed_from_cells", "cluster_programs_tsv_matches_lossless", "detection_fractions_match_attempt", "umap_ids_identical_order",
                    "calc_fixed_before_umap_and_posthoc") %in% ck$check))
  cmf <- names(jsonlite::fromJSON(file.path(out, "calc_manifest.json"))$files); expect_false(any(grepl("posthoc|figure", cmf)))
  man <- jsonlite::fromJSON(file.path(out, "manifest.json"))$files
  expect_true(all(vapply(names(man), function(f) identical(digest::digest(file.path(out, f), algo = "sha256", file = TRUE), man[[f]]), TRUE)))
  ft <- acl_rd(out, "figure_text.csv"); expect_false(any(ft$text %in% c("NA", "NaN", "Inf"))); expect_lte(length(readLines(file.path(out, "README_ko.md"))), 20L)
  px <- acl_rd(out, "posthoc_cluster_author_context.csv"); expect_equal(px$top_author_label[px$partition == "p92" & px$cluster == "A"], "GMP-Neut")
})

test_that("실행기: 세포 ID 가 뒤섞이면 그림 전에 멈추고, 덮어쓰기 거절 · 고정 hash 불일치 · 인자 누락", {
  w <- acl_world(shuffle_ids = TRUE); expect_equal(acl_run(w), 4L); expect_true(any(grepl("세포 ID", readLines(file.path(w$out, "run.log")))))
  expect_false(file.exists(file.path(w$out, "figures", "figA_cluster_top_candidate_umap.png")))
  w2 <- acl_world(); expect_equal(acl_run(w2), 0L); expect_equal(acl_run(w2), 6L)
  w3 <- acl_world(); bad <- w3$pins; bad$attempt <- strrep("0", 64); expect_equal(acl_run(w3, pins = bad), 4L); expect_true(any(grepl("hash", readLines(file.path(w3$out, "run.log")))))
  w4 <- acl_world(); expect_equal(acl_run(w4, argv = w4$argv[1:2]), 2L)
})

test_that("실제 attempt_02: lossless 의 marker 프로그램 점수가 cluster_programs.tsv 와 숫자로 같고, 세포 평균에서 다시 계산해도 같다", {
  att <- file.path(ROOT, "results", "neutrophil_atlas_rna_2026-10-19", "attempt_02"); skip_if_not(file.exists(file.path(att, "lossless.rds")), "attempt_02 없음")
  ls <- readRDS(.mbvq_manifest_file(att, "manifest.json", "lossless.rds")); tsv <- utils::read.delim(.mbvq_manifest_file(att, "manifest.json", "cluster_programs.tsv"), stringsAsFactors = FALSE)
  sc <- atlas_cl_scores(ls$programs$clusters); k <- startsWith(tsv$program_id, "marker_"); tsv <- tsv[k, ]
  j <- match(paste(tsv$partition, tsv$cluster, sub("^marker_", "", tsv$program_id)), paste(sc$partition, sc$cluster, sc$stage))
  expect_false(anyNA(j)); expect_equal(nrow(sc), nrow(tsv)); expect_equal(sc$mean_log_pooled[j], tsv$mean_log_pooled, tolerance = 1e-14)
  expect_equal(sc$mean_log_equal_block[j], tsv$mean_log_equal_block, tolerance = 1e-14); expect_equal(sc$n_cells[j], tsv$n_cells)
  rc <- atlas_cl_recompute_pooled(ls$metadata, ls$programs$cells$mean_log); m <- match(paste(rc$partition, rc$cluster, rc$stage), paste(sc$partition, sc$cluster, sc$stage))
  expect_lt(max(abs(rc$mean_cells - sc$mean_log_pooled[m])), 1e-9); expect_equal(nrow(sc), (92 + 51) * 6)
})
