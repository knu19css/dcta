# =============================================================================
# 73_gse253355_umap.R — GSE253355 새 UMAP 표시 (표시 전용): 저장된 PCA 30 · kNN 20 을 그대로 uwot::umap 에 넣어 좌표를 계산하고, 저자 주석으로 그림 2 장
# =============================================================================
#
# 계약: Codex 전달 (2026-10-19, 사용자 새 요청 "전체 세포 UMAP 1 개 + HSC 부터 Neutrophil / Monocyte 계통만 추린 UMAP 1 개"). Ruling: claude_umap_progress.md
#   - 입력: gse253355_02 lossless.rds (sha256 고정) 의 graph_build$pca_scores (82,742 x 30) 와 graph_build$knn (자기 포함 k 20). PC 재척도화 · HVG 재선택 · 그래프 재적합 ·
#     점수 재계산을 하지 않는다. 이웃 거리는 저장 PCA 로 정확히 다시 계산한다(자기 자신 = 첫 열, 거리 0).
#   - 계산: uwot::umap(X = PCA, nn_method = list(idx, dist), n_neighbors 20, min_dist 0.3, spread 1, n_components 2, n_epochs 200, metric euclidean, scale FALSE,
#     init spectral, n_threads 2, n_sgd_threads 1, fast_sgd FALSE, seed 42) + set.seed(42). 실제 설치된 uwot 의 인자 이름을 확인하고, 필수 인자가 없으면 멈춘다.
#     uwot 이 없으면 멈춘다 (설치 승인 전에는 계산하지 않는다; DrL 레이아웃을 UMAP 이라 부르는 대체 없음). marker · DAG · 저자 라벨은 좌표 입력에 넣지 않는다.
#   - 표시: 저자 주석 cluster_anno_l2 (로컬 객체, sha256 고정; 세포 ID 정확 정렬)를 계산 뒤 색에만 쓴다. 그림 1 = 82,742 세포 전부 + 35 라벨 범례,
#     그림 2 = 지정 계통 라벨 세포만, 같은 전체 좌표의 확대(재학습 없음) + 라벨 · 범례 · 세포 수. 단위 축 비율, 연결선 · 방향 화살표 · 성숙 단계 가정 없음.
#   - --preflight TRUE: UMAP 없이 입력 검사만 기록한다.
# =============================================================================

GSU_SOURCES <- c(GSX_SOURCES_BASE <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                 "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                 "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R", "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R"),
                 "73_gse253355_umap.R")
GSU_LOSSLESS <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_02/lossless.rds"
GSU_OBJECT <- "data/external/gse253355/S2_steady_GSE253355.local_copy.rds"
GSU_PUBLIC <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_umap_01"
GSU_PINNED <- list(lossless = "c423ecfe33be76951acfc3e71ea0173e949c908d289839a19828e48b3d6bb215", object = "83fd17576c48e710e5a03a60be485beedbc088f82c6a97fa79d7b80d5367186a",
                   calc_manifest = "5842fd04cf6e0fa943538ac671c11bbc778e4439db71662f38e4fd46871c6955")                    # gse253355_02 calc_manifest.json (보존 확인)
GSU_SUBSET <- c("HSC", "MPP", "Cycling HSPC", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil", "Monocyte")
GSU_UMAP_ARGS <- list(n_neighbors = 20L, min_dist = 0.3, spread = 1, n_components = 2L, n_epochs = 200L, metric = "euclidean", scale = FALSE, init = "spectral",
                      n_threads = 2L, n_sgd_threads = 1L, fast_sgd = FALSE, seed = 42L)                         # uwot 0.2.x seed 인자 (Codex 인계서)
GSU_SEED <- 42L
GSU_UMAP_DEPS <- c("uwot", "Rcpp", "RcppAnnoy", "RcppProgress", "irlba", "FNN", "RSpectra", "dqrng", "RcppEigen", "Matrix")


# ── 1. 계산 (순수) ───────────────────────────────────────────────────────────
#' 저장 kNN 번호마다 저장 PCA 로 정확한 유클리드 거리 (n x k, 자기 자신 첫 열 = 0)
gsu_knn_dist <- function(Y, idx) {
  if (nrow(idx) != nrow(Y)) stop("kNN 행 수가 PCA 행 수와 다르다.")
  d <- vapply(seq_len(ncol(idx)), function(j) sqrt(rowSums((Y - Y[idx[, j], , drop = FALSE])^2)), numeric(nrow(Y))); if (is.null(dim(d))) d <- matrix(d, nrow(Y))
  d
}

#' uwot::umap 호출. umap_fun 이 없으면 설치된 uwot 을 쓴다 (없으면 멈춘다). 실제 API 에 필수 인자가 모두 있어야 한다.
gsu_umap_call <- function(Y, idx, dist, args, seed, umap_fun = NULL) {
  ver <- "injected test function"
  if (is.null(umap_fun)) { if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot 이 설치되어 있지 않다: 설치 승인 전에는 UMAP 을 계산하지 않는다 (DrL 등 다른 좌표로 대체하지 않는다).")
    umap_fun <- uwot::umap; ver <- as.character(utils::packageVersion("uwot")) }
  api <- names(formals(umap_fun)); need <- c("X", "nn_method", names(args)); miss <- setdiff(need, api)
  if (length(miss)) stop("UMAP 함수 API 에 필수 인자가 없다: ", paste(miss, collapse = ", "))
  call_args <- c(list(X = Y, nn_method = list(idx = idx, dist = dist)), args); if ("verbose" %in% api) call_args$verbose <- FALSE
  set.seed(seed); U <- do.call(umap_fun, call_args); U <- as.matrix(U)
  if (!identical(dim(U), c(nrow(Y), as.integer(args$n_components)))) stop("UMAP 좌표 크기가 다르다: ", paste(dim(U), collapse = " x "))
  rownames(U) <- rownames(Y); colnames(U) <- paste0("umap_", seq_len(ncol(U)))
  list(U = U, api = api, version = ver, passed = names(call_args), has_seed_arg = "seed" %in% api)
}

#' 부분집합: 정확한 저자 라벨만 (없는 라벨은 멈춘다). 반환: mask, counts (라벨 순서 = 지정 순서)
gsu_subset <- function(labels, subset) {
  miss <- setdiff(subset, labels); if (length(miss)) stop("지정한 저자 라벨이 자료에 없다 (정확한 이름 확인): ", paste(miss, collapse = ", "))
  list(mask = labels %in% subset, counts = data.frame(author_label = subset, n_cells = vapply(subset, function(l) sum(labels == l), 0L), row.names = NULL, stringsAsFactors = FALSE))
}

#' 입력 검사: PCA 유한 · 크기, 세포 ID 정확 정렬 (PCA 행 = lossless cells = 객체 meta.data), kNN 크기 · 자기 첫 이웃 · 범위, 라벨 결측 없음, 거리 유한 · 0 · 비감소
gsu_checks_inputs <- function(Y, idx, d, cells, meta_ids, labels) {
  n <- length(cells); r <- list(); add <- function(k, ok, det) r[[length(r) + 1L]] <<- data.frame(check = k, pass = isTRUE(ok), detail = det, stringsAsFactors = FALSE)
  add("pca_finite_dims", all(is.finite(Y)) && nrow(Y) == n, sprintf("PCA %d x %d, all finite", nrow(Y), ncol(Y)))
  add("cell_ids_exact_alignment", identical(rownames(Y), cells) && identical(as.character(meta_ids), cells) && !anyDuplicated(cells),
      sprintf("PCA rows == lossless cells == object meta.data rows, %d unique IDs in the same order", n))
  add("knn_self_first", is.integer(idx) && nrow(idx) == n && all(idx[, 1] == seq_len(n)) && all(idx >= 1L & idx <= n), sprintf("kNN %d x %d, self first, indices in 1..%d", nrow(idx), ncol(idx), n))
  add("knn_dist_exact", all(is.finite(d)) && all(d[, 1] == 0) && all(apply(d, 1, function(z) all(diff(z) >= -1e-9))), "exact Euclidean distances from stored PCA; self 0; non-decreasing per row")
  add("labels_complete", !anyNA(labels) && all(nzchar(labels)), sprintf("%d author labels, %d distinct", length(labels), length(unique(labels))))
  out <- do.call(rbind, r); if (!all(out$pass)) stop("입력 검사 실패: ", paste(out$check[!out$pass], collapse = ", ")); out
}


# ── 2. 그림 ───────────────────────────────────────────────────────────────────
GSU_LINEAGE_COLS <- c(HSC = "#F6222E", MPP = "#FEAF16", `Cycling HSPC` = "#85660D", GMP = "#3283FE", `Early Myeloid Progenitor` = "#AA0DFE", `Late Myeloid` = "#1CBE4F",
                      Neutrophil = "#FE00FA", Monocyte = "#1C7F93")                                       # 계통 8 라벨: 서로 뚜렷한 고정 색 (Polychrome 36 에서)
GSU_OTHER_COLS <- c("#5A5156", "#16FF32", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#B10DA1", "#FA0087", "#FC1CBF",
                    "#C075A6", "#782AB6", "#AAF400", "#822E1C", "#7ED7D1", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB", "#BCBD22", "#7F7F7F", "#BDCDFF", "#8C564B")   # 나머지: 흰 바탕에서 보이는 색만
#' 색 매핑: 계통 라벨은 고정 색, 나머지는 이름 순서로 남은 색 (두 그림이 같은 매핑을 쓴다; 색이 모자라면 멈춘다)
.gsu_cols <- function(labels, subset = GSU_SUBSET) {
  al <- sort(unique(labels)); lin <- intersect(names(GSU_LINEAGE_COLS), al); oth <- setdiff(al, lin)
  if (length(oth) > length(GSU_OTHER_COLS)) stop("라벨 수가 준비한 색보다 많다: ", length(al))
  cols <- c(GSU_LINEAGE_COLS[lin], stats::setNames(GSU_OTHER_COLS[seq_along(oth)], oth)); cols[al]
}
.gsu_order <- function(labels) { f <- table(labels); order(-as.integer(f[labels]), seq_along(labels)) }       # 큰 라벨부터 그려 작은 라벨이 위에 보이게

#' 그림 1. 전체 세포 (모든 점), 35 라벨 범례 (세포 수), 단위 축 비율
gsu_draw_full <- function(path, U, labels, cols, sub, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U); o <- .gsu_order(labels); tb <- table(labels)
  grDevices::png(path, width = 3200, height = 2450, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:2, 1), widths = c(3.3, 1.2)); graphics::par(oma = c(3.2, 0, 5.2, 0), family = fam)
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5)); mbvu_open_panel(b); graphics::points(U[o, 1], U[o, 2], pch = 16, cex = 0.24, col = cols[labels[o]]); graphics::box(col = "grey60")
    graphics::par(mar = c(0.5, 0, 0.5, 0.5)); graphics::plot.new(); lg <- sprintf("%s (%s)", names(cols), .nso_fmt_n(as.integer(tb[names(cols)])))
    graphics::legend("left", legend = lg, pch = 16, col = cols, pt.cex = 2.2, cex = 1.18, bty = "n", y.intersp = 1.12, title = if (ko) "저자 주석 cluster_anno_l2 (세포 수)" else "author annotation (cells)")
    for (k in lg) note("legend", "item", k)
    main <- if (ko) sprintf("GSE253355 골수 전체 %s 세포 — 새로 계산한 UMAP (저자 원 UMAP 아님)", .nso_fmt_n(nrow(U))) else sprintf("GSE253355 marrow, all %s cells - newly computed UMAP (not the author UMAP)", .nso_fmt_n(nrow(U)))
    graphics::mtext(main, side = 3, outer = TRUE, line = 2.6, cex = 2.0, font = 2); graphics::mtext(sub, side = 3, outer = TRUE, line = 0.8, cex = 1.15, col = "grey20")
    note("figure", "title", main); note("figure", "sub", sub)
    foot <- if (ko) "색 = 저자 세포 주석 (계산 뒤 표시에만 사용, 정답 라벨 아님). DCTA 단계 후보가 아니다. 좌표 입력에 marker · DAG · 저자 라벨을 넣지 않았다."
            else "colour = author cell annotation (display only, after computing; not ground truth). Not DCTA stage candidates. No marker, DAG or author label entered the embedding."
    graphics::mtext(foot, side = 1, outer = TRUE, line = 1.4, cex = 1.1, col = "grey20"); note("footer", "footer", foot)
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 그림 2. 지정 계통 라벨 세포만, 같은 전체 좌표의 확대 (재학습 없음) + 라벨 · 범례 · 세포 수. 반환: 글자 표와 그린 좌표
gsu_draw_subset <- function(path, U, labels, mask, cols, subset, sub, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; Us <- U[mask, , drop = FALSE]; ls_ <- labels[mask]; pad <- function(r) r + c(-1, 1) * 0.03 * diff(r)   # 모든 부분집합 세포가 창 안에
  b <- list(xlim = pad(range(Us[, 1])), ylim = pad(range(Us[, 2]))); o <- .gsu_order(ls_); tb <- table(factor(ls_, levels = subset))
  grDevices::png(path, width = 3200, height = 2450, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:2, 1), widths = c(3.3, 1.2)); graphics::par(oma = c(3.2, 0, 5.2, 0), family = fam)
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5)); mbvu_open_panel(b); graphics::points(Us[o, 1], Us[o, 2], pch = 16, cex = 0.42, col = cols[ls_[o]]); graphics::box(col = "grey60")
    rp <- atlas_cl_rep_cells(Us, ls_); rp$cluster <- sprintf("%s (%s)", rp$cluster, .nso_fmt_n(as.integer(tb[rp$cluster]))); cx <- 1.3
    lp <- .acl_label_layout(rp, cx, b, rp$n_cells); .acl_draw_labels(lp, lp$text, cx); if (any(!lp$free)) note("map", "labels_overlapping", paste(lp$text[!lp$free], collapse = " | "))
    graphics::par(mar = c(0.5, 0, 0.5, 0.5)); graphics::plot.new(); lg <- sprintf("%s (%s)", subset, .nso_fmt_n(as.integer(tb[subset])))
    graphics::legend("left", legend = lg, pch = 16, col = cols[subset], pt.cex = 2.6, cex = 1.35, bty = "n", y.intersp = 1.3, title = if (ko) "저자 주석 (세포 수)" else "author annotation (cells)")
    for (k in lg) note("legend", "item", k)
    main <- if (ko) sprintf("HSC 부터 Neutrophil · Monocyte 까지 계통의 저자 라벨 %d 개 — 같은 새 UMAP 좌표의 확대 (재학습 없음)", length(subset))
            else sprintf("Author labels from HSC to Neutrophil / Monocyte (%d) - zoom of the same new UMAP (no re-embedding)", length(subset))   # 방향 화살표 글자 없음
    sub2 <- sprintf(if (ko) "%s · 표시 세포 %s / %s" else "%s · cells shown %s / %s", sub, .nso_fmt_n(nrow(Us)), .nso_fmt_n(nrow(U)))
    graphics::mtext(main, side = 3, outer = TRUE, line = 2.6, cex = 2.0, font = 2); graphics::mtext(sub2, side = 3, outer = TRUE, line = 0.8, cex = 1.15, col = "grey20")
    note("figure", "title", main); note("figure", "sub", sub2)
    foot <- if (ko) "연결선 · 방향 화살표 · 성숙 단계 가정을 넣지 않았다. 다른 라벨 세포는 숨겼다(좌표는 전체 UMAP 그대로). 색 = 저자 주석 (정답 라벨 아님, DCTA 단계 후보 아님)."
            else "No connecting lines, arrows or maturation assumptions. Other labels hidden (coordinates unchanged from the full UMAP). Colour = author annotation (not ground truth, not DCTA stage candidates)."
    graphics::mtext(foot, side = 1, outer = TRUE, line = 1.4, cex = 1.1, col = "grey20"); note("footer", "footer", foot)
  }, finally = grDevices::dev.off())
  list(text = nf$table(), coords = Us)
}


# ── 3. 실행기 ────────────────────────────────────────────────────────────────
gsu_main <- function(argv, root = NULL, pins = GSU_PINNED, umap_fun = NULL) {
  keys <- c("--lossless" = "lossless", "--object" = "object", "--out-dir" = "out", "--public-dir" = "public", "--subset-labels" = "subset", "--preflight" = "preflight",
            "--author-col" = "author_col"); a <- list(); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (is.null(a$out) || is.na(a$out) || !nzchar(a$out)) { cat("--out-dir 가 필요하다.\n"); return(2L) }
  if (is.null(root)) root <- mdu_find_root()
  dflt <- list(lossless = GSU_LOSSLESS, object = GSU_OBJECT, public = GSU_PUBLIC); for (k in names(dflt)) if (is.null(a[[k]])) a[[k]] <- file.path(root, dflt[[k]])
  a$subset <- if (is.null(a$subset)) GSU_SUBSET else trimws(strsplit(a$subset, ";", fixed = TRUE)[[1]]); a$preflight <- isTRUE(as.logical(a$preflight)); if (is.null(a$author_col)) a$author_col <- "cluster_anno_l2"
  chk_dirs <- c(a$out, if (!a$preflight) a$public)
  for (d in chk_dirs) if (dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", d, "\n"); return(6L) }
  out <- a$out; dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", GSU_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(GSU_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = GSU_SOURCES)
    get("gsu_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins, umap_fun = umap_fun)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

gsu_build <- function(a, out, snap, log, argv, pins, umap_fun) {
  t0 <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"); rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  files <- character(0); crt <- logical(0); wc <- function(df, f) { p0 <- file.path(out, f); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  hl <- .nrv_sha(a$lossless); ho <- .nrv_sha(a$object); cmf0 <- file.path(dirname(a$lossless), "calc_manifest.json"); hc0 <- if (file.exists(cmf0)) .nrv_sha(cmf0) else NA_character_
  if (!is.null(pins$calc_manifest) && !identical(hc0, pins$calc_manifest)) stop("기존 calc_manifest.json hash 가 고정값과 다르다: ", hc0)
  if (!is.null(pins$lossless) && !identical(hl, pins$lossless)) stop("lossless.rds hash 가 고정값과 다르다: ", hl)
  if (!is.null(pins$object) && !identical(ho, pins$object)) stop("객체 hash 가 고정값과 다르다: ", ho)
  ls <- readRDS(a$lossless); Y <- ls$graph_build$pca_scores; idx <- ls$graph_build$knn; cells <- ls$cells
  if (is.null(Y) || is.null(idx)) stop("lossless 에 graph_build$pca_scores / knn 이 없다.")
  ob <- gsx_open_rds(a$object)$object; md <- gsx_meta(ob); rm(ob); invisible(gc())
  if (!a$author_col %in% names(md)) stop("meta.data 에 저자 라벨 열이 없다: ", a$author_col)
  lab <- as.character(md[[a$author_col]]); d <- gsu_knn_dist(Y, idx)
  ci <- gsu_checks_inputs(Y, idx, d, cells, rownames(md), lab); ss <- gsu_subset(lab, a$subset)
  log("inputs checked:", nrow(Y), "cells x", ncol(Y), "PCs; kNN", ncol(idx), "| subset", sum(ss$mask), "cells in", length(a$subset), "labels")
  envr <- data.frame(item = c("R_executable", "R", "platform", "libPaths", "umap_packages", "umap_args", "seed", "command", "written_at"),
                     value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                               paste(GSU_UMAP_DEPS, vapply(GSU_UMAP_DEPS, function(p) if (nzchar(system.file(package = p))) as.character(utils::packageVersion(p)) else "absent", ""), sep = "=", collapse = " "),
                               paste(names(GSU_UMAP_ARGS), unlist(GSU_UMAP_ARGS), sep = "=", collapse = " "), paste("set.seed", GSU_SEED), paste(c("73_gse253355_umap.R", argv), collapse = " "), stamp()),
                     stringsAsFactors = FALSE)
  if (a$preflight) { wc(ci, "preflight_checks.csv"); wc(ss$counts, "preflight_subset_counts.csv"); wc(envr, "preflight_environment.csv")
    log("preflight only: UMAP not computed (uwot", if (nzchar(system.file(package = "uwot"))) "present" else "absent", ")"); return(invisible(TRUE)) }
  # ── UMAP (표시 좌표만) ──
  res <- gsu_umap_call(Y, idx, d, GSU_UMAP_ARGS, GSU_SEED, umap_fun); U <- res$U
  ck <- rbind(ci, data.frame(check = c("umap_finite_dims", "all_cells_kept", "umap_args_recorded"),
                             pass = c(all(is.finite(U)) && ncol(U) == 2L, identical(rownames(U), cells) && nrow(U) == length(cells),
                                      all(c("X", "nn_method", names(GSU_UMAP_ARGS)) %in% res$passed)),
                             detail = c(sprintf("%d x %d, all finite", nrow(U), ncol(U)), sprintf("all %d cells, same IDs and order; none added or dropped", nrow(U)),
                                        sprintf("uwot %s; passed %s; set.seed(%d); API has seed argument: %s", res$version, paste(res$passed, collapse = ","), GSU_SEED, res$has_seed_arg)),
                             stringsAsFactors = FALSE))
  log("UMAP computed:", res$version)
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; cols <- .gsu_cols(lab, a$subset); fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  sub <- sprintf(if (ko) "uwot::umap %s · 입력 = 저장 PCA %d · kNN %d (자기 포함, 정확 거리) · min_dist %s · n_epochs %d · seed %d · DrL 레이아웃과 다른 새 좌표" else
                   "uwot::umap %s · input = stored PCA %d · kNN %d (self-including, exact distances) · min_dist %s · n_epochs %d · seed %d · new coordinates, not the DrL layout",
                 res$version, ncol(Y), ncol(idx), format(GSU_UMAP_ARGS$min_dist), GSU_UMAP_ARGS$n_epochs, GSU_SEED)
  f1 <- file.path(fig, "fig_umap_full_author.png"); t1 <- gsu_draw_full(f1, U, lab, cols, sub, ko, fam)
  f2 <- file.path(fig, "fig_umap_hsc_neutrophil_monocyte.png"); r2 <- gsu_draw_subset(f2, U, lab, ss$mask, cols, a$subset, sub, ko, fam); files <- c(files, f1, f2)
  ck <- rbind(ck, data.frame(check = "subset_coordinates_identical", pass = identical(unname(r2$coords), unname(U[ss$mask, , drop = FALSE])),
                             detail = sprintf("subset figure plots the %d subset cells at their full-UMAP coordinates (no re-embedding)", sum(ss$mask)), stringsAsFactors = FALSE))
  wc(rbind(t1, r2$text), "figure_text.csv")
  ac <- data.frame(author_label = names(cols), n_cells = as.integer(table(lab)[names(cols)]), in_subset = names(cols) %in% a$subset, colour = unname(cols), stringsAsFactors = FALSE); wc(ac, "annotation_counts.csv")
  p0 <- file.path(out, "umap_coordinates.rds")
  saveRDS(list(cells = cells, U = U, new_cluster = ls$clusters, author_label = lab, subset_labels = a$subset, subset_mask = ss$mask, umap_args = GSU_UMAP_ARGS, seed = GSU_SEED,
               uwot_version = res$version, uwot_api = res$api, input_sha256 = c(lossless = hl, object = ho)), p0, version = 3); files <- c(files, p0)
  p0 <- file.path(out, "umap_coordinates_local.csv.gz"); con <- gzfile(p0, "w")
  utils::write.csv(data.frame(cell_id = cells, umap_1 = sprintf("%.17g", U[, 1]), umap_2 = sprintf("%.17g", U[, 2]), new_cluster = ls$clusters, author_label = lab), con, row.names = FALSE); close(con)
  files <- c(files, p0)
  wc(rbind(envr, data.frame(item = c("uwot_version", "uwot_api"), value = c(res$version, paste(res$api, collapse = ",")), stringsAsFactors = FALSE)), "environment.csv")
  keep_ok <- identical(.nrv_sha(a$lossless), hl) && identical(.nrv_sha(a$object), ho) && (is.na(hc0) || identical(.nrv_sha(cmf0), hc0))
  ck <- rbind(ck, data.frame(check = c("existing_inputs_results_unchanged", "tables_roundtrip_ok"), pass = c(keep_ok, all(crt)),
                             detail = c(sprintf("lossless %s…, object %s…, existing calc_manifest %s unchanged after the run", substr(hl, 1, 12), substr(ho, 1, 12), if (is.na(hc0)) "(none)" else paste0(substr(hc0, 1, 12), "…")),
                                        sprintf("%d CSVs re-read identically", length(crt))), stringsAsFactors = FALSE)); wc(ck, "checks.csv")
  sm <- c("# GSE253355 새 UMAP (표시 전용) — 요약", "",
          sprintf("- 입력: gse253355_02 lossless.rds (sha256 %s…) 의 저장 PCA %d · kNN %d, 로컬 객체 (sha256 %s…) 의 저자 주석 %s.", substr(hl, 1, 12), ncol(Y), ncol(idx), substr(ho, 1, 12), a$author_col),
          sprintf("- 계산: uwot::umap %s, %s, set.seed(%d). 이웃 거리는 저장 PCA 로 정확 계산(자기 0). 재척도화 · HVG · 그래프 · 점수 재계산 없음.", res$version,
                  paste(names(GSU_UMAP_ARGS), unlist(GSU_UMAP_ARGS), sep = "=", collapse = ", "), GSU_SEED),
          sprintf("- 그림 1: 전체 %s 세포, 저자 라벨 %d 개 범례. 그림 2: %s (%s 세포), 같은 좌표의 확대.", .nso_fmt_n(nrow(U)), length(cols), paste(a$subset, collapse = " · "), .nso_fmt_n(sum(ss$mask))),
          paste0("- 계통 라벨 세포 수: ", paste(sprintf("%s %s", ss$counts$author_label, .nso_fmt_n(ss$counts$n_cells)), collapse = " · "), "."),
          "- 이 좌표는 새로 계산한 UMAP 이다: 저자 원 UMAP 도, gse253355_02 의 DrL 레이아웃도 아니다. 연결선 · 방향 · 성숙 단계 가정을 넣지 않았다.",
          "- 색은 저자 주석(정답 라벨 아님)이며 DCTA 단계 후보와 다르다. 세포 좌표 · 군집 번호는 results 의 umap_coordinates*.{rds,csv.gz} 에만 둔다(공개 · git 제외).",
          sprintf("- 검사: %s.", paste(ck$check, ck$pass, sep = "=", collapse = ", ")))
  p0 <- file.path(out, "summary_ko.md"); writeLines(sm, p0); files <- c(files, p0)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), inputs = list(lossless = hl, object = ho), umap_args = GSU_UMAP_ARGS, seed = GSU_SEED, uwot_version = res$version,
                                   note = "new UMAP display coordinates; not the author UMAP; differs from the DrL layout of gse253355_02",
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "manifest.json"))
  if (!all(ck$pass)) stop("확인 실패: ", paste(ck$check[!ck$pass], collapse = ", "))
  pub <- a$public; dir.create(pub, recursive = TRUE, showWarnings = FALSE); pf <- c(f1, f2, file.path(out, c("summary_ko.md", "annotation_counts.csv")))
  if (!all(file.copy(pf, file.path(pub, basename(pf)), overwrite = FALSE))) stop("공개 사본 복사 실패")
  writeLines(jsonlite::toJSON(list(source_output = normalizePath(out), files = as.list(stats::setNames(vapply(file.path(pub, basename(pf)), .nrv_sha, ""), basename(pf)))), auto_unbox = TRUE, pretty = TRUE),
             file.path(pub, "public_manifest.json"))
  log("done | checks:", paste(ck$check, ck$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .gsu_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.gsu_bb, .f))
  quit(status = gsu_main(commandArgs(trailingOnly = TRUE)))
}
