# =============================================================================
# 74_gse253355_marker_umap.R — 실제 R UMAP (gse253355_umap_03) 에 고정 6 단계 RNA 패널 군집 점수 (gse253355_02) 대응 (표시 · 집계 전용)
#   그림 1 · 2 = 저자 주석 vs 군집별 최고 마커 단계 (전체 · 계통 부분집합), 그림 3 = 군집 x 6 단계 원점수 heatmap, 군집 x 단계 x 유전자 기여 표
# =============================================================================
#
# 계약: Codex 인계서 codex_marker_umap_instruction.md (2026-10-19, 사용자 "R로 그려봐 · 우리 마커들이랑 대응시켜봐"). Ruling: claude_marker_umap_progress.md
#   - 좌표: gse253355_umap_03 umap_coordinates.rds (sha256 고정) 의 cells · U · new_cluster · author_label · subset_mask 를 그대로 쓴다 (UMAP 재계산 · DrL 재사용 없음).
#   - 점수 · 후보: gse253355_02 의 cluster_scores.csv · cluster_candidates.csv · stage_peak_clusters.csv (calc_manifest.json 기록 hash 와 같아야 한다) 그대로.
#   - 유전자 기여: 원 count (로컬 객체, sha256 고정, 읽기 전용) 에서 R69 gsx_panel_map · gsx_cell_scores 로 log1p(1e4 · count / 모든 유전자 library) 를 다시 계산하고,
#     패널 대응표 · 세포 점수 · library 가 lossless 와 identical 일 때만 군집 x 단계 x 유전자 평균 · 검출률 · 기여 (= 평균 / 그 단계 목록 유전자 수, 실제 패널 10) 를 쓴다.
#     없는 유전자는 발현 · 검출 NA, 기여 0. 기여 합 = 기존 군집 점수 (허용 1e-10). 겹치는 유전자는 단계마다 한 행 + listed_in_stages 진단.
#   - 저자 라벨은 계산 고정 (calc_manifest.json) 뒤 표시 · 사후 대조에만 쓴다 (정답 · 정확도 주장 없음). 군집 안 1 위 (winner) 와 단계별 전체 군집 최고 (peak) 를 구분한다.
#   - 그림은 R base graphics (png). 저자 색 = R73 `.gsu_cols`, 단계 색 = R67 `ATLAS_CL_STAGE_COL`. 연결선 · 화살표 없음.
#   - 결과 폴더 · 공개 폴더가 비어 있지 않으면 쓰지 않는다 (exit 6). 공개 사본 = PNG · 집계 표 · 요약 (세포 단위 · 좌표 없음).
#
# 실행 (W 에서, dcta 환경만)
#   run-r dcta R/branch_biology/74_gse253355_marker_umap.R --out-dir results/neutrophil_atlas_rna_2026-10-19/gse253355_marker_umap_01 [--public-dir …] [--focus-clusters "6;12;15;1;25"]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

GSM_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                 "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                 "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R", "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R",
                 "73_gse253355_umap.R", "74_gse253355_marker_umap.R")
GSM_UMAP <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_umap_03/umap_coordinates.rds"
GSM_SCORES_DIR <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_02"
GSM_OBJECT <- "data/external/gse253355/S2_steady_GSE253355.local_copy.rds"
GSM_PANELS <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/atlas_stage_markers.csv"
GSM_HGNC <- "data/external/granulopoiesis_atlas_2026-10-19/hgnc_complete_set.txt"
GSM_PUBLIC <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_marker_umap_01"
GSM_PINNED <- list(umap = "1f91bc928f06b736980bd56e18a595346ff1ddc430a96a021e7537149ab3001f", lossless = "c423ecfe33be76951acfc3e71ea0173e949c908d289839a19828e48b3d6bb215",
                   object = "83fd17576c48e710e5a03a60be485beedbc088f82c6a97fa79d7b80d5367186a", calc_manifest = "5842fd04cf6e0fa943538ac671c11bbc778e4439db71662f38e4fd46871c6955")
GSM_FOCUS <- c("6", "12", "15", "1", "25")                                                             # 인계서의 주요 mismatch 군집
GSM_TOP_GENES <- 5L
GSM_LABEL_MIN_FRAC <- 0.01                                                                             # Ruling M1: 그림 2 번호 = 부분집합 세포의 1 % 이상 군집
GSM_CONTRIB_TOL <- 1e-10
GSM_FIXED_FILES <- c("cluster_scores.csv", "cluster_candidates.csv", "stage_peak_clusters.csv", "lossless.rds")   # gse253355_02 calc_manifest.json 기록과 대조
GSM_PUBLIC_FILES <- c("summary_ko.md", "cluster_marker_summary.csv", "cluster_stage_gene_contributions.csv", "focus_cluster_top_genes.csv", "panel_gene_diagnostics.csv",
                      "contribution_check.csv", "posthoc_author_composition.csv", "checks.csv", "figure_text.csv")


# ── 1. 계산 (순수) ───────────────────────────────────────────────────────────
#' 군집 x 단계 x 유전자 표. Lg / Dg = 측정 패널 유전자 x 세포 (log1p CP10k / 검출 0 · 1, 행 이름 = symbol), pm = gsx_panel_map(), cluster = 세포 군집 (Lg 열 순서),
#' clus = 기존 군집 점수 긴 표 (cluster · stage · mean_all_cells; 군집 · 단계 순서 = 이 표). 기여 = 평균 / 그 단계 목록 유전자 수. 없는 유전자: 발현 · 검출 NA, 기여 0
gsm_gene_table <- function(Lg, Dg, pm, cluster, clus) {
  cluster <- as.character(cluster); if (length(cluster) != ncol(Lg) || !identical(dim(Dg), dim(Lg))) stop("세포 수 · 유전자 행렬 크기가 맞지 않는다.")
  u <- unique(as.character(clus$cluster)); st <- unique(as.character(clus$stage)); if (!all(cluster %in% u)) stop("점수 표에 없는 군집의 세포가 있다.")
  f <- factor(cluster, levels = u); nc <- tabulate(f, nbins = length(u)); if (any(nc == 0L)) stop("세포가 없는 군집이 점수 표에 있다: ", paste(u[nc == 0L], collapse = ", "))
  ind <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), length(u)))
  ML <- sweep(as.matrix(Lg %*% ind), 2, nc, "/"); MD <- sweep(as.matrix(Dg %*% ind), 2, nc, "/"); dimnames(ML) <- dimnames(MD) <- list(rownames(Lg), u)
  pm <- pm[pm$stage %in% st, , drop = FALSE]; meas <- pm$status == "exact"; gi <- match(pm$symbol, rownames(ML))
  if (any(meas & is.na(gi))) stop("측정된 패널 유전자가 유전자 행렬에 없다: ", paste(pm$symbol[meas & is.na(gi)], collapse = ", "))
  nl <- as.integer(table(factor(pm$stage, levels = st))[pm$stage]); li <- vapply(pm$symbol, function(g) paste(st[st %in% pm$stage[pm$symbol == g]], collapse = "|"), "", USE.NAMES = FALSE)
  sc <- stats::setNames(as.numeric(clus$mean_all_cells), paste(clus$cluster, clus$stage, sep = "\r"))
  d <- do.call(rbind, lapply(seq_along(u), function(j) { ml <- ifelse(meas, ML[cbind(gi, j)], NA_real_); md <- ifelse(meas, MD[cbind(gi, j)], NA_real_)
    data.frame(cluster = u[j], n_cells = nc[j], stage = pm$stage, rank = as.integer(pm$rank), gene = pm$symbol, ensg = if ("ensg" %in% names(pm)) as.character(pm$ensg) else "",
               gene_status = ifelse(meas, "measured", "missing"), listed_in_stages = li, n_listed_in_stage = nl, mean_log1p_cp10k = ml, frac_cells_detected = md,
               contribution = ifelse(meas, ml / nl, 0), stage_score = unname(sc[paste(u[j], pm$stage, sep = "\r")]), stringsAsFactors = FALSE) }))
  d$share_of_stage_score <- ifelse(!is.na(d$stage_score) & d$stage_score > 0, d$contribution / d$stage_score, NA_real_)
  d <- d[order(match(d$cluster, u), match(d$stage, st), d$rank), , drop = FALSE]; rownames(d) <- NULL
  d$rank_in_stage <- as.integer(stats::ave(seq_len(nrow(d)), d$cluster, d$stage, FUN = function(ix) { o <- order(-d$contribution[ix], d$rank[ix]); r <- integer(length(ix)); r[o] <- seq_along(ix); r }))
  if (!any(nzchar(d$ensg))) d$ensg <- NULL
  d
}

#' 기여 합 검사: 군집 x 단계마다 유전자 기여 합과 기존 군집 점수 (mean_all_cells) 의 차이 (허용 tol)
gsm_contribution_check <- function(gt, clus, tol = GSM_CONTRIB_TOL) {
  s <- tapply(gt$contribution, paste(gt$cluster, gt$stage, sep = "\r"), sum); cs <- as.numeric(s[paste(clus$cluster, clus$stage, sep = "\r")])
  d <- data.frame(cluster = as.character(clus$cluster), stage = as.character(clus$stage), cluster_score = as.numeric(clus$mean_all_cells), contribution_sum = cs, stringsAsFactors = FALSE)
  d$abs_diff <- abs(d$contribution_sum - d$cluster_score); d$pass <- !is.na(d$abs_diff) & d$abs_diff <= tol; d
}

#' 저자 구성 (사후): 군집마다 많은 순서 (R72 posthoc_author_composition 과 같은 정렬 · 결측 이름)
gsm_author_composition <- function(cluster, author, clusters) {
  cluster <- as.character(cluster); author <- as.character(author); author[is.na(author) | !nzchar(author)] <- "(no author label)"
  do.call(rbind, lapply(as.character(clusters), function(k) { w <- cluster == k; tb <- sort(table(author[w]), decreasing = TRUE)
    data.frame(cluster = k, author_label = names(tb), n_cells = as.integer(tb), frac_of_cluster = as.numeric(tb) / sum(w), stringsAsFactors = FALSE) }))
}

#' 군집 요약: 저자 구성 (사후) + 6 단계 점수 + 군집 안 1 위 (winner) · 2 위 · 차이 · donor 동일 가중 일치 + 단계별 전체 군집 최고 (peak) 인 단계
gsm_cluster_summary <- function(cand, peaks, comp, subset_n) {
  cl <- as.character(cand$cluster); sp <- function(x) { z <- strsplit(x, "|", fixed = TRUE)[[1]]; z[nzchar(z)] }
  pk <- lapply(cl, function(k) peaks$stage[vapply(peaks$peak_clusters, function(z) k %in% sp(z), TRUE)]); win <- lapply(cand$score_top_candidate, sp)
  au <- lapply(cl, function(k) comp[comp$cluster == k, , drop = FALSE]); pick <- function(i, col, na) vapply(au, function(z) if (nrow(z) >= i) z[[col]][i] else na, na)
  d <- data.frame(cluster = cl, n_cells = as.integer(cand$n_cells), n_donors = as.integer(cand$n_donors), n_subset_cells = as.integer(subset_n[cl]),
                  author_top_label = pick(1L, "author_label", ""), author_top_n = pick(1L, "n_cells", NA_integer_), author_top_frac = pick(1L, "frac_of_cluster", NA_real_),
                  author_second_label = pick(2L, "author_label", ""), author_second_frac = pick(2L, "frac_of_cluster", NA_real_), author_n_labels = vapply(au, nrow, 0L),
                  author_composition = vapply(au, function(z) { z <- utils::head(z, 3L); paste(sprintf("%s %d (%d%%)", z$author_label, z$n_cells, round(100 * z$frac_of_cluster)), collapse = " | ") }, ""),
                  stringsAsFactors = FALSE)
  for (s in NSO_STAGES) d[[paste0("score_", s)]] <- as.numeric(cand[[paste0("score_", s)]])
  d$winner_status <- cand$score_top_status; d$winner_stage <- cand$score_top_candidate; d$winner_score <- as.numeric(cand$score_top_value)
  d$second_stage <- cand$score_second_candidate; d$second_score <- as.numeric(cand$score_second_value); d$winner_minus_second <- as.numeric(cand$score_top_minus_second)
  d$donor_equal_top_candidate <- cand$donor_equal_top_candidate; d$donor_equal_agreement <- cand$donor_equal_agreement
  d$stage_peak_of <- vapply(pk, paste, "", collapse = "|"); d$winner_is_stage_peak <- vapply(seq_along(cl), function(i) any(win[[i]] %in% pk[[i]]), TRUE)
  d$note <- rep(paste("winner = within-cluster top stage of the cluster-mean panel score (fixed denominator); stage_peak_of = stages for which this cluster has the highest score",
                      "among all clusters; author labels are post hoc display only (not truth, no accuracy claim)"), nrow(d))
  d
}

#' mismatch 군집 기여 유전자: 군집마다 1 위 · 2 위 · peak 단계에서 기여 > 0 상위 n_top 개 (기여 큰 순서, 같으면 패널 순위), 단계 점수 대비 비율 · 누적 비율
gsm_focus_genes <- function(gt, summ, focus, n_top = GSM_TOP_GENES) {
  miss <- setdiff(focus, summ$cluster); if (length(miss)) stop("지정한 군집이 없다: ", paste(miss, collapse = ", "))
  sp <- function(x) { z <- strsplit(x, "|", fixed = TRUE)[[1]]; z[nzchar(z)] }; rows <- list()
  for (k in focus) { i <- match(k, summ$cluster); roles <- list(winner = sp(summ$winner_stage[i]), second = sp(summ$second_stage[i]), peak = sp(summ$stage_peak_of[i]))
    for (s in unique(unlist(roles))) { z <- gt[gt$cluster == k & gt$stage == s & gt$contribution > 0, , drop = FALSE]; if (!nrow(z)) next
      z <- utils::head(z[order(-z$contribution, z$rank), , drop = FALSE], n_top)
      rows[[length(rows) + 1L]] <- data.frame(cluster = k, author_top_label = summ$author_top_label[i], author_top_frac = summ$author_top_frac[i],
                                              role = paste(names(roles)[vapply(roles, function(r) s %in% r, TRUE)], collapse = "+"), stage = s, stage_score = z$stage_score,
                                              gene_order = seq_len(nrow(z)), gene = z$gene, listed_in_stages = z$listed_in_stages, mean_log1p_cp10k = z$mean_log1p_cp10k,
                                              frac_cells_detected = z$frac_cells_detected, contribution = z$contribution, share_of_stage_score = z$share_of_stage_score,
                                              cumulative_share = cumsum(z$share_of_stage_score), stringsAsFactors = FALSE) } }
  if (!length(rows)) stop("지정한 군집에 기여 > 0 유전자가 없다.")
  d <- do.call(rbind, rows); rownames(d) <- NULL; d
}


# ── 2. 그림 (R base graphics) ────────────────────────────────────────────────
.gsm_stage_cols <- function(cand) stats::setNames(.acl_top_cols(cand), as.character(cand$cluster))
.gsm_short <- function(cand) stats::setNames(ifelse(cand$score_top_status == "score_top_candidate", cand$score_top_candidate,
                                                    ifelse(cand$score_top_status == "tied", gsub("|", "/", cand$score_top_candidate, fixed = TRUE), "0")), as.character(cand$cluster))

#' 그림 1 · 2. 왼쪽 = 저자 주석 (R73 색), 오른쪽 = 군집 1 위 단계 색 + 군집 번호. mask = NULL 이면 전체 세포, 아니면 그 세포만 같은 좌표로 (창 = 그 세포 범위 + 3 %).
#' 번호: 전체 = 모든 군집, 부분집합 = 부분집합 세포의 label_min_frac 이상 군집에 "번호 · 단계". 반환: 글자 표 · 그린 좌표 · 군집별 색 확인 · 번호 군집 · 번호 없는 군집
gsm_draw_pair <- function(path, U, author, acols, cluster, cand, mask = NULL, subset_labels = NULL, label_min_frac = 0, sub = "", ko = TRUE, fam = "") {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; cluster <- as.character(cluster); full <- is.null(mask); keep <- if (full) rep(TRUE, nrow(U)) else as.logical(mask)
  Us <- U[keep, , drop = FALSE]; au <- as.character(author[keep]); cl <- cluster[keep]; ns <- nrow(Us); fmt <- .nso_fmt_n; pad <- function(r) r + c(-1, 1) * 0.03 * diff(r)
  b <- if (full) mbvu_bounds(U) else list(xlim = pad(range(Us[, 1])), ylim = pad(range(Us[, 2])))
  scol <- .gsm_stage_cols(cand); if (!all(cl %in% names(scol))) stop("후보 표에 없는 군집의 세포가 있다."); pcol <- unname(scol[cl]); short <- .gsm_short(cand)
  cnt <- table(factor(cl, levels = as.character(cand$cluster))); shown <- names(cnt)[cnt > 0]
  numbered <- if (full) shown else shown[as.integer(cnt[shown]) >= label_min_frac * ns]; unn <- setdiff(shown, numbered)
  colours <- data.frame(cluster = shown, winner = cand$score_top_candidate[match(shown, cand$cluster)], status = cand$score_top_status[match(shown, cand$cluster)],
                        n_cells = as.integer(cnt[shown]), n_colours = vapply(shown, function(k) length(unique(pcol[cl == k])), 0L, USE.NAMES = FALSE),
                        colour = vapply(shown, function(k) pcol[cl == k][1], "", USE.NAMES = FALSE), stringsAsFactors = FALSE)
  pcex <- if (full) 0.28 else 0.42                                                                     # marker_umap_02 그림 검토: 점 · 글자가 작았다 (Ruling M9)
  grDevices::png(path, width = 4400, height = 2100, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:4, 1), widths = c(3, 1.4, 3, 1.4)); graphics::par(oma = c(5.4, 0, 7.6, 0), family = fam)
    # 왼쪽: 저자 주석
    graphics::par(mar = c(0.5, 0.5, 5.6, 0.5)); mbvu_open_panel(b); o <- .gsu_order(au)
    graphics::points(Us[o, 1], Us[o, 2], pch = 16, cex = pcex, col = acols[au[o]]); graphics::box(col = "grey60")
    if (!full) { rp <- atlas_cl_rep_cells(Us, au); cx <- 1.7; lp <- .acl_label_layout(rp, cx, b, rp$n_cells); .acl_draw_labels(lp, lp$text, cx)
      for (k in lp$text) note("author", "label", k); if (any(!lp$free)) note("author", "labels_overlapping", paste(lp$text[!lp$free], collapse = " | ")) }
    t1 <- if (ko) "저자 주석 (cluster_anno_l2)" else "Author annotation (cluster_anno_l2)"
    s1 <- sprintf(if (ko) "라벨 %d 개 · 사후 대조용 (정답 아님)" else "%d labels · post hoc comparison (not truth)", length(unique(au)))
    graphics::title(main = t1, cex.main = 2.15, line = 2.9); graphics::mtext(s1, side = 3, line = 0.8, cex = 1.4, col = "grey20"); note("author", "title", t1); note("author", "sub", s1)
    graphics::par(mar = c(0.5, 0, 5.6, 0.5)); graphics::plot.new(); labs <- if (full) names(acols)[names(acols) %in% au] else subset_labels; tb <- table(au)
    lg <- sprintf("%s (%s)", labs, fmt(as.integer(tb[labs])))
    graphics::legend("left", legend = lg, pch = 16, col = acols[labs], pt.cex = if (full) 2.4 else 3.0, cex = if (full) 1.7 else 1.95, bty = "n", y.intersp = if (full) 1.0 else 1.3,
                     title = if (ko) "저자 라벨 (세포 수)" else "author label (cells)", title.adj = 0)
    for (k in lg) note("author_legend", "item", k)
    # 오른쪽: 군집 1 위 단계 + 번호
    graphics::par(mar = c(0.5, 0.5, 5.6, 0.5)); mbvu_open_panel(b); o2 <- .gsu_order(cl)
    graphics::points(Us[o2, 1], Us[o2, 2], pch = 16, cex = pcex, col = pcol[o2]); graphics::box(col = "grey60")
    rc <- atlas_cl_rep_cells(Us, cl); rc <- rc[rc$cluster %in% numbered, , drop = FALSE]; if (!full) rc$cluster <- sprintf("%s · %s", rc$cluster, short[rc$cluster])
    cx2 <- if (full) 1.8 else 1.95; lc <- .acl_label_layout(rc, cx2, b, rc$n_cells); .acl_draw_labels(lc, lc$text, cx2)
    for (k in lc$text) note("marker", "label", k); if (any(!lc$free)) note("marker", "labels_overlapping", paste(lc$text[!lc$free], collapse = " | "))
    t2 <- if (ko) "군집별 최고 마커 단계 (6 단계 RNA 패널)" else "Top marker stage per cluster (6-stage RNA panel)"
    s2 <- if (full) sprintf(if (ko) "군집 %d 개 번호 · 색 = 군집 평균 점수 1 위 단계" else "%d clusters numbered · colour = top stage of the cluster-mean score", length(shown))
          else if (ko) sprintf("번호 = 이 세포의 %s%% 이상인 군집 %d 개 · 번호 없는 군집 %d 개 (%s 세포) 는 색만", format(100 * label_min_frac), length(numbered), length(unn), fmt(sum(cnt[unn])))
          else sprintf("numbered = %d clusters with >= %s%% of these cells · %d clusters (%s cells) colour only", length(numbered), format(100 * label_min_frac), length(unn), fmt(sum(cnt[unn])))
    graphics::title(main = t2, cex.main = 2.15, line = 2.9); graphics::mtext(s2, side = 3, line = 0.8, cex = 1.4, col = "grey20"); note("marker", "title", t2); note("marker", "sub", s2)
    graphics::par(mar = c(0.5, 0, 5.6, 0.5)); graphics::plot.new()
    stl <- c(NSO_STAGES, "tied", "zero"); stc <- c(ATLAS_CL_STAGE_COL[NSO_STAGES], tied = ATLAS_CL_TIE_COL, zero = ATLAS_CL_NONE_COL)
    grp <- ifelse(colours$status == "score_top_candidate", colours$winner, ifelse(colours$status == "tied", "tied", "zero"))
    keep_s <- stl[stl %in% NSO_STAGES | stl %in% grp]
    items <- vapply(keep_s, function(s) { k <- colours$cluster[grp == s]; nn <- sum(colours$n_cells[grp == s]); nm <- if (s == "tied") (if (ko) "동률" else "tied") else if (s == "zero") (if (ko) "0 점" else "zero") else s
      if (full) sprintf(if (ko) "%s  %d 군집 · %s 세포" else "%s  %d clusters · %s cells", nm, length(k), fmt(nn))
      else { kn <- intersect(k, numbered); sprintf(if (ko) "%s  %s · %s 세포" else "%s  %s · %s cells", nm, if (length(kn)) paste0(if (ko) "군집 " else "clusters ", paste(kn, collapse = ",")) else (if (ko) "번호 없음" else "unnumbered"), fmt(nn)) } }, "")
    graphics::legend("left", legend = items, pch = 15, col = stc[keep_s], pt.cex = 3.6, cex = 1.7, bty = "n", y.intersp = 1.45, title = if (ko) "최고 단계 (군집 · 세포)" else "top stage (clusters · cells)", title.adj = 0)
    for (k in items) note("marker_legend", "item", k)
    main <- if (full) sprintf(if (ko) "GSE253355 골수 전체 %s 세포 — 저자 주석과 군집별 최고 마커 단계" else "GSE253355 marrow, all %s cells - author annotation vs top marker stage per cluster", fmt(ns))
            else sprintf(if (ko) "HSC 부터 Neutrophil · Monocyte 까지 저자 라벨 %d 개 (%s 세포) — 같은 UMAP" else "Author labels from HSC to Neutrophil / Monocyte (%d labels, %s cells) - same UMAP", length(subset_labels), fmt(ns))
    graphics::mtext(main, side = 3, outer = TRUE, line = 4.5, cex = 2.5, font = 2); graphics::mtext(sub, side = 3, outer = TRUE, line = 2.0, cex = 1.4, col = "grey20")
    note("figure", "title", main); note("figure", "sub", sub)
    foot <- if (ko) c("오른쪽 색 = 군집 평균 점수 (패널 10 유전자 log1p(CP10k) 합 / 10, gse253355_02 그대로) 가 가장 높은 단계. 후보일 뿐 확률 · 확정 세포 유형 · 분화 판정이 아니다.",
                      if (full) "왼쪽 색 = 저자 주석 (R73 과 같은 색). 좌표 = 저장된 R uwot UMAP 그대로 (재계산 없음). 연결선 · 화살표 없음."
                      else "선택한 세포를 모두 같은 전체 좌표에 그렸다 (재학습 · 재군집 · 재점수 없음, 점수 = 전체 자료 군집 점수). 연결선 · 화살표 없음.")
            else c("Right colour = stage with the highest cluster-mean score (sum of 10 panel genes' log1p(CP10k) / 10, as in gse253355_02). A candidate, not a probability, cell type or differentiation call.",
                   if (full) "Left colour = author annotation (same colours as R73). Coordinates = stored R uwot UMAP (not recomputed). No lines or arrows."
                   else "All selected cells at their full-UMAP coordinates (no re-embedding, re-clustering or re-scoring; full-data cluster scores). No lines or arrows.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.9 + (k - 1) * 1.75, cex = 1.35, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  list(text = nf$table(), coords = Us, colours = colours, numbered = numbered, unnumbered = data.frame(cluster = unn, n_cells = as.integer(cnt[unn]), stringsAsFactors = FALSE))
}

#' 그림 3. 군집 x 6 단계 원점수 heatmap (선형 한 척도, 칸 숫자), ● = 군집 안 1 위 (winner), 굵은 테 = 단계별 전체 군집 최고 (peak). 행 = 군집 번호 순서 · 저자 최다 라벨 · 세포 수
gsm_draw_heatmap <- function(path, S, cand, peaks, summ, sub = "", ko = TRUE, fam = "") {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; cl <- rownames(S); st <- colnames(S); nr <- length(cl); sp <- function(x) { z <- strsplit(x, "|", fixed = TRUE)[[1]]; z[nzchar(z)] }
  vmax <- if (any(S > 0, na.rm = TRUE)) max(S, na.rm = TRUE) else 1; miss <- if (ko) "결측" else "miss."
  colv <- function(v) ifelse(is.na(v), "white", ifelse(v <= 0, NSO_ZERO_COL, NSO_SEQ_COLS[pmax(1L, pmin(100L, ceiling(v / vmax * 100)))]))
  fv <- function(v) if (is.na(v)) miss else if (v == 0) "0" else formatC(v, format = "f", digits = 2)
  ci <- match(cl, cand$cluster); win <- lapply(cand$score_top_candidate[ci], sp); pkl <- stats::setNames(lapply(st, function(s) { z <- peaks$peak_clusters[peaks$stage == s]; if (length(z)) sp(z) else character(0) }), st)
  mw <- do.call(rbind, lapply(seq_len(nr), function(r) if (length(win[[r]])) data.frame(cluster = cl[r], stage = win[[r]], mark = "winner", stringsAsFactors = FALSE)))
  mp <- do.call(rbind, lapply(st, function(s) if (length(pkl[[s]])) data.frame(cluster = pkl[[s]], stage = s, mark = "peak", stringsAsFactors = FALSE)))
  si <- match(cl, summ$cluster)
  grDevices::png(path, width = 2900, height = 2800, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::par(mar = c(9.5, 30, 8.6, 17), oma = c(3.8, 0, 5.4, 0), family = fam)
    graphics::plot(NA, xlim = c(0.5, length(st) + 0.5), ylim = c(nr + 0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")
    for (r in seq_len(nr)) for (k in seq_along(st)) { v <- S[cl[r], st[k]]; isw <- st[k] %in% win[[r]]; tc <- if (!is.na(v) && v / vmax > 0.62) "white" else "grey5"
      graphics::rect(k - 0.5, r - 0.5, k + 0.5, r + 0.5, col = colv(v), border = "white", lwd = 2)
      graphics::text(k + if (isw) 0.08 else 0, r, fv(v), cex = 1.2, font = if (isw) 2 else 1, col = tc); if (isw) graphics::points(k - 0.3, r, pch = 16, cex = 1.6, col = tc) }
    if (!is.null(mp)) for (i in seq_len(nrow(mp))) { k <- match(mp$stage[i], st); r <- match(mp$cluster[i], cl); if (!is.na(r)) graphics::rect(k - 0.46, r - 0.42, k + 0.46, r + 0.42, border = "black", lwd = 4.5) }
    rl <- sprintf("%s   %s %d%% · n %s", cl, summ$author_top_label[si], round(100 * summ$author_top_frac[si]), .nso_fmt_n(summ$n_cells[si]))
    graphics::axis(2, at = seq_len(nr), labels = rl, las = 1, tick = FALSE, cex.axis = 1.18, line = -0.6); for (k in rl) note("heatmap", "row", k)
    gap <- cand$score_top_minus_second[ci]; dif <- summ$donor_equal_agreement[si] == "different_top"
    wl <- ifelse(cand$score_top_status[ci] == "score_top_candidate", sprintf("%s  +%s%s", cand$score_top_candidate[ci], .nso_fmt(gap, ko, 2), ifelse(dif, " *", "")),
                 ifelse(cand$score_top_status[ci] == "tied", gsub("|", "/", cand$score_top_candidate[ci], fixed = TRUE), if (ko) "0 점" else "zero"))
    wc <- ifelse(cand$score_top_status[ci] == "score_top_candidate", ATLAS_CL_STAGE_COL[cand$score_top_candidate[ci]], "grey30")
    graphics::text(length(st) + 0.62, seq_len(nr), wl, adj = c(0, 0.5), cex = 1.18, font = 2, col = wc, xpd = NA); for (k in wl) note("heatmap", "winner", k)
    wh <- if (ko) "1 위  (2 위와 차)" else "top  (gap to 2nd)"; graphics::text(length(st) + 0.62, 0.5 - 0.9, wh, adj = c(0, 0.5), cex = 1.15, xpd = NA); note("heatmap", "winner_header", wh)
    graphics::axis(3, at = seq_along(st), labels = st, tick = FALSE, cex.axis = 1.4, font = 2, line = 1.7); for (k in st) note("heatmap", "column", k)
    ph <- vapply(st, function(s) if (length(pkl[[s]])) sprintf(if (ko) "최고 군집 %s" else "peak %s", paste(pkl[[s]], collapse = ", ")) else (if (ko) "최고 없음" else "no peak"), "")
    graphics::mtext(ph, side = 3, at = seq_along(st), line = 0.5, cex = 1.05); for (k in ph) note("heatmap", "peak_header", k)
    xs <- seq(0.5, length(st) + 0.5, length.out = 101); yb <- nr + 0.75; graphics::rect(xs[-101], yb, xs[-1], yb + 0.4, col = NSO_SEQ_COLS, border = NA, xpd = NA)
    graphics::rect(0.5, yb, length(st) + 0.5, yb + 0.4, border = "grey40", xpd = NA); graphics::rect(0.05, yb, 0.35, yb + 0.4, col = NSO_ZERO_COL, border = "grey40", xpd = NA)
    graphics::text(0.2, yb + 0.85, "0", xpd = NA, cex = 1.15); tv <- c(0, vmax / 2, vmax); graphics::text(0.5 + length(st) * tv / vmax, yb + 0.85, formatC(tv, format = "f", digits = 2), xpd = NA, cex = 1.15)
    cb <- sprintf(if (ko) "색 = 원점수 (선형, 모든 칸 같은 척도 0 – %s)" else "colour = raw score (linear, one scale 0 - %s)", formatC(vmax, format = "f", digits = 2))
    graphics::text(0.5 + length(st) / 2, yb + 1.45, cb, xpd = NA, cex = 1.15); note("colorbar", "scale", cb)
    main <- sprintf(if (ko) "군집 x 6 단계 RNA 패널 점수 (군집 %d 개 전부)" else "Cluster x 6-stage RNA panel scores (all %d clusters)", nr)
    lgd <- if (ko) "● = 군집 안 1 위 (winner) · 굵은 테 = 그 단계 점수가 모든 군집 중 최고 (peak) · * = donor 동일 가중 1 위가 다름"
           else "● = within-cluster top (winner) · thick box = highest cluster for that stage (peak) · * = donor-equal top differs"
    graphics::mtext(main, side = 3, outer = TRUE, line = 2.9, cex = 2.0, font = 2); graphics::mtext(lgd, side = 3, outer = TRUE, line = 1.1, cex = 1.15)
    note("figure", "title", main); note("figure", "legend", lgd)
    foot <- if (ko) c("점수 = 군집 안 모든 세포의 패널 10 유전자 log1p(CP10k) 합 / 10 의 평균 (gse253355_02 그대로). 행 = 군집 번호 · 저자 최다 라벨과 비율 (사후, 정답 아님) · 세포 수.",
                      paste0(sub, " · 유전자별 기여 = cluster_stage_gene_contributions.csv"))
            else c("Score = mean over the cluster's cells of the summed log1p(CP10k) of 10 panel genes / 10 (as in gse253355_02). Row = cluster · most frequent author label and share (post hoc, not truth) · cells.",
                   paste0(sub, " · per-gene contributions = cluster_stage_gene_contributions.csv"))
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.6 + (k - 1) * 1.45, cex = 1.08, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  list(text = nf$table(), marks = rbind(mw, mp))
}


# ── 3. 실행기 ────────────────────────────────────────────────────────────────
gsm_main <- function(argv, root = NULL, pins = GSM_PINNED) {
  keys <- c("--umap" = "umap", "--scores-dir" = "scores_dir", "--object" = "object", "--panels" = "panels", "--hgnc" = "hgnc", "--out-dir" = "out", "--public-dir" = "public",
            "--focus-clusters" = "focus", "--author-col" = "author_col", "--assay" = "assay"); a <- list(); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (is.null(a$out) || is.na(a$out) || !nzchar(a$out)) { cat("--out-dir 가 필요하다.\n"); return(2L) }
  if (is.null(root)) root <- mdu_find_root()
  dflt <- list(umap = GSM_UMAP, scores_dir = GSM_SCORES_DIR, object = GSM_OBJECT, panels = GSM_PANELS, hgnc = GSM_HGNC, public = GSM_PUBLIC)
  for (k in names(dflt)) if (is.null(a[[k]])) a[[k]] <- file.path(root, dflt[[k]])
  a$focus <- if (is.null(a$focus)) GSM_FOCUS else trimws(strsplit(a$focus, ";", fixed = TRUE)[[1]])
  if (is.null(a$author_col)) a$author_col <- "cluster_anno_l2"; if (is.null(a$assay)) a$assay <- "RNA"
  for (d in c(a$out, a$public)) if (dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", d, "\n"); return(6L) }
  out <- a$out; dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", GSM_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(GSM_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = GSM_SOURCES)
    get("gsm_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

gsm_build <- function(a, out, snap, log, argv, pins) {
  ck <- .mbvq_checks(); t0 <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"); rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  files <- character(0); crt <- logical(0); fmt <- .nso_fmt_n
  wc <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  rdc <- function(p) utils::read.csv(p, colClasses = "character", na.strings = character(0), check.names = FALSE); num <- function(x) suppressWarnings(as.numeric(x))
  sd <- a$scores_dir; f_ls <- file.path(sd, "lossless.rds"); f_cm <- file.path(sd, "calc_manifest.json"); f_man <- file.path(sd, "manifest.json"); f_in <- file.path(sd, "input_sha256.csv")
  f_ac <- file.path(dirname(a$umap), "annotation_counts.csv")
  # ── (1) 입력 고정 ──
  if (grepl("\\.part$", a$object)) stop("다운로드 중인 .part 입력은 쓰지 않는다: ", a$object)
  for (p in c(a$umap, f_ls, f_cm, f_man, f_in, f_ac, a$object, a$panels, a$hgnc)) if (!file.exists(p)) stop("입력 파일이 없다: ", p)
  for (k in c("umap", "lossless", "object", "calc_manifest")) if (is.null(pins[[k]])) stop("hash 고정값이 없다: ", k)
  h0 <- c(umap = .nrv_sha(a$umap), lossless = .nrv_sha(f_ls), object = .nrv_sha(a$object), calc_manifest = .nrv_sha(f_cm), panels = .nrv_sha(a$panels), hgnc = .nrv_sha(a$hgnc))
  for (k in c("umap", "lossless", "object", "calc_manifest")) if (!identical(h0[[k]], pins[[k]])) stop(k, " hash 가 고정값과 다르다: ", h0[[k]])
  cm <- jsonlite::fromJSON(f_cm)$files; man <- jsonlite::fromJSON(f_man)$files; isha <- rdc(f_in); hx <- vapply(GSM_FIXED_FILES, function(f) .nrv_sha(file.path(sd, f)), "")
  bad <- GSM_FIXED_FILES[!vapply(GSM_FIXED_FILES, function(f) identical(hx[[f]], cm[[f]]), TRUE)]; if (length(bad)) stop("고정 점수 파일 hash 가 calc_manifest.json 기록과 다르다: ", paste(bad, collapse = ", "))
  hpc <- .nrv_sha(file.path(sd, "posthoc_author_composition.csv")); if (!identical(hpc, man[["posthoc_author_composition.csv"]])) stop("posthoc_author_composition.csv hash 가 manifest.json 기록과 다르다.")
  rec <- function(k) isha$sha256[isha$input == k]
  if (!identical(rec("panels"), h0[["panels"]]) || !identical(rec("hgnc"), h0[["hgnc"]]) || !identical(rec("input"), h0[["object"]])) stop("패널 · HGNC · 객체 hash 가 gse253355_02 input_sha256.csv 기록과 다르다.")
  ck$add("inputs_pinned", TRUE, sprintf("umap %s, lossless %s, object %s, calc_manifest %s == pinned; %s == calc_manifest.json; panels %s, hgnc %s, object == gse253355_02 input_sha256.csv",
                                        substr(h0[["umap"]], 1, 12), substr(h0[["lossless"]], 1, 12), substr(h0[["object"]], 1, 12), substr(h0[["calc_manifest"]], 1, 12),
                                        paste(GSM_FIXED_FILES, collapse = ", "), substr(h0[["panels"]], 1, 12), substr(h0[["hgnc"]], 1, 12)))
  # ── (2) 고정 좌표 · 점수 ──
  um <- readRDS(a$umap); ls <- readRDS(f_ls); cells <- ls$cells; cl <- as.character(ls$clusters); n <- length(cells)
  miss_u <- setdiff(c("cells", "U", "new_cluster", "author_label", "subset_labels", "subset_mask"), names(um)); if (length(miss_u)) stop("좌표 RDS 에 항목이 없다: ", paste(miss_u, collapse = ", "))
  U <- um$U; mask <- um$subset_mask; au <- as.character(um$author_label); sl <- as.character(um$subset_labels)
  ok_u <- identical(um$cells, cells) && identical(as.character(um$new_cluster), cl) && is.matrix(U) && identical(dim(U), c(n, 2L)) && all(is.finite(U)) && identical(rownames(U), cells) &&
          is.logical(mask) && length(mask) == n && identical(mask, au %in% sl)
  if (!ok_u) stop("좌표 RDS 의 세포 · 군집 · 좌표 · 부분집합이 lossless 와 맞지 않는다.")
  clus <- rdc(file.path(sd, "cluster_scores.csv")); for (k in c("n_cells", "n_donors", "n_donors_min_cells")) clus[[k]] <- as.integer(clus[[k]])
  for (k in c("mean_all_cells", "mean_donor_equal", "frac_cells_ge1", "frac_cells_ge2", "mean_available_genes")) if (k %in% names(clus)) clus[[k]] <- num(clus[[k]])
  cand <- rdc(file.path(sd, "cluster_candidates.csv")); peaks <- rdc(file.path(sd, "stage_peak_clusters.csv"))
  for (k in c(paste0("score_", NSO_STAGES), "score_top_value", "score_second_value", "score_top_minus_second")) cand[[k]] <- num(cand[[k]])
  for (k in c("n_cells", "n_donors")) cand[[k]] <- as.integer(cand[[k]])
  dsc <- if (length(clus$mean_all_cells) == length(ls$scores$mean_all_cells)) max(abs(clus$mean_all_cells - ls$scores$mean_all_cells)) else Inf   # CSV 17 자리 → R 소수 변환은 1 ulp 어긋날 수 있다 (Ruling M7)
  ok_sc <- identical(clus$cluster, ls$scores$cluster) && identical(clus$stage, ls$scores$stage) && !is.na(dsc) && dsc <= 1e-14 &&
           identical(cand$cluster, ls$candidates$cluster) && identical(cand$score_top_candidate, ls$candidates$score_top_candidate) &&
           identical(cand$score_top_status, ls$candidates$score_top_status) && identical(peaks$peak_clusters, ls$peaks$peak_clusters) && all(cl %in% cand$cluster)
  if (!ok_sc) stop("고정 점수 표가 lossless.rds 사본과 다르다 (최대 차 ", format(dsc), ").")
  ck$add("frozen_scores_equal_lossless", TRUE, sprintf("cluster_scores.csv mean_all_cells (%d rows) == lossless.rds copy within decimal round-trip (max |diff| %.3g <= 1e-14); winners, statuses and stage peaks identical; %d clusters x %d stages used as stored",
                                                       nrow(clus), dsc, nrow(cand), length(NSO_STAGES)))
  log("frozen inputs:", n, "cells,", nrow(cand), "clusters | subset", sum(mask), "cells in", length(sl), "labels")
  # ── (3) 원 count → 유전자 값 (고정 점수 재현 확인, Ruling M6) ──
  ro <- gsx_open_rds(a$object); obj <- ro$object; md <- gsx_meta(obj); cn <- gsx_counts(obj, a$assay); X <- cn$counts; rm(obj, ro); invisible(gc())
  if (!a$author_col %in% names(md)) stop("meta.data 에 저자 라벨 열이 없다: ", a$author_col)
  if (!identical(colnames(X), cells) || !identical(rownames(md), cells) || !identical(as.character(md[[a$author_col]]), au)) stop("객체 세포 ID · 저자 라벨이 좌표 RDS · lossless 와 다르다.")
  ck$add("cells_clusters_coordinates_unchanged", TRUE, sprintf("%d cell IDs identical in order (UMAP RDS, lossless, raw counts, meta.data); new_cluster == lossless clusters; author_label == meta.data %s; U %d x 2 finite, drawn as stored",
                                                               n, a$author_col, n))
  pan <- rdc(a$panels); kc <- c("stage", "rank", "ensg", "symbol"); pan$rank <- as.integer(pan$rank); pan <- pan[order(match(pan$stage, NSO_STAGES), pan$rank), kc]
  hg <- utils::read.delim(a$hgnc, quote = "", colClasses = "character", stringsAsFactors = FALSE)[, c("symbol", "prev_symbol", "alias_symbol")]
  pm <- gsx_panel_map(rownames(X), pan, hg); cs <- gsx_cell_scores(X, pm); rm(X); invisible(gc())
  if (!identical(pm, ls$panel_map) || !identical(cs$fixed, ls$fixed) || !identical(cs$library, ls$library)) stop("다시 계산한 패널 대응 · 세포 점수 · library 가 lossless 와 identical 하지 않다.")
  ck$add("gene_values_reproduce_frozen_scores", TRUE, sprintf("panel map (%d rows), per-cell stage scores (%d x %d) and libraries recomputed from raw counts (R69, log1p(1e4*count/all-gene library)) are identical to lossless.rds",
                                                              nrow(pm), nrow(cs$fixed), ncol(cs$fixed)))
  log("gene values recomputed:", nrow(cs$gene_log), "measured panel genes; frozen scores reproduced")
  # ── (4) 유전자 기여 표 (계산 고정) ──
  gt <- gsm_gene_table(cs$gene_log, cs$gene_det, pm, cl, clus); cc <- gsm_contribution_check(gt, clus, GSM_CONTRIB_TOL)
  ck$add("contribution_sum_equals_cluster_score", all(cc$pass), sprintf("sum of gene contributions == cluster_scores.csv mean_all_cells for %d cluster x stage; max |diff| %.3g <= %g; denominator = listed genes per stage (%s)",
                                                                        nrow(cc), max(cc$abs_diff), GSM_CONTRIB_TOL, paste(unique(gt$n_listed_in_stage), collapse = ",")))
  mg <- gt[gt$gene_status == "missing", , drop = FALSE]
  ok_m <- all(is.na(mg$mean_log1p_cp10k)) && all(is.na(mg$frac_cells_detected)) && all(mg$contribution == 0) && !anyNA(gt$mean_log1p_cp10k[gt$gene_status == "measured"])
  ck$add("missing_gene_na_and_zero", ok_m, sprintf("missing panel genes (%s): expression and detection NA, contribution 0 in %d rows; measured genes all have values",
                                                   if (nrow(mg)) paste(unique(mg$gene), collapse = ", ") else "none", nrow(mg)))
  dup <- vapply(pm$symbol, function(g) sum(pm$symbol == g), 0L, USE.NAMES = FALSE)
  pd <- data.frame(pm, listed_in_stages = vapply(pm$symbol, function(g) paste(NSO_STAGES[NSO_STAGES %in% pm$stage[pm$symbol == g]], collapse = "|"), "", USE.NAMES = FALSE),
                   n_stages_listed = dup, stringsAsFactors = FALSE)
  wc(gt, "cluster_stage_gene_contributions.csv"); wc(cc, "contribution_check.csv"); wc(pd, "panel_gene_diagnostics.csv")
  cmf <- file.path(out, "calc_manifest.json"); calc_files <- files
  writeLines(jsonlite::toJSON(list(phase = "gene contribution tables fixed before author labels are used", written_at = stamp(), inputs = as.list(h0),
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  log("calculation fixed:", length(calc_files), "files | contribution max |diff|", signif(max(cc$abs_diff), 3))
  # ── (5) 사후: 저자 구성 · 군집 요약 · mismatch 유전자 ──
  comp <- gsm_author_composition(cl, au, cand$cluster); ex <- rdc(file.path(sd, "posthoc_author_composition.csv"))
  ok_c <- identical(ex$cluster, comp$cluster) && identical(ex$author_label, comp$author_label) && identical(as.integer(ex$n_cells), comp$n_cells) &&
          isTRUE(max(abs(num(ex$frac_of_cluster) - comp$frac_of_cluster)) <= 1e-14)                     # 비율은 CSV 소수 변환 허용 (Ruling M7)
  ck$add("author_composition_matches_existing", ok_c, sprintf("cluster x author composition from the coordinate RDS (%d rows): clusters, labels and counts identical to gse253355_02 posthoc_author_composition.csv; fractions within 1e-14", nrow(comp)))
  comp$note <- rep("post hoc display and comparison only; not used for any score or candidate; not truth or accuracy", nrow(comp)); wc(comp, "posthoc_author_composition.csv")
  subn <- stats::setNames(as.integer(table(factor(cl[mask], levels = cand$cluster))), cand$cluster)
  summ <- gsm_cluster_summary(cand, peaks, comp, subn); wc(summ, "cluster_marker_summary.csv")
  fg <- gsm_focus_genes(gt, summ, a$focus, GSM_TOP_GENES); wc(fg, "focus_cluster_top_genes.csv")
  # ── (6) 그림 (R base graphics) ──
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; acols <- .gsu_cols(au, sl); fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  ac <- rdc(f_ac); ok_ac <- identical(sort(ac$author_label), sort(names(acols))) && identical(unname(acols[ac$author_label]), ac$colour)
  ck$add("author_colours_match_r73", ok_ac, sprintf("%d author-label colours identical to %s/annotation_counts.csv (R73 mapping)", length(acols), basename(dirname(a$umap))))
  cdesc <- if (!is.null(um$uwot_version)) sprintf("uwot %s · seed %s", um$uwot_version, um$seed) else "stored coordinates"
  sub <- sprintf(if (ko) "R %s base graphics · 좌표 = %s (%s, 재계산 없음) · 점수 = %s (군집 %d 개 x %d 단계)" else "R %s base graphics · coordinates = %s (%s, not recomputed) · scores = %s (%d clusters x %d stages)",
                 as.character(getRversion()), basename(dirname(a$umap)), cdesc, basename(sd), nrow(cand), length(NSO_STAGES))
  f1 <- file.path(fig, "fig1_umap_full_author_vs_marker.png"); r1 <- gsm_draw_pair(f1, U, au, acols, cl, cand, sub = sub, ko = ko, fam = fam)
  f2 <- file.path(fig, "fig2_umap_lineage_author_vs_marker.png")
  r2 <- gsm_draw_pair(f2, U, au, acols, cl, cand, mask = mask, subset_labels = sl, label_min_frac = GSM_LABEL_MIN_FRAC, sub = sub, ko = ko, fam = fam)
  S <- .gsx_wide(clus, "mean_all_cells"); f3 <- file.path(fig, "fig3_cluster_stage_heatmap.png"); r3 <- gsm_draw_heatmap(f3, S, cand, peaks, summ, sub = sub, ko = ko, fam = fam)
  files <- c(files, f1, f2, f3); wc(rbind(r1$text, r2$text, r3$text), "figure_text.csv"); log("figures written")
  colok <- function(r) { e <- ifelse(r$colours$status == "score_top_candidate", unname(ATLAS_CL_STAGE_COL[r$colours$winner]), ifelse(r$colours$status == "tied", ATLAS_CL_TIE_COL, ATLAS_CL_NONE_COL))
    all(r$colours$n_colours == 1L) && identical(unname(e), r$colours$colour) }
  ck$add("figure_colours_match_winner", colok(r1) && colok(r2) && identical(r1$colours$cluster, cand$cluster),
         sprintf("every cell drawn in the colour of its cluster's cluster_candidates.csv winner (one colour per cluster; %d clusters full, %d subset)", nrow(r1$colours), nrow(r2$colours)))
  same_xy <- function(A, B) identical(as.vector(A), as.vector(B)) && identical(dim(A), dim(B)) && identical(dimnames(A), dimnames(B))   # 값 · 차원 · 세포 ID (uwot 의 scaled:center 속성은 행 선택에서 빠진다, Ruling M8)
  ck$add("subset_same_coordinates", same_xy(r1$coords, U) && same_xy(r2$coords, U[mask, , drop = FALSE]) && nrow(r2$coords) == sum(mask) && identical(rownames(r2$coords), cells[mask]),
         sprintf("figure 1 draws all %d cells and figure 2 all %d subset cells at the stored UMAP coordinates (values, dimensions and cell IDs identical; no re-embedding, re-clustering or re-scoring); figure 2 numbers %d clusters, %d clusters (%d cells) colour only",
                 n, sum(mask), length(r2$numbered), nrow(r2$unnumbered), sum(r2$unnumbered$n_cells)))
  # ── (7) 기록 · 요약 · 공개 ──
  allinp <- c(umap = a$umap, annotation_counts = f_ac, lossless = f_ls, calc_manifest = f_cm, gse253355_02_manifest = f_man, gse253355_02_input_sha256 = f_in,
              cluster_scores = file.path(sd, "cluster_scores.csv"), cluster_candidates = file.path(sd, "cluster_candidates.csv"), stage_peak_clusters = file.path(sd, "stage_peak_clusters.csv"),
              posthoc_author_composition = file.path(sd, "posthoc_author_composition.csv"), object = a$object, panels = a$panels, hgnc = a$hgnc)
  h1 <- vapply(names(allinp), function(k) .nrv_sha(allinp[[k]]), "")
  keep_ok <- identical(h1[["umap"]], h0[["umap"]]) && identical(h1[["lossless"]], h0[["lossless"]]) && identical(h1[["object"]], h0[["object"]]) && identical(h1[["calc_manifest"]], h0[["calc_manifest"]]) &&
             identical(h1[["panels"]], h0[["panels"]]) && identical(h1[["hgnc"]], h0[["hgnc"]]) && identical(unname(h1[c("cluster_scores", "cluster_candidates", "stage_peak_clusters")]), unname(hx[1:3])) &&
             identical(h1[["posthoc_author_composition"]], hpc)
  ck$add("existing_inputs_unchanged", keep_ok, sprintf("%d existing inputs re-hashed after the run: unchanged (coordinates, lossless, object, calc_manifest, panels, hgnc, frozen score tables)", length(allinp)))
  wc(data.frame(input = names(allinp), path = unname(allinp), sha256 = unname(h1), stringsAsFactors = FALSE), "input_sha256.csv"); wc(snap, "source_snapshot.csv")
  pkg <- c("Matrix", "digest", "jsonlite")
  wc(data.frame(item = c("R_executable", "R", "platform", "libPaths", "packages", "graphics", "title_font", "command", "working_dir", "started_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pkg, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), "R base graphics (grDevices::png, res 100)",
                          if (ko) MBVR_TITLE_FAMILY else "default", paste(c("74_gse253355_marker_umap.R", argv), collapse = " "), getwd(), format(t0, "%Y-%m-%dT%H:%M:%S%z"), stamp()),
                stringsAsFactors = FALSE), "environment.csv")
  # 요약 (25 줄 이하)
  atop <- stats::setNames(summ$author_top_label, summ$cluster); tbw <- table(factor(summ$winner_stage[summ$winner_status == "score_top_candidate"], levels = NSO_STAGES))
  pkt <- paste(vapply(seq_len(nrow(peaks)), function(i) { z <- strsplit(peaks$peak_clusters[i], "|", fixed = TRUE)[[1]]
    if (!length(z) || !nzchar(peaks$peak_clusters[i])) sprintf("%s -", peaks$stage[i]) else sprintf("%s %s", peaks$stage[i], paste(sprintf("%s (%s)", z, atop[z]), collapse = "/")) }, ""), collapse = " · ")
  pnw <- summ[nzchar(summ$stage_peak_of) & !summ$winner_is_stage_peak, , drop = FALSE]
  fl <- vapply(a$focus, function(k) { z <- fg[fg$cluster == k & grepl("winner", fg$role), , drop = FALSE]; z <- utils::head(z, 3L)
    sprintf("  - %s (%s %d%%) 1 위 %s %.2f: %s", k, atop[[k]], round(100 * summ$author_top_frac[summ$cluster == k]), summ$winner_stage[summ$cluster == k], summ$winner_score[summ$cluster == k],
            paste(sprintf("%s %.0f%%", z$gene, 100 * z$share_of_stage_score), collapse = " · ")) }, "")
  sm <- c("# GSE253355 실제 R UMAP 에 고정 6 단계 RNA 패널 대응 (표시 · 집계 전용) — 요약", "",
          sprintf("- 그림 3 장은 모두 R (%s, base graphics). 좌표 = %s 의 R uwot UMAP 그대로, 점수 = %s 군집 점수 그대로 (재계산 없음).", R.version.string, basename(dirname(a$umap)), basename(sd)),
          sprintf("- 입력 sha256: 좌표 %s… · lossless %s… · 객체 %s… · calc_manifest %s… (고정값 일치, 실행 뒤에도 같음).", substr(h0[["umap"]], 1, 12), substr(h0[["lossless"]], 1, 12),
                  substr(h0[["object"]], 1, 12), substr(h0[["calc_manifest"]], 1, 12)),
          sprintf("- 군집 %d 개의 1 위 단계 (winner): %s · 동률 %d · 0 점 %d.", nrow(summ), paste(sprintf("%s %d", NSO_STAGES, as.integer(tbw)), collapse = " · "),
                  sum(summ$winner_status == "tied"), sum(summ$winner_status == "unassigned_zero")),
          paste0("- 단계별 전체 군집 최고 (peak, 저자 최다 라벨): ", pkt, "."),
          paste0("- peak 이지만 그 군집 1 위는 다른 단계: ", if (nrow(pnw)) paste(sprintf("군집 %s (peak %s, 1 위 %s)", pnw$cluster, pnw$stage_peak_of, pnw$winner_stage), collapse = " · ") else "없음", "."),
          "- 주요 mismatch 군집의 1 위 단계 기여 상위 유전자 (단계 점수 대비 비율):", fl,
          sprintf("- 부분집합 %s 세포 (저자 라벨 %d 개): 번호 군집 %s; 번호 없는 군집 %d 개 %s 세포는 색만.", fmt(sum(mask)), length(sl), paste(r2$numbered, collapse = " · "), nrow(r2$unnumbered),
                  fmt(sum(r2$unnumbered$n_cells))),
          sprintf("- 유전자 기여 = 군집 평균 log1p(CP10k) / %s, 기여 합 = 군집 점수 (최대 차 %.2g, %d 군집 x 단계). 없는 유전자 %s 은 발현 NA · 기여 0.", paste(unique(gt$n_listed_in_stage), collapse = ","),
                  max(cc$abs_diff), nrow(cc), if (nrow(mg)) paste(unique(mg$gene), collapse = ", ") else "없음"),
          paste0("- 두 단계에 겹치는 패널 유전자: ", if (any(dup > 1L)) paste(unique(sprintf("%s (%s)", pd$symbol[dup > 1L], gsub("|", " · ", pd$listed_in_stages[dup > 1L], fixed = TRUE))), collapse = ", ") else "없음",
                 " — 두 단계 점수에 모두 들어간다."),
          "- 저자 라벨은 계산 고정 뒤 표시 · 사후 대조에만 썼다. 후보는 확률 · 확정 세포 유형 · 분화 판정이 아니며 정확도를 주장하지 않는다.",
          "- 표: cluster_marker_summary · cluster_stage_gene_contributions · focus_cluster_top_genes · panel_gene_diagnostics · contribution_check. 그림: figures/fig1 · fig2 · fig3 (PNG).")
  pub_list <- c(GSM_PUBLIC_FILES, file.path("figures", basename(c(f1, f2, f3))))
  ck$add("tables_roundtrip_ok", all(crt), sprintf("%d CSVs re-read with identical dimensions, names and cell text", length(crt)))
  ck$add("public_copy_aggregate_only", !any(grepl("lossless|coordinates|\\.rds$|\\.gz$", pub_list)), sprintf("%d aggregate tables, summary and 3 PNG figures for the public copy; no per-cell or coordinate files", length(GSM_PUBLIC_FILES)))
  cks <- ck$table(); sm <- c(sm, sprintf("- 검사 %d / %d 통과: %s.", sum(cks$pass), nrow(cks), paste(cks$check, collapse = ", ")))
  p0 <- file.path(out, "summary_ko.md"); writeLines(sm, p0); files <- c(files, p0); wc(cks, "checks.csv"); allf <- c(files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), inputs = as.list(h0), files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  pub <- a$public; dir.create(file.path(pub, "figures"), recursive = TRUE, showWarnings = FALSE); src <- file.path(out, pub_list); dest <- file.path(pub, pub_list)
  if (!all(file.copy(src, dest, overwrite = FALSE))) stop("공개 사본 복사 실패")
  writeLines(jsonlite::toJSON(list(source_output = normalizePath(out), inputs = as.list(h0), files = as.list(stats::setNames(vapply(dest, .nrv_sha, ""), pub_list))), auto_unbox = TRUE, pretty = TRUE),
             file.path(pub, "public_manifest.json"))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .gsm_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.gsm_bb, .f))
  quit(status = gsm_main(commandArgs(trailingOnly = TRUE)))
}
