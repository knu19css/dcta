# =============================================================================
# 67_atlas_cluster_view.R — 군집 단위 atlas 6 단계 패널 그림: 그림 A (군집 내부 1 위 후보, p92 · p51) · 그림 B (단계별 군집 평균 점수와 단계별 최고 군집)
# =============================================================================
#
# 계약: codex_cluster_instruction.md 5 항. 원 UMAP 은 표시 좌표일 뿐이다 — 새 UMAP · edge · 군집을 만들지 않고, 표시 위치를 계산에 쓰지 않는다.
#   그림 A: 군집 전체를 그 군집 1 위 후보 색으로 칠한다. 회색 = unassigned_zero / unavailable, 보라 = tied. 군집 번호는 대표 세포(UMAP 중앙값에 가장 가까운 실제 세포)에.
#   그림 B: 6 단계 패널, 군집 평균 점수(mean_log_pooled)를 그 군집 모든 세포에 같은 색으로, 6 패널 같은 척도. 단계별 최고 군집 = 진한 테두리 표식 + 군집 번호 · 세포 수.
#   그림 글자는 figure_text 로 남긴다. 결측은 "결측" / "miss." 로 쓰고 글자 "NA" 를 쓰지 않는다. R34 · R36 · R50 · R61 · R62 도우미를 수정 없이 재사용한다.
# =============================================================================

ATLAS_CL_STAGE_COL <- c(Myeloblast = "#0072B2", PM = "#56B4E9", MY = "#009E73", MM = "#E69F00", BN = "#D55E00", SN = "#CC79A7")   # 색각 이상에 덜 겹치는 6 색
ATLAS_CL_TIE_COL <- "#4B2C83"
ATLAS_CL_NONE_COL <- "grey84"

#' 군집 번호 자리: 입력 순서대로(호출자가 우선순위로 정렬) 대표 위치 중심에 먼저 놓고, 이미 놓은 번호와 겹치면 둘레 고리로 옮긴다.
#' 자리가 없으면 중심에 두고 free = FALSE 로 남긴다 (그림 글자 표에 기록).
atlas_cl_place_labels <- function(x, y, w, h, xlim, ylim) {
  n <- length(x); R <- matrix(NA_real_, n, 4); moved <- logical(n); free <- logical(n); ang <- seq(0, 2 * pi, length.out = 13)[-13]
  ok <- function(r) { if (r[1] < xlim[1] || r[2] > xlim[2] || r[3] < ylim[1] || r[4] > ylim[2]) return(FALSE)
    for (j in which(!is.na(R[, 1]))) if (!(r[2] <= R[j, 1] || R[j, 2] <= r[1] || r[4] <= R[j, 3] || R[j, 4] <= r[3])) return(FALSE)
    TRUE }
  for (i in seq_len(n)) { cand <- c(list(c(x[i], y[i])), unlist(lapply(c(1, 1.6, 2.4, 3.4, 4.8, 6.5), function(d) lapply(ang, function(a) c(x[i] + cos(a) * d * w[i], y[i] + sin(a) * d * h[i] * 1.3))),
                                                            recursive = FALSE))
    for (k in seq_along(cand)) { r <- c(cand[[k]][1] - w[i] / 2, cand[[k]][1] + w[i] / 2, cand[[k]][2] - h[i] / 2, cand[[k]][2] + h[i] / 2)
      if (ok(r)) { R[i, ] <- r; moved[i] <- k > 1L; free[i] <- TRUE; break } }
    if (!free[i]) R[i, ] <- c(x[i] - w[i] / 2, x[i] + w[i] / 2, y[i] - h[i] / 2, y[i] + h[i] / 2) }
  data.frame(x = x, y = y, lx = (R[, 1] + R[, 2]) / 2, ly = (R[, 3] + R[, 4]) / 2, x0 = R[, 1], x1 = R[, 2], y0 = R[, 3], y1 = R[, 4], moved = moved, free = free)
}

#' 군집마다 그림 A 색: 단일 1 위 = 단계 색, 동률 = 보라, 0 · 판정불가 = 회색
.acl_top_cols <- function(cand) unname(ifelse(cand$score_top_status == "score_top_candidate", ATLAS_CL_STAGE_COL[cand$score_top_candidate],
                                              ifelse(cand$score_top_status == "tied", ATLAS_CL_TIE_COL, ATLAS_CL_NONE_COL)))

#' 번호 상자 그리기 (옮긴 번호는 대표 위치까지 지시선)
.acl_draw_labels <- function(lp, text, cex) {
  for (i in which(lp$moved)) graphics::segments(lp$x[i], lp$y[i], lp$lx[i], lp$ly[i], col = "grey20", lwd = 1.1)
  graphics::rect(lp$x0, lp$y0, lp$x1, lp$y1, col = grDevices::adjustcolor("white", 0.86), border = "grey25", lwd = 1)
  graphics::text(lp$lx, lp$ly, text, cex = cex, font = 2, col = "grey5")
}

.acl_label_layout <- function(rp, cex, b, order_by) {
  span <- diff(b$xlim); txt <- rp$cluster; o <- order(-order_by, seq_along(order_by))
  w <- graphics::strwidth(txt, cex = cex) * 1.25 + 0.006 * span; h <- graphics::strheight(txt, cex = cex) * 1.7 + 0.004 * span
  usr <- graphics::par("usr"); lp <- atlas_cl_place_labels(rp$x[o], rp$y[o], w[o], h[o], usr[1:2], usr[3:4]); lp$text <- txt[o]; lp
}

#' 그림 A. 군집 내부 1 위 후보 (p92 · p51 나란히). 반환: 글자 표
atlas_cl_draw_top <- function(path, U, md, cand, reps, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U)
  grDevices::png(path, width = 3400, height = 2000, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(10, 1.25)); graphics::par(oma = c(3.6, 0, 4, 0), family = fam)
    for (p in c("p92", "p51")) { cl <- if (p == "p92") md$cluster92 else md$cluster51; z <- cand[cand$partition == p, , drop = FALSE]; col <- .acl_top_cols(z)
      graphics::par(mar = c(0.4, 0.4, 5.6, 0.4)); mbvu_open_panel(b)
      graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.18, col = col[match(cl, z$cluster)])
      rp <- reps[[p]]; rp <- rp[match(z$cluster, rp$cluster), , drop = FALSE]; cex <- if (p == "p92") 0.95 else 1.1
      lp <- .acl_label_layout(rp, cex, b, rp$n_cells); .acl_draw_labels(lp, lp$text, cex)
      if (any(!lp$free)) note(p, "labels_overlapping", paste(lp$text[!lp$free], collapse = " "))
      t1 <- table(factor(z$score_top_candidate[z$score_top_status == "score_top_candidate"], levels = NSO_STAGES))
      ttl <- sprintf(if (ko) "%s · 군집 %d 개" else "%s · %d clusters", p, nrow(z))
      sub1 <- paste(sprintf("%s %d", NSO_STAGES, as.integer(t1)), collapse = " · ")
      sub2 <- sprintf(if (ko) "동률 %d · 0(미검출) %d · 판정불가 %d · 작은 군집(경고) %d · 동일 block 가중에서 1 위가 다름 %d" else "tied %d · zero %d · unavailable %d · small (warned) %d · equal-block top differs %d",
                      sum(z$score_top_status == "tied"), sum(z$score_top_status == "unassigned_zero"), sum(z$score_top_status == "unavailable"), sum(z$small_cluster),
                      sum(z$equal_block_agreement == "different_top"))
      graphics::box(col = "grey60"); graphics::title(main = ttl, cex.main = 2.0, line = 3.6)
      graphics::mtext(paste0(if (ko) "1 위 후보 군집 수: " else "clusters by top candidate: ", sub1), side = 3, line = 1.9, cex = 1.12); graphics::mtext(sub2, side = 3, line = 0.5, cex = 1.05)
      note(p, "title", ttl); note(p, "sub1", sub1); note(p, "sub2", sub2) }
    graphics::par(mar = c(0.5, 2, 0.5, 2)); graphics::plot.new(); graphics::plot.window(c(0, 8), c(0, 1))
    labs <- c(NSO_STAGES, if (ko) "동률 (tied)" else "tied", if (ko) "0 · 판정불가" else "zero / unavailable"); cols <- c(ATLAS_CL_STAGE_COL, ATLAS_CL_TIE_COL, ATLAS_CL_NONE_COL)
    xs <- seq(0.15, 7.15, length.out = 8); graphics::rect(xs, 0.3, xs + 0.22, 0.75, col = cols, border = "grey30"); graphics::text(xs + 0.28, 0.52, labs, adj = c(0, 0.5), cex = 1.35)
    for (k in seq_along(labs)) note("legend", "item", labs[k])
    main <- if (ko) "군집마다 atlas 6 단계 패널 점수가 가장 높은 단계 (score_top_candidate) — 군집 평균, 같은 원 UMAP" else "Per cluster: atlas stage panel with the highest cluster-mean score (score_top_candidate) - same original UMAP"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.4, cex = 2.05, font = 2); note("figure", "title", main)
    foot <- if (ko) c("점수 = 군집 안 모든 세포의 패널 유전자 logRNA 평균 (mean_log_pooled). 군집 번호는 그 군집 세포 중 UMAP 중앙값에 가장 가까운 세포에 둔다 (표시 전용).",
                      "score_top_candidate 는 확률 · 확정 cell type · 분화 판정이 아니다. 패널 유전자는 단핵구 · T · 적혈구계 · 전구세포에서도 발현된다.")
            else c("Score = mean logRNA of panel genes over all cells of the cluster (mean_log_pooled). Cluster IDs sit at the cell nearest the cluster's UMAP median (display only).",
                   "score_top_candidate is not a probability, a confirmed cell type or a differentiation call. Panel genes are also expressed in monocytes, T, erythroid and progenitor cells.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.6 + (k - 1) * 1.35, cex = 1.12, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 그림 B. 단계별 군집 평균 점수 (6 패널, 같은 척도)와 단계별 최고 군집. S = 군집 x 단계 pooled 행렬, peaks = 그 분할의 stage_peak_clusters. 반환: 글자 표
atlas_cl_draw_stage <- function(path, U, cluster, S, peaks, rp, cand, partition, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U); cluster <- as.character(cluster); ci <- match(cluster, rownames(S))
  vmax <- if (any(S > 0, na.rm = TRUE)) max(S, na.rm = TRUE) else 1
  colv <- function(v) ifelse(is.na(v) | v <= 0, NSO_ZERO_COL, NSO_SEQ_COLS[pmax(1L, pmin(100L, ceiling(v / vmax * 100)))])
  grDevices::png(path, width = 3000, height = 2650, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(c(1:6, 7, 7, 7), 3, byrow = TRUE), heights = c(5, 5, 1.05)); graphics::par(oma = c(4.4, 0, 3.8, 0), family = fam)
    for (s in NSO_STAGES) { pk <- peaks[peaks$stage == s, , drop = FALSE]
      graphics::par(mar = c(0.3, 0.3, 6.2, 0.3)); mbvu_open_panel(b)
      graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.2, col = colv(S[ci, s]))
      ttl <- sprintf(if (ko) "%s 패널 · 군집 평균" else "%s panel · cluster mean", s)
      if (pk$peak_status %in% c("peak", "tied")) {
        pc <- strsplit(pk$peak_clusters, "|", fixed = TRUE)[[1]]; j <- match(pc, rp$cluster)
        graphics::points(rp$x[j], rp$y[j], pch = 21, cex = 3.6, lwd = 3.6, col = "black", bg = colv(S[pc, s]))
        lab <- sprintf(if (ko) "군집 %s · n %s" else "cluster %s · n %s", pc, .nso_fmt_n(rp$n_cells[j])); cx <- 1.05
        w <- graphics::strwidth(lab, cex = cx) * 1.15 + 0.01 * diff(b$xlim); h <- graphics::strheight(lab, cex = cx) * 1.8
        right <- rp$x[j] + 0.03 * diff(b$xlim) + w <= b$xlim[2]; lx <- ifelse(right, rp$x[j] + 0.03 * diff(b$xlim) + w / 2, rp$x[j] - 0.03 * diff(b$xlim) - w / 2)
        ly <- rp$y[j] + (seq_along(j) - 1) * 1.3 * h
        graphics::segments(rp$x[j], rp$y[j], lx - ifelse(right, w / 2, -w / 2), ly, lwd = 1.6)
        graphics::rect(lx - w / 2, ly - h / 2, lx + w / 2, ly + h / 2, col = grDevices::adjustcolor("white", 0.9), border = "black", lwd = 1.6)
        graphics::text(lx, ly, lab, cex = cx, font = 2); for (k in seq_along(lab)) note(s, "peak_label", lab[k])
        sub1 <- sprintf(if (ko) "단계별 최고 군집 %s (점수 %s) · 2 위 군집 %s (%s)" else "stage peak cluster %s (score %s) · runner-up %s (%s)", .nso_dots(pk$peak_clusters),
                        .nso_fmt(pk$peak_score, ko, 3), if (nzchar(pk$second_cluster)) pk$second_cluster else "-", .nso_fmt(pk$second_score, ko, 3))
        sub2 <- sprintf(if (ko) "그 군집의 군집 내부 1 위 후보: %s · 이 단계가 1 위 후보인 군집 %d 개" else "within-cluster top of the peak cluster: %s · clusters with this stage as top: %d",
                        gsub("|", " · ", pk$peak_cluster_top_candidate, fixed = TRUE), pk$n_clusters_with_stage_top)
      } else {
        sub1 <- if (pk$peak_status == "all_zero") (if (ko) "모든 군집 점수 0 (미검출)" else "all cluster scores are 0 (not detected)") else (if (ko) "판정불가 (패널 없음)" else "unavailable (no panel)")
        sub2 <- sprintf(if (ko) "이 단계가 1 위 후보인 군집 %d 개" else "clusters with this stage as top: %d", pk$n_clusters_with_stage_top)
      }
      graphics::box(col = "grey60"); graphics::title(main = ttl, cex.main = 1.95, line = 4.1); graphics::mtext(sub1, side = 3, line = 2.2, cex = 1.05); graphics::mtext(sub2, side = 3, line = 0.7, cex = 1.0)
      note(s, "title", ttl); note(s, "sub1", sub1); note(s, "sub2", sub2) }
    graphics::par(mar = c(1.6, 22, 1.2, 22)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    xs <- seq(0, 1, length.out = 101); graphics::rect(xs[-101], 0.25, xs[-1], 0.95, col = NSO_SEQ_COLS, border = NA); graphics::rect(0, 0.25, 1, 0.95, border = "grey40")
    graphics::rect(-0.09, 0.25, -0.03, 0.95, col = NSO_ZERO_COL, border = "grey40", xpd = NA); graphics::text(-0.06, 0.6, if (ko) "0" else "0", pos = 2, xpd = NA, cex = 1.1)
    graphics::text(c(0, 0.5, 1), 0.05, .nso_fmt(c(0, vmax / 2, vmax), ko, 3), xpd = NA, cex = 1.1)
    cb <- if (ko) sprintf("색 = 군집 평균 점수 (군집 안 모든 세포에 같은 색), 6 패널 같은 척도 0 – %s (최댓값)", .nso_fmt(vmax, ko, 3))
          else sprintf("colour = cluster-mean score (same colour for every cell of a cluster), one scale for all 6 panels 0 - %s (maximum)", .nso_fmt(vmax, ko, 3))
    graphics::mtext(cb, side = 3, line = 0.1, cex = 1.08); note("colorbar", "scale", cb)
    main <- if (ko) sprintf("atlas 6 단계 패널의 군집 평균 점수 (%s) — 단계별 최고 군집 표시, 같은 원 UMAP", partition)
            else sprintf("Cluster-mean scores of the 6 atlas stage panels (%s) - stage peak clusters marked, same original UMAP", partition)
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.3, cex = 2.0, font = 2); note("figure", "title", main)
    foot <- if (ko) c("검은 테 표식 = 그 단계 점수가 모든 군집 중 가장 높은 군집 (동률이면 모두). 표식 위치 = 그 군집 세포 중 UMAP 중앙값에 가장 가까운 세포 (표시 전용).",
                      "단계별 최고 군집은 군집 내부 1 위(그림 A)와 다른 비교다. 점수는 확률 · 확정 cell type 이 아니며, 패널 유전자는 다른 계열에서도 발현된다.")
            else c("Black-ringed marker = cluster with the highest score for this stage among all clusters (all if tied); placed at the cell nearest the cluster's UMAP median (display only).",
                   "A stage peak is a different comparison from the within-cluster top (figure A). Scores are not probabilities or confirmed cell types; panel genes are also expressed in other lineages.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.7 + (k - 1) * 1.4, cex = 1.08, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}
