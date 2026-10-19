# =============================================================================
# 71_gse253355_view.R — GSE253355 그림 (PNG + PDF): 1) 표시 좌표 위 저자 주석 vs 군집별 단계 1 위 후보, 2) 과립구 라벨 확대 + 실제 군집 edge · 단계별 최고 군집 + 6 단계 점수 heatmap
# =============================================================================
#
# 계약: codex_gse253355_instruction.md 그림 항. 모두 계산 고정(calc_manifest.json) 뒤에 그린다. 저자 라벨은 사후 표시뿐이다(패널 · 기준 조정에 쓰지 않는다).
#   좌표가 새 DrL 레이아웃이면 "새 표시 좌표 (UMAP 아님)" 으로 적는다. 후보는 정답이 아니다. 결측 글자는 "결측" / "miss." 로 쓰고 "NA" 를 쓰지 않는다.
#   R34 mbvu_bounds · mbvu_open_panel, R36 한글 글꼴, R50 글자 기록기, R61 색, R67 단계 색 · 번호 배치를 수정 없이 재사용한다.
# =============================================================================

GSX_AUTHOR_COLS <- c(grDevices::hcl.colors(12, "Dark 3"), grDevices::hcl.colors(12, "Set 2"), grDevices::hcl.colors(12, "Dynamic"))

#' 장치 열기: .png = 비트맵, .pdf = quartz pdf (macOS, 한글 글꼴 이름 공유) 또는 cairo_pdf
.gsx_open_device <- function(path, w, h) {
  if (grepl("\\.pdf$", path)) { if (capabilities("aqua") || identical(getOption("bitmapType"), "quartz")) grDevices::quartz(file = path, type = "pdf", width = w / 100, height = h / 100, pointsize = 12)
                                else grDevices::cairo_pdf(path, width = w / 100, height = h / 100, pointsize = 12) }
  else grDevices::png(path, width = w, height = h, units = "px", res = 100, pointsize = 12)
}

#' 한 그림을 PNG 와 PDF 로 같은 함수로 그린다. draw(note) 는 글자를 note 로 남긴다 (PNG 때만 기록)
.gsx_draw_both <- function(base, w, h, draw) {
  nf <- .mbvp_note_factory(basename(paste0(base, ".png"))); out <- character(0)
  for (ext in c("png", "pdf")) { p <- paste0(base, ".", ext); .gsx_open_device(p, w, h)
    tryCatch(draw(if (ext == "png") nf$note else function(...) invisible(NULL)), finally = grDevices::dev.off()); out <- c(out, p) }
  list(files = out, text = nf$table())
}

.gsx_rep <- function(U, group, keep = NULL) { if (!is.null(keep)) { U <- U[keep, , drop = FALSE]; group <- group[keep] }
  r <- atlas_cl_rep_cells(U, group); r }

#' 그림 1. 왼쪽 = 저자 주석(사후), 오른쪽 = 군집별 단계 1 위 후보(색) + 군집 번호. coord_label = 좌표 이름(원 UMAP / 새 레이아웃)
gsx_draw_fig1 <- function(base, U, cluster, author, cand, coord_label, ko, fam) {
  b <- mbvu_bounds(U); al <- sort(unique(ifelse(is.na(author) | !nzchar(author), "(none)", author))); acol <- stats::setNames(rep(GSX_AUTHOR_COLS, length.out = length(al)), al)
  au <- ifelse(is.na(author) | !nzchar(author), "(none)", author); ccol <- .acl_top_cols(cand); names(ccol) <- cand$cluster
  .gsx_draw_both(base, 3400, 1950, function(note) {
    graphics::layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(10, 1.3)); graphics::par(oma = c(3.4, 0, 4.2, 0), family = fam)
    graphics::par(mar = c(0.4, 0.4, 4.6, 0.4)); mbvu_open_panel(b); graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.32, col = acol[au])
    rp <- .gsx_rep(U, au); cx <- 1.0; lp <- .acl_label_layout(rp, cx, b, rp$n_cells); .acl_draw_labels(lp, lp$text, cx)
    if (any(!lp$free)) note("author", "labels_overlapping", paste(lp$text[!lp$free], collapse = " | "))
    t1 <- if (ko) "저자 세포 주석 (사후 표시 · 정답 라벨 아님)" else "Author cell annotation (post hoc display; not ground truth)"
    graphics::box(col = "grey60"); graphics::title(main = t1, cex.main = 1.9, line = 2.6); graphics::mtext(coord_label, side = 3, line = 0.8, cex = 1.05, col = "grey20")
    note("author", "title", t1); note("author", "coords", coord_label)
    graphics::par(mar = c(0.4, 0.4, 4.6, 0.4)); mbvu_open_panel(b); graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.32, col = ccol[as.character(cluster)])
    rc <- .gsx_rep(U, as.character(cluster)); cx2 <- 1.25; lc <- .acl_label_layout(rc, cx2, b, rc$n_cells); .acl_draw_labels(lc, lc$text, cx2)
    if (any(!lc$free)) note("candidate", "labels_overlapping", paste(lc$text[!lc$free], collapse = " "))
    t2 <- if (ko) "군집별 atlas 단계 패널 1 위 후보 (군집 평균, 후보일 뿐 정답 아님)" else "Per-cluster top atlas stage panel (cluster mean; candidate, not truth)"
    tb <- table(factor(cand$score_top_candidate[cand$score_top_status == "score_top_candidate"], levels = NSO_STAGES))
    s2 <- paste0(if (ko) "1 위 후보 군집 수: " else "clusters by top: ", paste(sprintf("%s %d", NSO_STAGES, as.integer(tb)), collapse = " · "),
                 sprintf(if (ko) " · 동률 %d · 0 %d" else " · tied %d · zero %d", sum(cand$score_top_status == "tied"), sum(cand$score_top_status == "unassigned_zero")))
    graphics::box(col = "grey60"); graphics::title(main = t2, cex.main = 1.9, line = 2.6); graphics::mtext(s2, side = 3, line = 0.8, cex = 1.05)
    note("candidate", "title", t2); note("candidate", "sub", s2)
    graphics::par(mar = c(0.3, 2, 0.3, 2)); graphics::plot.new(); graphics::plot.window(c(0, 8), c(0, 1))
    labs <- c(NSO_STAGES, if (ko) "동률" else "tied", if (ko) "0 · 판정불가" else "zero / unavailable"); cols <- c(ATLAS_CL_STAGE_COL, ATLAS_CL_TIE_COL, ATLAS_CL_NONE_COL)
    xs <- seq(0.15, 7.15, length.out = 8); graphics::rect(xs, 0.3, xs + 0.2, 0.8, col = cols, border = "grey30"); graphics::text(xs + 0.26, 0.55, labs, adj = c(0, 0.5), cex = 1.3)
    for (k in seq_along(labs)) note("legend", "item", labs[k])
    main <- if (ko) "GSE253355 골수: 저자 주석과 atlas 단계 패널 군집 후보 (같은 표시 좌표)" else "GSE253355 marrow: author annotation and atlas stage-panel cluster candidates (same display coordinates)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.6, cex = 2.0, font = 2); note("figure", "title", main)
    foot <- if (ko) c("점수 = 군집 안 모든 세포의 패널 유전자 logRNA 평균 (패널 10 개 고정 분모). 후보는 확률 · 확정 cell type · 분화 판정이 아니다.",
                      "저자 주석은 계산을 고정한 뒤에 겹쳐 그린 사후 맥락이다. 패널 유전자는 호중구 밖 계열에서도 발현된다.")
            else c("Score = mean logRNA of panel genes over all cells of the cluster (fixed 10-gene denominator). Candidates are not probabilities, confirmed cell types or differentiation calls.",
                   "Author annotation is post hoc context drawn after calculations were fixed. Panel genes are also expressed outside the neutrophil lineage.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.4 + (k - 1) * 1.35, cex = 1.08, col = "grey20"); note("footer", "footer", foot[k]) }
  })
}

#' 그림 2. 왼쪽 = 과립구 라벨 세포 확대 (군집 1 위 후보 색) + 과립구 유도 그래프의 실제 군집 edge + 단계별 최고 군집, 오른쪽 = 군집 x 6 단계 점수 heatmap
gsx_draw_fig2 <- function(base, U, cluster, gran, S, cand, peaks, gedges, author_top, gran_labels, coord_label, ko, fam, min_frac = 0.5) {
  cluster <- as.character(cluster); ccol <- .acl_top_cols(cand); names(ccol) <- cand$cluster
  sub_n <- stats::setNames(as.integer(table(factor(cluster[gran], levels = cand$cluster))), cand$cluster); sub_f <- sub_n / cand$n_cells
  pkc <- unique(unlist(strsplit(peaks$peak_clusters[nzchar(peaks$peak_clusters)], "|", fixed = TRUE)))
  shown <- cand$cluster[(sub_f >= min_frac | cand$cluster %in% pkc) & sub_n > 0]                      # 표시 규칙: 과립구 라벨이 군집의 min_frac 이상이거나 단계별 최고 군집
  Ug <- U[gran, , drop = FALSE]; pad <- function(r) r + c(-1, 1) * 0.08 * diff(r)
  b <- if (nrow(Ug) >= 2L) list(xlim = pad(stats::quantile(Ug[, 1], c(0.005, 0.995))), ylim = pad(stats::quantile(Ug[, 2], c(0.005, 0.995)))) else mbvu_bounds(U)
  win <- U[, 1] >= b$xlim[1] & U[, 1] <= b$xlim[2] & U[, 2] >= b$ylim[1] & U[, 2] <= b$ylim[2]
  rg <- if (any(gran)) .gsx_rep(U, cluster, keep = gran) else NULL; vmax <- if (any(S > 0, na.rm = TRUE)) max(S, na.rm = TRUE) else 1
  hcol <- function(v) ifelse(is.na(v), "white", ifelse(v <= 0, NSO_ZERO_COL, NSO_SEQ_COLS[pmax(1L, pmin(100L, ceiling(sqrt(v / vmax) * 100)))]))   # 제곱근 색 척도
  .gsx_draw_both(base, 3400, 2100, function(note) {
    graphics::layout(matrix(c(1, 2), 1), widths = c(1.15, 1)); graphics::par(oma = c(4.2, 0, 4.2, 0), family = fam)
    graphics::par(mar = c(0.4, 0.4, 5.2, 0.4)); mbvu_open_panel(b)
    graphics::points(U[win & !gran, 1], U[win & !gran, 2], pch = 16, cex = 0.18, col = "grey88")
    if (any(gran)) graphics::points(U[gran, 1], U[gran, 2], pch = 16, cex = 0.45, col = ccol[cluster[gran]])
    if (!is.null(rg) && nrow(gedges)) { m <- match(gedges$cluster_a, rg$cluster); n2 <- match(gedges$cluster_b, rg$cluster); ok <- !is.na(m) & !is.na(n2) & gedges$cluster_a %in% shown & gedges$cluster_b %in% shown
      lw <- 0.8 + 2.2 * log10(gedges$n_cell_edges + 1) / max(1, log10(max(gedges$n_cell_edges) + 1))
      graphics::segments(rg$x[m[ok]], rg$y[m[ok]], rg$x[n2[ok]], rg$y[n2[ok]], lwd = lw[ok], col = grDevices::adjustcolor("grey15", 0.55)) }
    if (!is.null(rg)) { rs <- rg[rg$cluster %in% shown, , drop = FALSE]; rs$cluster <- sprintf("%s (%s)", rs$cluster, .nso_fmt_n(sub_n[rs$cluster])); cx <- 1.1
      lp <- .acl_label_layout(rs, cx, b, rs$n_cells); .acl_draw_labels(lp, lp$text, cx) }
    pk <- peaks[nzchar(peaks$peak_clusters), , drop = FALSE]; pl <- character(0)
    for (i in seq_len(nrow(pk))) for (k in strsplit(pk$peak_clusters[i], "|", fixed = TRUE)[[1]]) { j <- if (is.null(rg)) NA else match(k, rg$cluster)
      lab <- sprintf(if (ko) "%s 최고 = 군집 %s (과립구 라벨 %s / %s)" else "%s peak = cluster %s (granulocyte labels %s / %s)", pk$stage[i], k, .nso_fmt_n(sub_n[[k]]), .nso_fmt_n(cand$n_cells[cand$cluster == k]))
      if (!is.na(j)) { graphics::points(rg$x[j], rg$y[j], pch = 21, cex = 4.2, lwd = 3.6, col = ATLAS_CL_STAGE_COL[[pk$stage[i]]], bg = NA); pl <- c(pl, lab) }
      else pl <- c(pl, paste0(lab, if (ko) " (과립구 라벨 세포 없음 → 확대 밖)" else " (no granulocyte-label cells; outside zoom)")) }
    t1 <- if (ko) "과립구 · 호중구 관련 저자 라벨 세포 확대 (선택된 부분집합, 사후)" else "Zoom on granulocyte / neutrophil-associated author labels (selected subset, post hoc)"
    s1 <- paste0(if (ko) "포함 라벨: " else "labels: ", paste(gran_labels, collapse = " · "))
    s2 <- if (ko) sprintf("선 = 과립구 라벨 세포끼리의 실제 edge (굵기 ∝ log edge 수) · 번호 = 군집 (과립구 라벨 세포 수) · 표시 = 과립구 라벨 %d%% 이상 군집 + 단계별 최고 군집", round(100 * min_frac))
          else sprintf("lines = actual edges among granulocyte-label cells (width ~ log edges) · label = cluster (granulocyte-label cells) · shown = clusters >= %d%% granulocyte labels + stage peaks", round(100 * min_frac))
    graphics::box(col = "grey60"); graphics::title(main = t1, cex.main = 1.75, line = 3.6); graphics::mtext(s1, side = 3, line = 2.0, cex = 1.0); graphics::mtext(s2, side = 3, line = 0.6, cex = 0.95)
    note("zoom", "title", t1); note("zoom", "labels", s1); note("zoom", "sub", s2); for (k in seq_along(pl)) note("zoom", "peak", pl[k])
    if (length(pl)) graphics::legend("bottomleft", legend = pl, pch = 21, pt.cex = 2.4, pt.lwd = 2.5, col = ATLAS_CL_STAGE_COL[sub(" .*$", "", pl)], bty = "n", cex = 1.2, inset = 0.01)
    # heatmap
    ord <- order(match(cand$score_top_candidate, NSO_STAGES), -cand$score_top_value, method = "radix"); cl <- cand$cluster[ord]; nr <- length(cl)
    graphics::par(mar = c(7, 19, 5.2, 2)); graphics::plot(NA, xlim = c(0.5, 6.5), ylim = c(nr + 0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    for (r in seq_len(nr)) for (k in seq_along(NSO_STAGES)) { v <- S[cl[r], NSO_STAGES[k]]
      graphics::rect(k - 0.5, r - 0.5, k + 0.5, r + 0.5, col = hcol(v), border = "white") }
    tops <- strsplit(cand$score_top_candidate[ord], "|", fixed = TRUE)
    for (r in seq_len(nr)) { k <- match(tops[[r]], NSO_STAGES); k <- k[!is.na(k)]; if (length(k)) graphics::points(k, rep(r, length(k)), pch = 16, cex = 0.9, col = "black") }
    for (i in seq_len(nrow(pk))) for (k in strsplit(pk$peak_clusters[i], "|", fixed = TRUE)[[1]]) { r <- match(k, cl); if (!is.na(r)) graphics::rect(match(pk$stage[i], NSO_STAGES) - 0.5, r - 0.5, match(pk$stage[i], NSO_STAGES) + 0.5, r + 0.5, border = "black", lwd = 3) }
    rl <- sprintf("%s · %s (n %s)", cl, author_top[cl], .nso_fmt_n(cand$n_cells[ord])); graphics::axis(2, at = seq_len(nr), labels = rl, las = 1, cex.axis = if (nr > 40) 0.62 else 0.85, tick = FALSE)
    graphics::axis(3, at = seq_along(NSO_STAGES), labels = NSO_STAGES, tick = FALSE, cex.axis = 1.15)
    t2 <- if (ko) "군집 x 6 단계 패널 평균 점수 (모든 군집, 전체 자료 분모)" else "Cluster x 6 stage-panel mean scores (all clusters, whole-data denominators)"
    s3 <- sprintf(if (ko) "검은 점 = 군집 내부 1 위 · 굵은 테 = 단계별 최고 군집 · 색 = 제곱근 척도 0 – %s (6 단계 같은 척도) · 행 = 군집 · 저자 최다 라벨(사후)"
                  else "dot = within-cluster top · thick box = stage peak · colour = square-root scale 0 - %s (one scale) · row = cluster · author top label (post hoc)", .nso_fmt(vmax, ko, 3))
    xs <- seq(0.5, 6.5, length.out = 101); yb <- nr + 0.9; graphics::rect(xs[-101], yb, xs[-1], yb + 0.45, col = NSO_SEQ_COLS, border = NA, xpd = NA)
    tv <- c(0, vmax / 4, vmax); graphics::text(0.5 + 6 * sqrt(tv / vmax), yb + 0.95, .nso_fmt(tv, ko, 3), xpd = NA, cex = 1.05)
    graphics::title(main = t2, cex.main = 1.75, line = 3.6); graphics::mtext(s3, side = 3, line = 2.1, cex = 0.95); note("heatmap", "title", t2); note("heatmap", "sub", s3)
    main <- if (ko) "GSE253355: 과립구 라벨 확대 · 실제 군집 연결 · 6 단계 점수 (후보, 정답 아님)" else "GSE253355: granulocyte zoom, actual cluster connections, 6-stage scores (candidates, not truth)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.6, cex = 2.0, font = 2); note("figure", "title", main)
    foot <- if (ko) c(paste0(coord_label, ". 확대 · 표시 규칙은 저자 라벨로 고른 부분집합의 서술이며 패널 · 기준 · 점수 계산에 쓰지 않았다. 모든 군집 쌍 edge 는 cluster_edges.csv."),
                      "선이 있어도 분화 · 성숙 경로의 증거가 아니고, 선이 없다는 것만으로 생물학적 부재를 말하지 않는다. 다른 계열의 높은 점수 군집도 표에 그대로 남는다.")
            else c(paste0(coord_label, ". The zoom is a description of a subset chosen by author labels; labels were not used for panels, thresholds or scores."),
                   "A line is not evidence of differentiation; absence of a line alone is not biological absence. High-scoring clusters of other lineages stay in the tables.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.8 + (k - 1) * 1.35, cex = 1.05, col = "grey20"); note("footer", "footer", foot[k]) }
  })
}
