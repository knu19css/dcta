# =============================================================================
# 50_mbv_mono_phase_view.R — monocyte 단계별 발현 구간: 그림
# =============================================================================
#
# 그림 1 fig_mono_phase_two_panel.png: 원 UMAP 같은 좌표 · 같은 창. 왼쪽 λ0 (저장 구조 기준 경로), 오른쪽 ordered 단계 경로, 같은 대표 key (기존 규칙).
#   세포 점 = 원저자 세부 라벨 사후 색, 이름표 = 군집 번호 + 최다 라벨, 화살표 = 우리 관측 경로 (탐색 순서), 화살표 이름표 = 배정 단계.
#   그 단계 근거가 있으면 단계 색, 결측이면 회색 "미해결" (동률 배치와 무관).
# 그림 2 fig_mono_phase_heatmap.png: 같은 key 의 두 경로. 행 = 선, 열 = 단계, 칸 색 = 학습 원점수 s, 굵은 테두리 = 배정, 아래 = 단계 손실 · 상태.
# 그림 3 fig_mono_phase_internal.png: ordered 경로 내부 B 마다 입구 / 출구 단계 유전자의 세 구간 평균 (학습 block, 세포 가중). 위상 상태와 평가 불가 사유.
#   기울기 · 구간 평균은 시간 속도 · flow · fate 가 아니다.
# 그린 글자는 모두 표 (figure_text) 로 돌려준다. 원저자 라벨은 계산 뒤에만 읽는다 (R49 단계 순서).
# =============================================================================

MBVP_PHASE_COLS <- c(E01 = "#1B9E77", E02 = "#D95F02", E03 = "#7570B3", M01 = "#E7298A", M02 = "#2C7FB8")
MBVP_MISSING_COL <- "grey55"
MBVP_FIG_TXT <- list(
  ko = list(left = "왼쪽: λ0 (저장 구조 기준 경로)", right = "오른쪽: ordered (단계 순서 발현 구간)", same = "두 경로가 같다", diff = "경로가 다르다",
            unres = "미해결", arrow = "화살표 = 고정 양끝 사이 우리 관측 경로의 탐색 순서 (분화 방향 · flow 검증 아님)", start = "출발 S", end = "끝 T",
            cells = "세포 색 = 원저자 세부 라벨 (사후, 경로 군집의 최다 라벨만)", other = "그 밖의 라벨", no_route = "경로 없음"),
  en = list(left = "Left: lambda0 (stored structural route)", right = "Right: ordered phase segments", same = "identical routes", diff = "routes differ",
            unres = "unresolved", arrow = "Arrows = search order of our observed route between fixed endpoints (not verified direction or flow)", start = "start S",
            end = "end T", cells = "Cell colours = author fine label (post hoc; route cluster top labels)", other = "other labels", no_route = "no route"))

.mbvp_note_factory <- function(fig) { txt <- list()
  list(note = function(panel, element, text) txt[[length(txt) + 1L]] <<- data.frame(figure = fig, panel = panel, element = element, text = text, stringsAsFactors = FALSE),
       table = function() do.call(rbind, txt)) }

#' 두 패널 그림. left / right = list(path, title, sub, phase (선마다 배정 단계), scored (선마다 그 단계 근거 여부)). 반환: 글자 표, 패널 범위
mbvp_draw_two_panel <- function(path, U, cluster, author, left, right, lab_of, top_of, main, foot, ko, fam) {
  cluster <- as.character(cluster); tx <- if (ko) MBVP_FIG_TXT$ko else MBVP_FIG_TXT$en; nf <- .mbvp_note_factory(basename(path)); note <- nf$note
  allc <- unique(c(left$path, right$path)); xy <- if (length(allc)) mbvu_node_xy(U, cluster, allc) else NULL
  P <- function(k) as.numeric(xy[match(k, xy$cluster), c("x_median", "y_median")])
  inz <- cluster %in% allc; pad <- function(r, f) r + c(-1, 1) * f * diff(r)
  b <- if (any(inz)) list(xlim = pad(stats::quantile(U[inz, 1], c(0.01, 0.99)), 0.14), ylim = pad(stats::quantile(U[inz, 2], c(0.01, 0.99)), 0.14)) else mbvu_bounds(U)
  pal <- mbvq_label_palette(top_of(allc)); grp <- ifelse(!is.na(author) & author %in% pal$author_label, author, "Other")
  win <- U[, 1] >= b$xlim[1] & U[, 1] <= b$xlim[2] & U[, 2] >= b$ylim[1] & U[, 2] <= b$ylim[2]; bounds <- list()
  grDevices::png(path, width = 2800, height = 1800, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(5.6, 1.45)); graphics::par(oma = c(0, 0, 3.6, 0), family = fam)
    for (side in c("left", "right")) { pn <- if (side == "left") left else right
      graphics::par(mar = c(0.4, 0.4, 4.6, 0.4)); bounds[[side]] <- mbvu_open_panel(b)
      o <- win & grp == "Other"; graphics::points(U[o, 1], U[o, 2], pch = 16, cex = 0.28, col = grDevices::adjustcolor("grey82", 0.6))
      n <- table(grp[win & grp != "Other"]); for (l in names(sort(n, decreasing = TRUE))) { i <- win & grp == l
        graphics::points(U[i, 1], U[i, 2], pch = 16, cex = 0.32, col = grDevices::adjustcolor(pal$color[pal$author_label == l], 0.75)) }
      p <- pn$path
      if (length(p) >= 2L) for (s in seq_len(length(p) - 1L)) .mbvq_arrow(P(p[s]), P(p[s + 1L]), "grey15", 5, b)
      if (length(p)) {
        bg <- pal$color[match(top_of(p), pal$author_label)]; bg[is.na(bg)] <- "white"
        q <- t(vapply(p, P, numeric(2))); graphics::points(q[, 1], q[, 2], pch = ifelse(seq_along(p) == 1L, 22, ifelse(seq_along(p) == length(p), 23, 21)), cex = 3.4, bg = bg,
                                                          col = "grey10", lwd = 2.5)
        if (length(p) >= 2L) for (s in seq_len(length(p) - 1L)) { m <- (P(p[s]) + P(p[s + 1L])) / 2; sc <- isTRUE(pn$scored[s]); ph <- pn$phase[s]
          t1 <- if (sc) ph else paste(ph, tx$unres); cl <- if (sc) MBVP_PHASE_COLS[[ph]] else MBVP_MISSING_COL
          w <- graphics::strwidth(t1, cex = 1.25, font = 2); h <- graphics::strheight(t1, cex = 1.25, font = 2)
          graphics::rect(m[1] - 0.6 * w, m[2] - 0.9 * h, m[1] + 0.6 * w, m[2] + 0.9 * h, col = if (sc) cl else "white", border = cl, lwd = 2.5, lty = if (sc) 1 else 2)
          graphics::text(m[1], m[2], t1, cex = 1.25, font = 2, col = if (sc) "white" else "grey30"); note(side, "edge_phase", paste(p[s], p[s + 1L], t1)) }
        lt <- paste0(p, "  ", vapply(p, function(k) { v <- lab_of(k); if (is.na(v)) "(no author label)" else v }, ""),
                     ifelse(seq_along(p) == 1L, paste0("  (", tx$start, ")"), ifelse(seq_along(p) == length(p), paste0("  (", tx$end, ")"), "")))
        bc <- pal$color[match(top_of(p), pal$author_label)]; bc[is.na(bc)] <- "grey40"
        lp <- .mbvt_draw_labels(q[, 1], q[, 2], lt, bc, 1.45, b); .mbvt_draw_boxes(lp, 1.45); for (t in lt) note(side, "node_label", t)
      } else { graphics::text(mean(b$xlim), mean(b$ylim), tx$no_route, cex = 2, col = "grey30"); note(side, "empty", tx$no_route) }
      graphics::box(col = "grey60"); graphics::title(main = pn$title, cex.main = 1.75, line = 2.6); graphics::mtext(pn$sub, side = 3, line = 0.6, cex = 1.15)
      note(side, "panel_title", pn$title); note(side, "panel_sub", pn$sub) }
    graphics::par(mar = c(0, 1, 0, 1)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    if (nrow(pal)) { graphics::legend(0.0, 1.0, legend = c(pal$author_label, tx$other), pch = 15, pt.cex = 2.4, col = c(pal$color, "grey82"), ncol = min(4L, nrow(pal) + 1L),
                                      cex = 1.2, bty = "n", title = tx$cells); for (t in c(pal$author_label, tx$other)) note("legend", "legend", t) }
    graphics::legend(0.70, 1.0, legend = c(names(MBVP_PHASE_COLS), tx$unres), fill = c(MBVP_PHASE_COLS, "white"), border = c(MBVP_PHASE_COLS, MBVP_MISSING_COL), cex = 1.2,
                     ncol = 3, bty = "n", title = if (ko) "화살표 이름표 = 배정 단계 (근거 있음 / 결측)" else "Arrow tag = assigned phase (scored / missing)")
    graphics::text(0.5, c(0.40, 0.23, 0.06), foot, cex = 1.25); for (t in foot) note("footer", "footer", t)
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.2, cex = 2.0, font = 2); note("figure", "title", main)
  }, finally = grDevices::dev.off())
  list(text = nf$table(), bounds = bounds)
}

#' 단계 heatmap. blocks = list(list(title, edge_label, smat (선 x 단계 원점수, 정순), assign (선마다 단계), status (단계마다), phase_loss (단계마다), L_phase))
mbvp_draw_phase_heatmap <- function(path, blocks, main, foot, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; ph <- MBVP_PHASES$phase
  ramp <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101); vcol <- function(v) ifelse(is.na(v), "grey88", ramp[round(50 * pmax(-1, pmin(1, v))) + 51])
  nr <- sum(vapply(blocks, function(k) length(k$edge_label) + 3L, 0L)); H <- max(1200L, 420L + 82L * nr)
  grDevices::png(path, width = 2800, height = H, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::par(mar = c(0.5, 0.5, 5, 0.5), family = fam); graphics::plot.new(); graphics::plot.window(c(0, 100), c(-4.2, nr + 0.8))
    x0 <- 40; cw <- 10.5; y <- nr
    for (k in blocks) {
      graphics::text(0.5, y, k$title, adj = c(0, 0.5), cex = 1.35, font = 2); note("table", "block_title", k$title)
      graphics::text(x0 + cw * (seq_along(ph) - 0.5), y, ph, cex = 1.3, font = 2, col = MBVP_PHASE_COLS[ph]); y <- y - 1
      for (r in seq_along(k$edge_label)) {
        graphics::text(0.5, y, k$edge_label[r], adj = c(0, 0.5), cex = 1.12); note("table", "edge_label", k$edge_label[r])
        for (j in seq_along(ph)) { v <- k$smat[r, j]; asg <- identical(k$assign[r], ph[j])
          graphics::rect(x0 + cw * (j - 1) + 0.3, y - 0.42, x0 + cw * j - 0.3, y + 0.42, col = vcol(v), border = if (asg) "black" else "white", lwd = if (asg) 4 else 1.5)
          t1 <- if (is.na(v)) (if (ko) "결측" else "miss.") else sprintf("%+.2f", v); t1 <- if (asg) paste0("[", t1, "]") else t1   # 맨 'NA' 는 표에서 결측으로 읽히므로 쓰지 않는다
          graphics::text(x0 + cw * (j - 0.5), y, t1, cex = 1.05, font = if (asg) 2 else 1); note("table", "cell", t1) }
        y <- y - 1 }
      lt <- vapply(seq_along(ph), function(j) sprintf("%s %.2f", if (k$status[j] == "scored") (if (ko) "근거" else "scored") else (if (ko) "미해결" else "unres."), k$phase_loss[j]), "")
      for (j in seq_along(ph)) { graphics::text(x0 + cw * (j - 0.5), y, lt[j], cex = 1.05, col = if (k$status[j] == "scored") MBVP_PHASE_COLS[[ph[j]]] else MBVP_MISSING_COL, font = 2)
        note("table", "phase_loss", lt[j]) }
      lab <- sprintf(if (ko) "단계 손실 (L_phase = %.3f)" else "phase loss (L_phase = %.3f)", k$L_phase); graphics::text(0.5, y, lab, adj = c(0, 0.5), cex = 1.12, font = 2); note("table", "L_phase", lab)
      y <- y - 2 }
    graphics::text(0.5, c(-1.3, -2.2, -3.1), foot, adj = c(0, 0.5), cex = 1.2); for (t in foot) note("footer", "footer", t)
    graphics::mtext(main, side = 3, line = 2, cex = 1.9, font = 2); note("figure", "title", main)
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 내부 B 그림. panels = list(list(title, sub, lines = data.frame(gene, phase, role (entry/exit), m1, m2, m3), note))
mbvp_draw_internal <- function(path, panels, main, foot, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; np <- max(1L, length(panels)); nc <- min(3L, np); nr <- ceiling(np / nc)
  grDevices::png(path, width = 2800, height = 300 + 760 * nr, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(rbind(matrix(seq_len(nr * nc), nr, nc, byrow = TRUE), nr * nc + 1L), heights = c(rep(5, nr), 1.2)); graphics::par(oma = c(0, 0, 3.2, 0), family = fam)
    for (i in seq_len(nr * nc)) { graphics::par(mar = c(3.8, 5.0, 6.2, 1))                    # 제목 · 부제 · 사유 줄이 겹치지 않게 위 여백을 넉넉히
      if (i > length(panels)) { graphics::plot.new(); next }
      pn <- panels[[i]]; L <- pn$lines; yl <- if (!is.null(L) && nrow(L) && any(is.finite(as.matrix(L[, c("m1", "m2", "m3")])))) range(as.matrix(L[, c("m1", "m2", "m3")]), na.rm = TRUE) else c(0, 1)
      yl <- yl + c(-0.08, 0.15) * max(diff(yl), 0.2)
      graphics::plot(NA, xlim = c(0.7, 3.3), ylim = yl, xaxt = "n", xlab = "", ylab = if (ko) "평균 log1p 발현 (학습 block)" else "mean log1p expression (training blocks)", cex.lab = 1.25, cex.axis = 1.15)
      graphics::axis(1, at = 1:3, labels = if (ko) c("A 쪽", "가운데", "C 쪽") else c("A side", "middle", "C side"), cex.axis = 1.25)
      if (!is.null(L) && nrow(L)) for (r in seq_len(nrow(L))) { cl <- MBVP_PHASE_COLS[[L$phase[r]]]
        graphics::lines(1:3, c(L$m1[r], L$m2[r], L$m3[r]), col = cl, lwd = 3, lty = if (L$role[r] == "entry") 1 else 2, type = "b", pch = 19)
        yy <- utils::tail(stats::na.omit(c(L$m1[r], L$m2[r], L$m3[r])), 1); if (length(yy)) graphics::text(3.05, yy, paste(L$gene[r], L$phase[r]), adj = c(0, 0.5), cex = 0.95, col = cl) }
      if (nzchar(pn$note)) graphics::mtext(pn$note, side = 3, line = 0.4, cex = 1.0, col = "grey25")
      graphics::title(main = pn$title, cex.main = 1.45, line = 3.9); graphics::mtext(pn$sub, side = 3, line = 2.0, cex = 1.1)
      for (t in c(pn$title, pn$sub, pn$note)) note(paste0("panel_", i), "panel_text", t) }
    graphics::par(mar = c(0, 1, 0, 1)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    graphics::text(0.5, c(0.75, 0.45, 0.15), foot, cex = 1.2); for (t in foot) note("footer", "footer", t)
    graphics::mtext(main, side = 3, outer = TRUE, line = 1, cex = 1.9, font = 2); note("figure", "title", main)
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 실행기에서 부르는 그림 단계: audit_09 세포표 (hash 확인; 여기서 처음 세포 저자 라벨을 읽는다), 기존 대표 key, 세 그림
.mbvp_figures <- function(x, a, out, rc, kr, ir, au, ck, log) {
  d9 <- file.path(x$attempt, a$stage); c9 <- .mbvq_manifest_file(d9, "audit_manifest.json", "cell_fine_stage_colors.tsv.gz"); inp <- c(cell_table = c9)
  ct <- utils::read.delim(gzfile(c9), stringsAsFactors = FALSE, colClasses = c(cell_id = "character", cluster = "character", author_label = "character", stage_group = "NULL", color = "NULL"))
  if (!identical(ct$cell_id, x$cell_id) || !identical(ct$cluster, x$cell_cluster)) stop("그림용 세포표의 세포 순서 · 군집이 NRV 세포와 다르다.")
  U <- cbind(ct$umap_x, ct$umap_y); ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; tx <- if (ko) MBVP_FIG_TXT$ko else MBVP_FIG_TXT$en
  fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  base <- rc[rc$condition == "lambda0" & rc$status == "route_ok", , drop = FALSE]
  rk <- if (nrow(base)) mbvq_mode_key(base$key, base$stored_path) else list(path = NA_character_, key = NA_character_)
  row_of <- function(cn) { z <- rc[rc$condition == cn & rc$key %in% rk$key & rc$status == "route_ok", , drop = FALSE]; if (nrow(z)) z[1, ] else NULL }
  r0 <- row_of("lambda0"); r1 <- row_of("ordered"); pa <- kr$assignments
  asg <- function(cn) { if (is.null(pa)) return(NULL); z <- pa[pa$key %in% rk$key & pa$condition == cn, , drop = FALSE]; z[order(z$step), , drop = FALSE] }
  a0 <- asg("lambda0"); a1 <- asg("ordered"); sp <- function(r) if (is.null(r)) character(0) else .mbvq_sp(r$path)
  rep <- data.frame(rule = "mode of stored monocyte primaries (ties by path string); first key name in radix order among keys using it; fixed before results",
                    key = rk$key, stored_mode_path = rk$path, n_keys_with_mode = if (nrow(base)) sum(base$stored_path == rk$path) else 0L, n_keys = length(unique(rc$key)),
                    lambda0_path = paste(sp(r0), collapse = ">"), ordered_path = paste(sp(r1), collapse = ">"), routes_identical = identical(sp(r0), sp(r1)), stringsAsFactors = FALSE)
  sub_of <- function(r, lab) if (is.null(r)) tx$no_route else sprintf(if (ko) "%s · %d 단계 · Σc/C %.3f · L_phase %.3f · J(%s) %.3f · 근거 단계 %d/5" else "%s · %d hops · Σc/C %.3f · L_phase %.3f · J(%s) %.3f · scored phases %d/5",
                                                                    r$key, r$hops, r$c_sum / r$c_scale, r$L_ordered, lab, r$J, r$n_phases_scored)
  left <- list(path = sp(r0), title = tx$left, sub = sub_of(r0, "λ0"), phase = if (is.null(a0)) character(0) else a0$phase, scored = if (is.null(a0)) logical(0) else a0$evidence == "scored")
  right <- list(path = sp(r1), title = tx$right, sub = sub_of(r1, "ordered"), phase = if (is.null(a1)) character(0) else a1$phase, scored = if (is.null(a1)) logical(0) else a1$evidence == "scored")
  oc <- rc[rc$condition == "ordered", ]
  main <- sprintf(if (ko) "monocyte 경로: 구조 기준 대 단계 순서 발현 구간 (%s, 대표 key %s, %s)" else "Monocyte route: structural baseline vs ordered phase segments (%s, key %s, %s)",
                  x$partition, rk$key, if (rep$routes_identical) tx$same else tx$diff)
  foot <- if (ko) c(paste0(tx$arrow, ". 같은 key · 같은 세포 좌표 · 같은 창. 사각 = 저장 S, 마름모 = 저장 T (두 패널 같은 양끝)."),
                    sprintf("ordered 경로가 λ0 과 다른 key: %d / %d (no_call %d). 대표 key 는 결과를 보기 전 규칙 (저장 경로 최빈값을 쓰는 key 중 이름순 첫째). 회색 = 배정은 됐지만 그 단계 근거 결측 (미해결).",
                            sum(oc$changed_vs_lambda0 %in% TRUE), nrow(oc), sum(oc$status != "route_ok")),
                    "이름표 = 군집 번호 + 원저자 최다 라벨 (사후). 단계 점수는 탐색적 순서 대응이며 분화 증거 · 확률 · 정체성 배정이 아니다. 왼쪽 단계 이름표는 같은 정순 규칙으로 기준 경로를 배정한 것.")
          else c(paste0(tx$arrow, ". Same key, coordinates and window. Square = stored S, diamond = stored T."),
                 sprintf("Keys whose ordered route differs from lambda0: %d / %d (no_call %d). Representative key fixed by rule before results. Grey = assigned but phase evidence missing.",
                         sum(oc$changed_vs_lambda0 %in% TRUE), nrow(oc), sum(oc$status != "route_ok")),
                 "Labels = cluster + author top label (post hoc). Phase scores are an exploratory order alignment, not differentiation evidence or identity assignment.")
  top_of <- function(k) unname(au$lab[k]); lab_of <- function(k) unname(au$short(k))
  p_two <- file.path(fig, "fig_mono_phase_two_panel.png")
  d2 <- mbvp_draw_two_panel(p_two, U, ct$cluster, ct$author_label, left, right, lab_of, top_of, main, foot, ko, fam)
  ck$add("figure_panels_same_window", identical(d2$bounds$left, d2$bounds$right), "par('usr') of the lambda0 and ordered panels are identical (same window, same coordinates)")
  # 그림 2: heatmap (두 경로)
  sm_of <- function(key) { es <- kr$edge_scores[kr$edge_scores$key %in% key, , drop = FALSE]; mbvp_score_matrix(es, unique(es$edge)) }
  blk <- function(r, aa, lab) { if (is.null(r) || is.null(aa) || !nrow(aa)) return(NULL); sm <- sm_of(r$key); e <- aa$edge
    st <- vapply(MBVP_PHASES$phase, function(p) r[[paste0("status_", p)]], ""); pl <- vapply(MBVP_PHASES$phase, function(p) as.numeric(r[[paste0("loss_", p)]]), 0)
    list(title = paste0(lab, ": ", r$path), edge_label = sprintf("%s  (%s → %s)", gsub(">", " → ", e, fixed = TRUE), vapply(aa$from, function(k) { v <- lab_of(k); if (is.na(v)) "-" else v }, ""),
                                                                  vapply(aa$to, function(k) { v <- lab_of(k); if (is.na(v)) "-" else v }, "")),
         smat = sm[e, MBVP_PHASES$phase, drop = FALSE], assign = aa$phase, status = st, phase_loss = pl, L_phase = r$L_condition) }
  blocks <- Filter(Negate(is.null), list(blk(r0, a0, if (ko) "λ0 (구조 기준) 경로, 정순 배정" else "lambda0 route, forward assignment"), blk(r1, a1, if (ko) "ordered 경로" else "ordered route")))
  foot2 <- if (ko) c("칸 색 = 학습 block 원점수 s = mean tanh(부호 x d_g) (파랑 = 문헌 방향과 반대, 빨강 = 같은 방향, 회색 NA = 근거 결측). 굵은 테두리 [ ] = 배정.",
                     "선 손실 = (1 - s)/2, 결측 = 1, 빈 단계 = 1. 아래 줄 = 단계 손실과 상태 (근거 = 배정 선에 그 단계 점수 있음, 미해결 = 없음 또는 빈 단계).",
                     "E03 = MPL/CD164 감소 (탐색적), M01 = IRF8 단일 유전자, M02 = T anchor 의존. 순서 E01→M02 는 실험 가설이다. 원저자 라벨은 사후 표시.")
           else c("Cell colour = training raw score s (blue = against the literature direction, red = along, grey NA = missing). Bold [ ] = assigned.",
                  "Edge loss = (1 - s)/2, missing = 1, empty phase = 1. Bottom row = phase loss and status (scored / unresolved).",
                  "E03 = MPL/CD164 down (exploratory), M01 = IRF8 single gene, M02 = T-anchor dependent. Phase order is a hypothesis. Author labels post hoc.")
  p_hm <- file.path(fig, "fig_mono_phase_heatmap.png")
  t2 <- mbvp_draw_phase_heatmap(p_hm, blocks, sprintf(if (ko) "선 x 단계 원점수와 연속 구간 배정 (%s, key %s)" else "Edge x phase scores and segment assignment (%s, key %s)", x$partition, rk$key), foot2, ko, fam)
  # 그림 3: 내부 B (ordered 경로)
  panels <- list(); is <- if (!is.null(ir$structure)) ir$structure else NULL
  if (!is.null(a1) && nrow(a1) >= 2L) for (pos in 2:nrow(a1)) { A <- a1$from[pos - 1L]; B <- a1$from[pos]; C <- a1$to[pos]; id <- paste(B, A, C, sep = "|")
    ent <- a1$phase[pos - 1L]; exi <- a1$phase[pos]; s1 <- if (!is.null(is)) is[is$triple == id, , drop = FALSE] else NULL
    lines <- NULL; nt <- ""
    if (!is.null(ir$stats)) { trb <- kr$keep[[rk$key]]$train
      for (role in c("entry", "exit")) { ph <- if (role == "entry") ent else exi; g <- MBVP_PHASE_GENES$gene[MBVP_PHASE_GENES$phase == ph]
        z <- ir$stats[ir$stats$triple == id & ir$stats$gene %in% g & ir$stats$block %in% trb, , drop = FALSE]
        for (gg in g) { zz <- z[z$gene == gg, , drop = FALSE]; wm <- function(m, n) { k <- is.finite(m) & n > 0; if (!any(k)) NA_real_ else sum(m[k] * n[k]) / sum(n[k]) }
          if (nrow(zz)) lines <- rbind(lines, data.frame(gene = gg, phase = ph, role = role, m1 = wm(zz$mean_A_side, zz$n_A_side), m2 = wm(zz$mean_middle, zz$n_middle),
                                                         m3 = wm(zz$mean_C_side, zz$n_C_side), stringsAsFactors = FALSE)) } }
      if (is.null(lines)) nt <- if (ko) "학습 block 의 내부 세포 없음" else "no interior cells in training blocks"
    } else nt <- paste0("not_evaluated: ", ir$raw_reason)
    stt <- if (!is.null(s1) && nrow(s1)) sprintf("%s (block span_ok %d/%d)", s1$pooled_structure[1], s1$n_blocks_span_ok[1], s1$n_blocks_with_B[1]) else ""
    panels[[length(panels) + 1L]] <- list(title = sprintf("B = %s  %s", B, { v <- lab_of(B); if (is.na(v)) "" else v }), sub = sprintf("%s → %s → %s · %s %s / %s %s · %s", A, B, C,
      if (ko) "입구" else "entry", ent, if (ko) "출구" else "exit", exi, stt), lines = lines, note = nt) }
  foot3 <- if (ko) c("ordered 경로의 내부 군집 B 마다: B 안 상대 위치 t (A 쪽 경계 0 → C 쪽 경계 1, 기존 NRV 세포 edge, pooled) 의 세 구간 평균 (학습 block 세포 가중).",
                     "실선 = 입구 선에 배정된 단계 유전자, 점선 = 출구 선 단계 유전자. 두 단계를 하나로 합치거나 단조를 강제하지 않는다. 감소 단계 (E03) 는 원래 값 그대로 그렸다.",
                     "기울기 · 구간 평균은 시간 속도 · flow · fate 가 아니다. 보류 block 값과 block 별 기울기는 internal_expression.tsv · internal_structure.tsv 에 있다.")
           else c("Per interior cluster B of the ordered route: three bins of relative position t inside B (pooled NRV cell edges; training blocks, cell weighted).",
                  "Solid = genes of the phase assigned to the entry edge, dashed = exit edge. Not merged, no monotonicity forced. Decreasing phase (E03) drawn as raw values.",
                  "Slopes and bin means are not velocity, flow or fate. Held-out and per-block values are in internal_expression.tsv / internal_structure.tsv.")
  p_in <- file.path(fig, "fig_mono_phase_internal.png")
  t3 <- mbvp_draw_internal(p_in, panels, sprintf(if (ko) "내부 B 의 단계 유전자 세 구간 평균 (%s, key %s, ordered 경로)" else "Interior B: three-bin means of phase genes (%s, key %s, ordered route)",
                                                 x$partition, rk$key), foot3, ko, fam)
  log("figures written | representative key:", rk$key, "| identical routes:", rep$routes_identical, "| internal panels:", length(panels))
  list(files = c(p_two, p_hm, p_in), text = rbind(d2$text, t2, t3), inputs = inp, representative = rep, ko = ko)
}
