# =============================================================================
# 65_atlas_rna_view.R — atlas RNA 패널 그림: 같은 원 UMAP 위 atlas 6 단계 패널 · 기준선 NSO 패널과 크기 · 검출 · 순서 도달 비교
# =============================================================================
#
# 계약: plan Task 2 (panel 크기 · 연결성 · 원 UMAP 그림). 원 UMAP 은 표시 좌표일 뿐이다 — 추론 · edge · 경로를 만들지 않는다.
#   그림 1 fig01_atlas_panel_umap.png: 6 패널 = atlas 단계 양성 패널(관측 유전자 logRNA 평균), 6 패널 공통 척도. 회색 = 그 패널 유전자가 하나도 검출되지 않은 세포.
#     판정불가 패널(빈 패널 · target 에 유전자 없음)은 값 대신 이유를 적는다. 어떤 패널도 단계 정답 · 예측이 아니다.
#   그림 2 fig02_baseline_vs_atlas.png: 단계별 패널 크기(선형), ≥1 검출 세포 · ge1 순서 도달(전체 · 같은 block, log10(n + 1) 축). 판정불가는 막대 대신 글자.
#   그림 글자는 figure_text 표로 남긴다. 결측은 "판정불가" / "unavailable" 로 쓰고 글자 "NA" 를 쓰지 않는다.
#   R34 mbvu_bounds · mbvu_open_panel, R36 한글 제목 글꼴, R50 글자 기록기, R61 색 · 형식, R62 log 눈금을 수정 없이 재사용한다.
# =============================================================================

ATLAS_BASE_COL <- "grey62"
ATLAS_ATLAS_COL <- "#2F6DB5"
ATLAS_TARGET_COL <- "#173B66"

.atlas_unavail_txt <- function(ko) if (ko) "판정불가" else "unavailable"

#' 유전자 목록을 다섯 개씩 줄로 나눈다 (패널 위 글자가 칸을 넘지 않게)
.atlas_gene_lines <- function(genes, per = 5L) {
  g <- if (nzchar(genes)) strsplit(genes, "|", fixed = TRUE)[[1]] else character(0)
  if (!length(g)) return(character(0))
  vapply(split(g, ceiling(seq_along(g) / per)), paste, "", collapse = " · ")
}

#' 그림 1. atlas 6 단계 패널 (같은 원 UMAP, 공통 척도). det = detection_sets 표 (패널 상태). 반환: 글자 표
atlas_draw_panel_umap <- function(path, U, pe, nref, det, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U); n <- nrow(U)
  ids <- paste0("marker_", NSO_STAGES); def <- pe$programs; cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  pos <- cm[!is.na(cm) & cm > 0]; vmax <- if (length(pos)) stats::quantile(pos, 0.99, names = FALSE) else 1
  colv <- function(v) NSO_SEQ_COLS[pmax(1L, pmin(100L, ceiling(v / vmax * 100)))]
  sm <- nref$shared_markers
  shared <- if (nrow(sm)) paste(sprintf("%s (%s)", sm$gene, .nso_dots(sm$stages)), collapse = ", ") else if (ko) "없음" else "none"
  grDevices::png(path, width = 3000, height = 2600, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(c(1:6, 7, 7, 7), 3, byrow = TRUE), heights = c(5, 5, 1.1)); graphics::par(oma = c(4.6, 0, 3.8, 0), family = fam)
    for (k in seq_along(ids)) { id <- ids[k]; s <- NSO_STAGES[k]; d <- def[def$program_id == id, , drop = FALSE]; st <- nref$stages[nref$stages$stage_id == s, , drop = FALSE]
      dt <- det[det$stage == s, , drop = FALSE]
      graphics::par(mar = c(0.3, 0.3, 7.4, 0.3)); mbvu_open_panel(b)
      graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.16, col = NSO_ZERO_COL)
      ttl <- if (identical(st$full_name, s)) s else sprintf("%s · %s", s, st$full_name)
      gl <- .atlas_gene_lines(d$genes)
      if (d$status != "ok") {
        msg <- if (d$n_expected == 0L) { if (ko) "빈 패널: source 기준을 만족한 후보가 없다 (판정불가, 0 이 아니다)" else "empty panel: no source candidate met the rule (unavailable, not 0)" }
               else if (ko) sprintf("판정불가: target 에 패널 유전자가 없다 (%s)", .nso_dots(d$missing_genes)) else sprintf("unavailable: panel genes not in target (%s)", .nso_dots(d$missing_genes))
        graphics::text(mean(b$xlim), mean(b$ylim), msg, cex = 1.3, font = 2, col = "grey20"); note(s, "unavailable", msg)
        sub <- if (ko) sprintf("관측 유전자 0 / %d (source 패널 상태 %s)", d$n_expected, dt$panel_status) else sprintf("observed genes 0 / %d (source panel status %s)", d$n_expected, dt$panel_status)
      } else {
        v <- cm[, id]; dd <- !is.na(cn[, id]) & cn[, id] > 0; o <- which(dd)[order(v[dd])]; nd <- sum(dd)
        if (length(o)) graphics::points(U[o, 1], U[o, 2], pch = 16, cex = 0.34, col = colv(v[o]))
        miss <- if (nzchar(d$missing_genes)) sprintf(if (ko) " · target 에 없음 %s" else " · not in target %s", .nso_dots(d$missing_genes)) else ""
        sub <- if (ko) sprintf("관측 유전자 %d / %d%s · 1 개 이상 검출 세포 %s / %s (%s%%)", d$n_observed, d$n_expected, miss, .nso_fmt_n(nd), .nso_fmt_n(n), .nso_fmt(100 * nd / n, ko, 1))
               else sprintf("observed genes %d / %d%s · cells with >= 1 detected %s / %s (%s%%)", d$n_observed, d$n_expected, miss, .nso_fmt_n(nd), .nso_fmt_n(n), .nso_fmt(100 * nd / n, ko, 1))
      }
      graphics::box(col = "grey60"); graphics::title(main = ttl, cex.main = 1.9, line = 5.3); note(s, "title", ttl)
      for (j in seq_along(gl)) { graphics::mtext(gl[j], side = 3, line = 3.4 - (j - 1) * 1.25, cex = 0.98, col = "grey15"); note(s, paste0("genes_", j), gl[j]) }
      if (!length(gl)) { gtx <- if (ko) "(패널 유전자 없음)" else "(no panel genes)"; graphics::mtext(gtx, side = 3, line = 3.4, cex = 0.98, col = "grey15"); note(s, "genes_1", gtx) }
      graphics::mtext(sub, side = 3, line = 0.5, cex = 1.0); note(s, "sub", sub) }
    graphics::par(mar = c(1.6, 22, 1.2, 22)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    xs <- seq(0, 1, length.out = 101); graphics::rect(xs[-101], 0.25, xs[-1], 0.95, col = NSO_SEQ_COLS, border = NA); graphics::rect(0, 0.25, 1, 0.95, border = "grey40")
    graphics::rect(-0.09, 0.25, -0.03, 0.95, col = NSO_ZERO_COL, border = "grey40", xpd = NA)
    graphics::text(-0.06, 0.6, if (ko) "미검출" else "none", pos = 2, xpd = NA, cex = 1.05)
    graphics::text(c(0, 0.5, 1), 0.05, .nso_fmt(c(0, vmax / 2, vmax), ko), xpd = NA, cex = 1.05)
    cb <- if (ko) sprintf("색 = 관측 패널 유전자 logRNA 평균, 6 패널 공통 척도 0 – %s (그 위 = 같은 끝색, 양수 값의 99 백분위수)", .nso_fmt(vmax, ko))
          else sprintf("colour = mean logRNA of observed panel genes, one scale for all 6 panels 0 - %s (above = same end colour; 99th percentile of positive values)", .nso_fmt(vmax, ko))
    graphics::mtext(cb, side = 3, line = 0.1, cex = 1.05); note("colorbar", "scale", cb)
    main <- if (ko) "atlas RNA 단계 양성 패널 (GSE294330, source 전용 고정) — 같은 원 UMAP 위 RNA 근거 (단계 판정 아님)"
            else "Atlas RNA stage-positive panels (GSE294330, fixed from source only) - RNA evidence on the same original UMAP (not stage calls)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.3, cex = 2.0, font = 2); note("figure", "title", main)
    foot <- if (ko) c("패널 = bulk RNA atlas 에서 인접 단계보다 높은 유전자 (nominal code 중앙값 대비, DCTA 실험 규칙 · 저자 공식 패널 아님). logRNA = log1p(10^4 · count / 세포 전체 library), UMAP 은 표시 좌표일 뿐이다.",
                      sprintf("회색 = 그 패널 유전자가 하나도 검출되지 않은 세포. 여러 단계 공유 유전자: %s. 어떤 패널도 원저자 단계 정답이나 예측 단계가 아니다.", shared))
            else c("Panel = genes higher than adjacent stages in the bulk RNA atlas (nominal-code medians; DCTA experimental rule, not the authors' panel). logRNA = log1p(1e4 * count / whole-cell library); UMAP is display only.",
                   sprintf("Grey = none of the panel genes detected. Genes shared by stages: %s. No panel is an author stage answer or a predicted stage.", shared))
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.7 + (k - 1) * 1.4, cex = 1.05, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 묶은 막대 하나 (기준선 · atlas). y = 값 (log 축이면 log10(값 + 1)), 판정불가(atlas 결측)는 막대 대신 글자
.atlas_pair_bars <- function(x, yb, ya, lb, la, ymax, ko, note, panel, inner = NULL) {
  w <- 0.36
  graphics::rect(x - 0.4, 0, x - 0.4 + w, yb, col = ATLAS_BASE_COL, border = "grey40")
  graphics::text(x - 0.4 + w / 2, yb, lb, pos = 3, cex = 1.0, col = "grey25")
  ok <- !is.na(ya)
  if (any(ok)) { graphics::rect(x[ok] + 0.04, 0, x[ok] + 0.04 + w, ya[ok], col = ATLAS_ATLAS_COL, border = "grey20")
    if (!is.null(inner)) graphics::rect(x[ok] + 0.04 + 0.08, 0, x[ok] + 0.04 + w - 0.08, inner[ok], col = ATLAS_TARGET_COL, border = NA)
    graphics::text(x[ok] + 0.04 + w / 2, ya[ok], la[ok], pos = 3, cex = 1.0, col = ATLAS_TARGET_COL, font = 2) }
  if (any(!ok)) { u <- .atlas_unavail_txt(ko); graphics::text(x[!ok] + 0.04 + w / 2, ymax * 0.06, u, srt = 90, adj = c(0, 0.5), cex = 1.05, col = "#CB181D", font = 2)
    for (k in which(!ok)) note(panel, paste0("unavailable_", k), u) }
  for (k in seq_along(x)) { note(panel, paste0("baseline_", k), lb[k]); if (ok[k]) note(panel, paste0("atlas_", k), la[k]) }
}

#' 그림 2. 기준선 NSO 패널 vs atlas 패널: 패널 크기 · ≥1 검출 세포 · ge1 순서 도달(전체 · 같은 block). cmp = baseline_comparison, rcmp = reach_comparison
atlas_draw_comparison <- function(path, cmp, rcmp, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; st <- cmp$stage; x <- seq_along(st); fn <- function(v) .nca_fmt_n(v)
  lab <- function(v) ifelse(is.na(v), .atlas_unavail_txt(ko), fn(v))
  grDevices::png(path, width = 2800, height = 2100, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::par(mfrow = c(2, 2), oma = c(5.2, 0, 4.4, 0), mar = c(4.2, 9, 5, 1.5), family = fam)   # 왼쪽 여백: "1,000,000" 눈금과 축 제목이 겹치지 않게 (attempt_01 그림 검토)
    # A. 패널 크기 (선형)
    ymax <- max(c(cmp$baseline_panel_size, cmp$atlas_panel_size, 1)) * 1.25
    graphics::plot(NA, xlim = c(0.4, length(st) + 0.6), ylim = c(0, ymax), xaxt = "n", xlab = "", ylab = "", las = 1)
    graphics::title(ylab = if (ko) "패널 유전자 수" else "panel genes", line = 3.2, cex.lab = 1.15)
    graphics::axis(1, at = x, labels = st, cex.axis = 1.25)
    .atlas_pair_bars(x, cmp$baseline_panel_size, ifelse(cmp$atlas_panel_size > 0, cmp$atlas_panel_size, NA), fn(cmp$baseline_panel_size),
                     sprintf("%s (%s)", fn(cmp$atlas_panel_size), fn(cmp$atlas_panel_in_target)), ymax, ko, note, "panel_size", inner = cmp$atlas_panel_in_target)
    ta <- if (ko) "A. 단계 패널 크기" else "A. Stage panel size"; sa <- if (ko) "atlas 막대 안 짙은 칸 · 괄호 = target 에 있는 유전자 수" else "dark inner bar and brackets = genes present in target"
    graphics::title(main = ta, cex.main = 1.6, line = 2.4); graphics::mtext(sa, side = 3, line = 0.6, cex = 1.05); note("panel_size", "title", ta); note("panel_size", "sub", sa)
    # B–D. log10(n + 1) 축
    lg <- function(v) log10(v + 1)
    one <- function(yb, ya, ttl, sub, panel, fb = NULL) {
      ymax <- lg(max(c(yb, ya, 1), na.rm = TRUE)) * 1.2
      graphics::plot(NA, xlim = c(0.4, length(st) + 0.6), ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      graphics::title(ylab = if (ko) "세포 수 (축 = log10(n + 1), 눈금 = 원 세포 수)" else "cells (axis = log10(n + 1), ticks = raw counts)", line = 6.6, cex.lab = 1.15)
      graphics::axis(1, at = x, labels = st, cex.axis = 1.25); tk <- .nca_log_ticks(ymax); graphics::axis(2, at = tk$at, labels = tk$label, las = 1)
      .atlas_pair_bars(x, lg(yb), ifelse(is.na(ya), NA, lg(ya)), fn(yb), lab(ya), ymax, ko, note, panel)
      graphics::title(main = ttl, cex.main = 1.6, line = 2.4); graphics::mtext(sub, side = 3, line = 0.6, cex = 1.05, col = "grey15"); note(panel, "title", ttl); note(panel, "sub", sub) }
    one(cmp$baseline_n_ge1, cmp$atlas_n_ge1, if (ko) "B. 패널 유전자 1 개 이상 검출 세포" else "B. Cells with >= 1 panel gene detected",
        if (ko) "같은 전체 세포 분모 · 패널이 크면 검출 세포가 늘 수 있다" else "same all-cell denominator; a larger panel can detect more cells", "detected_ge1")
    for (gr in c("all_edges", "same_block")) { z <- rcmp[rcmp$criterion == "ge1" & rcmp$graph == gr, , drop = FALSE]; z <- z[match(st, z$stage), , drop = FALSE]
      fb_b <- z$baseline_first_break[1]; fb_a <- z$atlas_first_break[1]; kd <- z$atlas_break_kind[1]
      fbt <- function(v) if (nzchar(v)) v else if (ko) "끊김 없음" else "no break"
      kdt <- if (!nzchar(kd)) "" else if (kd == "unavailable_stage") (if (ko) " (판정불가 단계)" else " (unavailable stage)") else (if (ko) " (실제 edge 없음)" else " (no real edge)")
      sub <- sprintf(if (ko) "첫 끊김: 기준선 %s · atlas %s%s" else "first break: baseline %s · atlas %s%s", fbt(fb_b), fbt(fb_a), kdt)
      ttl <- if (gr == "all_edges") (if (ko) "C. ge1 순서 도달 · 모든 원 edge" else "C. ge1 ordered reach · all original edges") else (if (ko) "D. ge1 순서 도달 · 같은 block edge 만" else "D. ge1 ordered reach · same-block edges only")
      one(z$baseline_n_reachable_cells, z$atlas_n_reachable_cells, ttl, sub, paste0("reach_", gr)) }
    graphics::par(xpd = NA)
    main <- if (ko) "기준선 NSO 패널 vs atlas RNA 패널 — 패널 크기가 달라 검출 · 연결 증가는 단계 정확도가 아니다" else "Baseline NSO panels vs atlas RNA panels - panel sizes differ; more detection or connection is not stage accuracy"
    graphics::mtext(main, side = 3, outer = TRUE, line = 2.2, cex = 1.85, font = 2); note("figure", "title", main)
    leg <- if (ko) c("회색 = 기준선 NSO 패널 (attempt_02 · connectivity_audit_02)", "파랑 = atlas RNA 패널 (이번 실행)") else c("grey = baseline NSO panels (attempt_02, connectivity_audit_02)", "blue = atlas RNA panels (this run)")
    graphics::mtext(paste(leg, collapse = "     "), side = 3, outer = TRUE, line = 0.4, cex = 1.2); for (k in seq_along(leg)) note("figure", "legend", leg[k])
    foot <- if (ko) c("순서 도달 = R62 규칙 그대로: 첫 단계 모든 검출 세포에서 연속 단계 사이 실제 원 edge 로만 전진 (단계 유도 성분 단위, zero-hop 전진 없음). 판정불가 단계 뒤의 도달은 평가하지 않는다.",
                      "판정불가 = 빈 패널 · target 유전자 없음 · 기준보다 적은 관측 유전자 (0 이 아니다). 어떤 막대도 단계 정확도 · 확률 · 분화 방향이 아니다.")
            else c("Ordered reach = R62 rule: from all first-stage cells, advance only through real original edges between consecutive stages (stage-induced components, no zero-hop advance). Reach after an unavailable stage is not evaluated.",
                   "Unavailable = empty panel, genes missing in target or fewer observed genes than the criterion (not 0). No bar is stage accuracy, probability or differentiation direction.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 1.2 + (k - 1) * 1.4, cex = 1.05, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}
