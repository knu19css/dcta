# =============================================================================
# 61_nso_view.R — 호중구 성숙 단계 관측(NSO) 그림: 같은 원 UMAP 위 단계 marker 프로그램 6 개 · 관측 연결과 경계 · 내부 단절
# =============================================================================
#
# 계약: 계획 Task 3 (그림 2 개). 원 UMAP 은 표시 좌표일 뿐이다 — 추론 · edge · 경로를 만들지 않는다.
#   그림 1 fig01_neut_stage_evidence_umap.png: 6 패널 = 6 형태 단계 marker 프로그램(관측 유전자 logRNA 평균). 6 패널 공통 척도.
#     회색 = 그 패널 유전자가 하나도 검출되지 않은 세포. 입력에 유전자가 없는 패널은 값 대신 이유를 적는다. 어떤 패널도 단계 정답 · 예측이 아니다.
#   그림 2 fig02_observed_connections_umap.png: 실제 관측 군집 쌍(세포 edge ≥ 1)을 군집 중앙값 사이 직선으로 요약(굵기 ∝ q).
#     경계 세포(다른 군집과 edge ≥ 1) 강조. 붉은 테 = 군집 내부 단절(A–B–C 의 B 안 두 경계가 다른 연결 성분). 관측 쌍이 없으면 그렇게 적는다.
#   그림 글자는 figure_text 표로 남긴다. 결측은 "결측" / "miss." 로 쓰고 글자 "NA" 를 쓰지 않는다.
#   R34 mbvu_bounds · mbvu_open_panel, R36 한글 제목 글꼴, R50 글자 기록기를 수정 없이 재사용한다.
# =============================================================================

NSO_SEQ_COLS <- grDevices::hcl.colors(140, "YlOrRd", rev = TRUE)[41:140]                   # 가장 옅은 구간을 빼서 낮은 값도 회색 바탕과 구별되게 (attempt_01 그림 검토)
NSO_ZERO_COL <- "grey86"
NSO_BOUNDARY_COL <- "#4A90C2"
NSO_DISCONNECT_COL <- "#CB181D"

.nso_fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
.nso_fmt <- function(x, ko, d = 2) ifelse(is.na(x) | !is.finite(x), if (ko) "결측" else "miss.", formatC(x, format = "f", digits = d))
.nso_dots <- function(x) gsub("|", " · ", x, fixed = TRUE)

#' 경계 세포: 같은 분할에서 다른 군집 세포와 양수 edge 가 하나라도 있는 세포
nso_boundary_cells <- function(edges, cluster) {
  cluster <- as.character(cluster); x <- cluster[edges$i] != cluster[edges$j]
  out <- logical(length(cluster)); out[unique(c(edges$i[x], edges$j[x]))] <- TRUE; out
}

#' 그림 1. 6 형태 단계 marker 프로그램 (같은 원 UMAP, 공통 척도). 반환: 글자 표
nso_draw_stage_evidence <- function(path, U, pe, ref, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U); n <- nrow(U)
  ids <- paste0("marker_", NSO_STAGES); def <- pe$programs; cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  pos <- cm[!is.na(cm) & cm > 0]; vmax <- if (length(pos)) stats::quantile(pos, 0.99, names = FALSE) else 1
  colv <- function(v) NSO_SEQ_COLS[pmax(1L, pmin(100L, ceiling(v / vmax * 100)))]
  sm <- ref$shared_markers[ref$shared_markers$n_stages >= 2L, , drop = FALSE]
  shared <- if (nrow(sm)) paste(sprintf("%s (%s)", sm$gene, .nso_dots(sm$stages)), collapse = ", ") else if (ko) "없음" else "none"
  grDevices::png(path, width = 3000, height = 2400, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(c(1:6, 7, 7, 7), 3, byrow = TRUE), heights = c(5, 5, 1.1)); graphics::par(oma = c(4.2, 0, 3.8, 0), family = fam)
    for (k in seq_along(ids)) { id <- ids[k]; s <- NSO_STAGES[k]; d <- def[def$program_id == id, , drop = FALSE]; st <- ref$stages[ref$stages$stage_id == s, , drop = FALSE]
      graphics::par(mar = c(0.3, 0.3, 5.4, 0.3)); mbvu_open_panel(b)
      graphics::points(U[, 1], U[, 2], pch = 16, cex = 0.16, col = NSO_ZERO_COL)
      gtxt <- if (nzchar(d$genes)) .nso_dots(d$genes) else if (ko) "참고 유전자 없음" else "no reference genes"
      ttl <- if (startsWith(st$full_name, s)) st$full_name else sprintf("%s · %s", s, st$full_name)   # 단계 이름을 두 번 쓰지 않는다
      if (d$status != "ok") {
        msg <- if (ko) sprintf("표시할 값 없음: 입력에 유전자가 없다 (%s)", .nso_dots(d$missing_genes)) else sprintf("unavailable: genes not in input (%s)", .nso_dots(d$missing_genes))
        graphics::text(mean(b$xlim), mean(b$ylim), msg, cex = 1.35, font = 2, col = "grey20"); note(s, "unavailable", msg)
        sub <- if (ko) sprintf("marker %s · 관측 유전자 0 / %d", gtxt, d$n_expected) else sprintf("markers %s · observed genes 0 / %d", gtxt, d$n_expected)
      } else {
        v <- cm[, id]; det <- !is.na(cn[, id]) & cn[, id] > 0; o <- which(det)[order(v[det])]; nd <- sum(det)
        if (length(o)) graphics::points(U[o, 1], U[o, 2], pch = 16, cex = 0.34, col = colv(v[o]))            # 검출 세포는 바탕(0.16)보다 크게: 드문 단계 프로그램도 보이게
        sub <- if (ko) sprintf("marker %s · 관측 유전자 %d / %d · 1 개 이상 검출 세포 %s / %s (%s%%)", gtxt, d$n_observed, d$n_expected, .nso_fmt_n(nd), .nso_fmt_n(n), .nso_fmt(100 * nd / n, ko, 1))
               else sprintf("markers %s · observed genes %d / %d · cells with >= 1 detected %s / %s (%s%%)", gtxt, d$n_observed, d$n_expected, .nso_fmt_n(nd), .nso_fmt_n(n), .nso_fmt(100 * nd / n, ko, 1))
      }
      graphics::box(col = "grey60"); graphics::title(main = ttl, cex.main = 1.9, line = 3.2); graphics::mtext(sub, side = 3, line = 0.7, cex = 1.0)
      note(s, "title", ttl); note(s, "sub", sub) }
    graphics::par(mar = c(1.6, 22, 1.2, 22)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    xs <- seq(0, 1, length.out = 101); graphics::rect(xs[-101], 0.25, xs[-1], 0.95, col = NSO_SEQ_COLS, border = NA); graphics::rect(0, 0.25, 1, 0.95, border = "grey40")
    graphics::rect(-0.09, 0.25, -0.03, 0.95, col = NSO_ZERO_COL, border = "grey40", xpd = NA)
    graphics::text(-0.06, 0.6, if (ko) "미검출" else "none", pos = 2, xpd = NA, cex = 1.05)
    graphics::text(c(0, 0.5, 1), 0.05, .nso_fmt(c(0, vmax / 2, vmax), ko), xpd = NA, cex = 1.05)
    cb <- if (ko) sprintf("색 = 관측 marker 유전자 logRNA 평균, 6 패널 공통 척도 0 – %s (그 위 = 같은 끝색, 양수 값의 99 백분위수)", .nso_fmt(vmax, ko))
          else sprintf("colour = mean logRNA of observed marker genes, one scale for all 6 panels 0 - %s (above = same end colour; 99th percentile of positive values)", .nso_fmt(vmax, ko))
    graphics::mtext(cb, side = 3, line = 0.1, cex = 1.05); note("colorbar", "scale", cb)
    main <- if (ko) "호중구 형태 단계 marker 프로그램의 RNA 근거 — 같은 원 UMAP (단계 판정 아님)" else "RNA evidence of neutrophil morphological-stage marker programs - same original UMAP (not stage calls)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.3, cex = 2.0, font = 2); note("figure", "title", main)
    foot <- if (ko) c("logRNA = log1p(10^4 · count / 세포 전체 library). 원 RNA 고차원 값으로 계산했고 UMAP 은 표시 좌표일 뿐이다. 회색 = 그 패널 marker 가 하나도 검출되지 않은 세포.",
                      sprintf("여러 단계가 같이 쓰는 marker: %s → 패널끼리 겹칠 수 있다. 어떤 패널도 원저자 단계 정답이나 예측 단계가 아니다.", shared))
            else c("logRNA = log1p(1e4 * count / whole-cell library), computed from the original high-dimensional RNA; UMAP is display only. Grey = none of the panel's markers detected.",
                   sprintf("Markers shared by stages: %s, so panels can overlap. No panel is an author stage answer or a predicted stage.", shared))
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.6 + (k - 1) * 1.3, cex = 1.05, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}

#' 그림 2. 관측 연결 · 경계 세포 · 내부 단절 (같은 원 UMAP). pairs = 관측 군집 쌍 표, ic_clusters = nso_internal_continuity()$clusters
nso_draw_connections <- function(path, U, cluster, pairs, ic_clusters, boundary, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; b <- mbvu_bounds(U); cluster <- as.character(cluster)
  uc <- fdt_cluster_order_safe(cluster)
  xy <- t(vapply(uc, function(k) { i <- which(cluster == k); c(stats::median(U[i, 1]), stats::median(U[i, 2])) }, c(0, 0))); rownames(xy) <- uc
  disc <- ic_clusters$cluster[ic_clusters$n_internal_disconnected > 0]
  grDevices::png(path, width = 2700, height = 2800, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:2, 2), heights = c(6, 1.05)); graphics::par(oma = c(4.4, 0, 3.6, 0), family = fam)
    graphics::par(mar = c(0.3, 0.3, 0.3, 0.3)); mbvu_open_panel(b)
    graphics::points(U[!boundary, 1], U[!boundary, 2], pch = 16, cex = 0.16, col = NSO_ZERO_COL)
    graphics::points(U[boundary, 1], U[boundary, 2], pch = 16, cex = 0.16, col = grDevices::adjustcolor(NSO_BOUNDARY_COL, 0.55))
    note("map", "boundary_cells", sprintf("%s / %s cells have >= 1 edge to another cluster", .nso_fmt_n(sum(boundary)), .nso_fmt_n(length(boundary))))
    if (nrow(pairs)) { qmax <- max(pairs$q)
      for (r in order(pairs$q)) { a <- xy[pairs$cluster_a[r], ]; z <- xy[pairs$cluster_b[r], ]; f <- pairs$q[r] / qmax
        graphics::segments(a[1], a[2], z[1], z[2], lwd = 0.3 + 6 * f, col = grDevices::adjustcolor("grey10", 0.15 + 0.7 * f)) }
      note("map", "pairs", sprintf("%d observed cluster pairs drawn; q %s - %s", nrow(pairs), .nso_fmt(min(pairs$q), ko, 4), .nso_fmt(qmax, ko, 4)))
    } else {
      msg <- if (ko) "관측된 군집 쌍 없음 (세포 edge 가 모두 같은 군집 안에 있다)" else "no observed cluster pairs (all cell edges are within clusters)"
      graphics::text(mean(b$xlim), mean(b$ylim), msg, cex = 1.5, font = 2, col = "grey20"); note("map", "no_observed_pairs", msg)
    }
    graphics::points(xy[, 1], xy[, 2], pch = 21, bg = "white", col = "grey20", cex = 3.0, lwd = 1.2); graphics::text(xy[, 1], xy[, 2], uc, cex = 0.8, font = 2)
    if (length(disc)) { graphics::points(xy[disc, 1], xy[disc, 2], pch = 1, cex = 4.4, lwd = 3, col = NSO_DISCONNECT_COL)
      for (k in disc) { r <- ic_clusters[ic_clusters$cluster == k, , drop = FALSE]
        note("map", "internal_disconnected", sprintf("cluster %s: %d / %d neighbour pairs internally disconnected (%d components)", k, r$n_internal_disconnected, r$n_triples, r$n_components)) } }
    graphics::box(col = "grey60")
    lab <- if (ko) c("내부 세포 (다른 군집과 edge 없음)", "경계 세포 (다른 군집과 edge 1 개 이상)", "군집 쌍 요약선 (굵기 ∝ q, 세포가 이동한 길 아님)", "붉은 테 = 군집 내부 단절")
           else c("interior cell (no edge to another cluster)", "boundary cell (>= 1 edge to another cluster)", "cluster-pair summary line (width ~ q, not a cell path)", "red ring = internal disconnection")
    graphics::par(mar = c(0, 1, 0, 1)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
    graphics::legend(0.02, 0.95, legend = lab, pch = c(16, 16, NA, 1), lty = c(NA, NA, 1, NA), lwd = c(NA, NA, 4, 3), pt.cex = c(2, 2, 1, 2.6),
                     col = c(NSO_ZERO_COL, NSO_BOUNDARY_COL, "grey10", NSO_DISCONNECT_COL), ncol = 2, cex = 1.15, bty = "n")
    for (t in lab) note("legend", "legend", t)
    main <- if (ko) "관측 연결 · 경계 세포 · 내부 단절 (p92, 원 그래프 edge · 같은 원 UMAP)" else "Observed connections, boundary cells and internal breaks (p92, original graph edges, same UMAP)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.2, cex = 1.95, font = 2); note("figure", "title", main)
    foot <- if (ko) c(sprintf("선 = 실제 세포 edge 가 있는 군집 쌍 %d 개 전부를 군집 중앙값 사이 직선으로 요약한 것 (굵기 ∝ q = 상대 연결 비중). 세포가 이동한 길 · 분화 방향이 아니고 2 차원 거리로 만든 연결도 아니다.", nrow(pairs)),
                      "붉은 테 = 그 군집 안에서 서로 다른 두 이웃 군집 쪽 경계 세포가 다른 연결 성분에 있는 삼중쌍이 1 개 이상. 그래프 인접 = 내부 연속이 아니다.")
            else c(sprintf("Lines = all %d cluster pairs with real cell edges, summarised as straight lines between cluster medians (width ~ q, relative connection share); not cell paths, not directions, not 2-D distances.", nrow(pairs)),
                   "Red ring = at least one neighbour pair whose boundary cells lie in different connected components inside the cluster; graph adjacency is not internal continuity.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.6 + (k - 1) * 1.3, cex = 1.02, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}
