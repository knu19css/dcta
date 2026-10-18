# =============================================================================
# render_rna_dag_review.R — RNA 중심 DAG 검토안 그림 (원본 대 제안 · 핵심안 · 계통별 상세), PNG 와 PDF
# =============================================================================
# rna_dag_original_vs_proposed : 원본 76 노드 · 89 연결이 제안 항목 · 처리로 어떻게 대응하는지, HSPC 분리 개요
# rna_dag_core                 : 제안 view 전체. RNA 상태 후보 · 산물 · 보류는 계통 묶음에, assay 정의 · taxonomy 는 아래 참고 panel
# rna_dag_lineage_<group>      : 계통별 상세. 노드 ID · 원본 이름 · 배지, 왼쪽 참고 열(assay 정의), 오른쪽 taxonomy 열(하위 목록), 연결 표
# 읽는 법: 세로 층은 원본 연결의 위상 순서(앞선 상태가 위)일 뿐 분화 시간 · 속도가 아니다. 가로 위치는 겹치지 않게 놓은 배치다.
# 배치(spec)와 그리기를 나눈다. rdr_core_spec / rdr_lineage_spec 은 장치 없이 좌표 · 경로만 계산하므로 테스트가 직접 검사한다:
#   - 출발점 미해결 연결은 HSC · MPP 묶음 상자 아래 변 한 점에서 나와 어느 후보 상자도 지나지 않는다 (한쪽 후보에 붙은 것처럼 보이지 않게).
#   - 그리는 연결은 자기 끝 상자 말고 다른 상자를 지나지 않는다 (직선 -> 층 사이 틈으로 꺾은 선 -> 옆으로 휜 곡선 순서로 고른다).
#   - 범례 · 상자는 캔버스 안에 있다.
# 사용: Rscript R/rna_dag_review/render_rna_dag_review.R [--source DIR] [--proposal DIR] [--out DIR] [--lang ko|en]
#   --proposal 의 표가 원본에서 새로 만든 표와 다르면 멈춘다 (빌더를 먼저 다시 실행).
# =============================================================================

RDR_LINEAGE_GROUPS <- list(
  mk_ery = list(ko = "거핵구 · 적혈구", en = "Megakaryocyte · erythroid", modules = c("MK", "Ery")),
  granulocyte = list(ko = "호중구 · 호산구 · 호염기구 / 비만세포", en = "Neutrophil · eosinophil · basophil / mast", modules = c("Neut", "Eos", "BasoMast")),
  mono_dc = list(ko = "단핵구 / 대식세포 · 수지상세포", en = "Monocyte / macrophage · dendritic cells", modules = c("MonoMac", "DC")),
  b = list(ko = "B 세포", en = "B cells", modules = "B"),
  t = list(ko = "T 세포 (흉선 · 말초)", en = "T cells (thymus · periphery)", modules = c("T_thymic", "T_periph")),
  nk_ilc = list(ko = "NK · ILC", en = "NK · ILC", modules = "NK_ILC"))

RDR_NODE_STYLE <- list(
  rna_state_candidate = list(fill = "#FFFFFF", border = "#333333", lty = 1, lwd = 1.1),
  split = list(fill = "#FFF4D6", border = "#B07D00", lty = 1, lwd = 1.7),
  assay_defined_reference = list(fill = "#EFEFEF", border = "#8A8A8A", lty = 2, lwd = 1),
  product_reference = list(fill = "#FFFFFF", border = "#7A7A7A", lty = 3, lwd = 1.2),
  held_candidate = list(fill = "#FBE9E7", border = "#C0392B", lty = 3, lwd = 1.3),
  taxonomy_reference = list(fill = "#E6F4F4", border = "#2A8C8C", lty = 4, lwd = 1))
RDR_EDGE_STYLE <- list(
  keep = list(col = "#333333", lty = 1, lwd = 1.2, arrow = TRUE),
  unresolved = list(col = "#3B6FB6", lty = 2, lwd = 1.1, arrow = TRUE),
  reference = list(col = "#9C9C9C", lty = 1, lwd = 0.9, arrow = TRUE),
  product = list(col = "#7A7A7A", lty = 3, lwd = 1.2, arrow = TRUE),
  held = list(col = "#C0392B", lty = 3, lwd = 1.3, arrow = TRUE),
  membership = list(col = "#8A8A8A", lty = 3, lwd = 0.9, arrow = FALSE))
RDR_BADGE_COL <- "#E67E22"
RDR_CORE <- list(slot = 1.12, dy = 0.58, bw = 0.98, bh = 0.34, gap_g = 0.45, top = 2.3)
RDR_LIN <- list(slot = 1.55, dy = 0.9, bw = 1.36, bh = 0.56, ref_w = 2.25, tax_w = 2.1, tab_w = 7.6)
RDR_LEGEND_COLW <- c(1.9, 2.35, 3.05, 2.85)
RDR_LEGEND_ROWH <- 0.2

RDR_TXT <- list(
  ko = list(
    lang = "ko",
    treat = c(keep = "유지 제안", unresolved = "출발점 미해결", reference_only = "참고", held = "보류", taxonomy = "분류(is_a)", membership = "gate 포함"),
    core_title = "RNA 중심 조혈 DAG 개정안 — 핵심안 (사용자 검토 대기)",
    core_sub = "원본 v1.1.0 의 76 노드 · 89 연결을 모두 대응. 제안 77 항목: RNA 상태 후보 %d · assay 정의 %d · 산물 %d · 보류 %d · taxonomy %d. 모든 user_approval = pending",
    origin = "HSC · MPP 별도 후보 (원본 HSPC_continuum 을 나눔)",
    origin_note = "파란 점선 = 원본 HSPC 출발 연결. 묶음 상자 아래에서 나오며 어느 후보에서 출발하는지 정하지 않았고 두 후보에 복제하지 않음",
    layer_note = "세로 층 = 원본 연결의 위상 순서 (분화 시간 · 속도가 아님). 가로 위치 = 겹치지 않게 놓은 배치",
    lg_keep = "유지 제안 연결 (원본 방향 근거, 승인 대기)", lg_unres = "출발점 미해결 (HSC/MPP 중 미정)", lg_held = "원본 보류 (held)",
    lg_ref = "참고 연결 (assay 정의 끝점)", lg_product = "산물 생성 (produces)", lg_member = "gate 포함 관계 (발달 아님)",
    lg_rna = "RNA 상태 후보", lg_split = "HSC · MPP 후보", lg_assay = "assay 정의 (참고)", lg_prodnode = "산물", lg_heldnode = "보류 후보",
    lg_taxnode = "taxonomy (참고)", lg_badge = "단백질 · 기능 정의 의존 (배지)",
    ref_panel = "참고 panel — assay 정의 7 · taxonomy 4 (RNA 상태 후보와 분리, 삭제 아님)",
    ref_unres = "출발점 미해결", tax_children = "하위", membership = "포함 (발달 관계 아님)",
    lin_title = "계통별 상세 — %s (사용자 검토 대기)",
    lin_sub = "왼쪽 회색 열 = assay 정의 참고 항목, 오른쪽 청록 열 = taxonomy (하위 목록). 노드 둘째 줄 = 원본 이름. 주황 표시 = 단백질 · 기능 정의 의존",
    ref_col = "참고: assay 정의", tax_col = "참고: taxonomy",
    tab_head = c("연결", "원본 연결", "처리", "원본 support / adult", "원본 근거 ID"),
    tab_title = "이 그림의 원본 연결 %d 개 (제안 처리 · 원본 근거 등급은 원본 그대로)",
    cmp_title = "원본 v1.1.0 대 제안 — 무엇이 바뀌는가 (사용자 검토 대기)",
    cmp_nodes = "노드 76 → 제안 항목 77", cmp_edges = "연결 89 → 처리",
    cmp_left_nodes = "원본 node_kind", cmp_right_nodes = "제안 역할", cmp_left_edges = "원본 support_class", cmp_right_edges = "제안 처리",
    role = c(split = "HSC · MPP 후보 (분리)", rna_state_candidate = "RNA 상태 후보 (유지)", assay_defined_reference = "assay 정의 (참고)",
             product_reference = "산물", held_candidate = "보류", taxonomy_reference = "taxonomy (참고)"),
    split_line = "분리 결정 · 연결 미승인",
    hspc_left = "원본: HSPC_continuum 하나에서 coarse 관계 19 개", hspc_right = "제안: HSC_candidate · MPP_candidate 를 별도 후보로 두고 19 개 관계의 출발점은 미해결로 보존",
    hspc_note = "복제하지 않음 · HSC→MPP 연결은 원본에 없어 그리지 않음 · 마커 · 연결 승인은 별도",
    footer = "RNA 구별 상태는 평가하지 않음 (not_assessed; HSC · MPP 는 unresolved). 원본 layer · adult_support · 근거는 덮어쓰지 않음. 표: data/priors/rna_dag_review/proposal/"),
  en = list(
    lang = "en",
    treat = c(keep = "keep (pending)", unresolved = "origin unresolved", reference_only = "reference", held = "held", taxonomy = "is_a", membership = "gate membership"),
    core_title = "RNA-centred haematopoietic DAG proposal — core view (pending user review)",
    core_sub = "All 76 nodes and 89 edges of v1.1.0 are traced. 77 proposed items: RNA state candidates %d · assay-defined %d · product %d · held %d · taxonomy %d. user_approval = pending for all",
    origin = "HSC and MPP as separate candidates (split from HSPC_continuum)",
    origin_note = "Blue dashed = source HSPC-start relations. They leave the pair box; the origin is unresolved and not copied to both candidates",
    layer_note = "Vertical tiers = topological order of source edges (not time or speed). Horizontal position = layout only",
    lg_keep = "kept relation (source direction evidence, pending)", lg_unres = "origin unresolved (HSC or MPP)", lg_held = "held in source",
    lg_ref = "reference relation (assay-defined end)", lg_product = "product formation", lg_member = "gate membership (not ancestry)",
    lg_rna = "RNA state candidate", lg_split = "HSC / MPP candidate", lg_assay = "assay-defined (reference)", lg_prodnode = "product", lg_heldnode = "held candidate",
    lg_taxnode = "taxonomy (reference)", lg_badge = "protein / function-defined (badge)",
    ref_panel = "Reference panel — 7 assay-defined · 4 taxonomy (separated from RNA candidates, not deleted)",
    ref_unres = "origin unresolved", tax_children = "children", membership = "membership (not ancestry)",
    lin_title = "Lineage detail — %s (pending user review)",
    lin_sub = "Grey column = assay-defined reference, teal column = taxonomy (children listed). Second line = source label. Orange mark = protein / function-defined",
    ref_col = "reference: assay-defined", tax_col = "reference: taxonomy",
    tab_head = c("edge", "source relation", "treatment", "support / adult", "evidence IDs"),
    tab_title = "%d source edges in this figure (source support grades unchanged)",
    cmp_title = "Source v1.1.0 versus proposal — what changes (pending user review)",
    cmp_nodes = "76 nodes → 77 proposed items", cmp_edges = "89 edges → treatment",
    cmp_left_nodes = "source node_kind", cmp_right_nodes = "proposed role", cmp_left_edges = "source support_class", cmp_right_edges = "proposed treatment",
    role = c(split = "HSC / MPP candidates (split)", rna_state_candidate = "RNA state candidates (kept)", assay_defined_reference = "assay-defined (reference)",
             product_reference = "product", held_candidate = "held", taxonomy_reference = "taxonomy (reference)"),
    split_line = "split decided · edges pending",
    hspc_left = "Source: one HSPC_continuum with 19 coarse relations", hspc_right = "Proposal: HSC_candidate and MPP_candidate kept separate; the origin of the 19 relations stays unresolved",
    hspc_note = "not copied · no HSC→MPP edge (absent in source) · markers and edges need separate approval",
    footer = "RNA distinction not assessed (not_assessed; HSC and MPP unresolved). Source layer, adult_support and evidence are unchanged. Tables: data/priors/rna_dag_review/proposal/"))


# ── 장치 · 글꼴 ──────────────────────────────────────────────────────────────
.rdr_device_kind <- function() if (isTRUE(capabilities("aqua"))) "quartz" else "cairo"

.rdr_family <- function(kind) {
  if (kind == "quartz") {
    ok <- tryCatch({ grDevices::quartzFonts(rdrkorean = grDevices::quartzFont(c("AppleSDGothicNeo-Regular", "AppleSDGothicNeo-Bold",
                                                                               "AppleSDGothicNeo-Regular", "AppleSDGothicNeo-Bold"))); TRUE },
                   error = function(e) FALSE)
    if (ok) return("rdrkorean")
  }
  "sans"
}

#' 같은 그림을 PNG 와 PDF 로 그린다. draw() 는 인치 좌표 (1 단위 = 1 인치, 원점 왼쪽 아래) 의 한 panel 에 그린다.
.rdr_save <- function(out_dir, stem, W, H, draw, kind, fam) {
  paths <- file.path(out_dir, paste0(stem, c(".png", ".pdf")))
  for (p in paths) {
    if (endsWith(p, ".png")) grDevices::png(p, width = W, height = H, units = "in", res = 150, pointsize = 10, type = kind)
    else if (kind == "quartz") grDevices::quartz(type = "pdf", file = p, width = W, height = H, pointsize = 10)
    else grDevices::cairo_pdf(p, width = W, height = H, pointsize = 10)
    tryCatch({
      graphics::par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), family = fam, xpd = NA)
      graphics::plot.new(); graphics::plot.window(c(0, W), c(0, H), xaxs = "i", yaxs = "i")
      draw()
    }, finally = grDevices::dev.off())
  }
  paths
}


# ── 기본 도형 ────────────────────────────────────────────────────────────────
.rdr_wrap_id <- function(id, width) {
  if (nchar(id) <= width) return(id)
  p <- strsplit(id, "_", fixed = TRUE)[[1]]; out <- character(0); cur <- ""
  for (i in seq_along(p)) {
    piece <- if (i < length(p)) paste0(p[i], "_") else p[i]
    if (nchar(cur) && nchar(cur) + nchar(piece) > width) { out <- c(out, cur); cur <- piece } else cur <- paste0(cur, piece)
  }
  c(out, cur)
}

.rdr_wrap_text <- function(s, width, max_lines) {
  w <- strsplit(s, " ", fixed = TRUE)[[1]]; out <- character(0); cur <- ""
  for (x in w) { if (nchar(cur) && nchar(cur) + 1 + nchar(x) > width) { out <- c(out, cur); cur <- x } else cur <- if (nchar(cur)) paste(cur, x) else x }
  out <- c(out, cur)
  if (length(out) > max_lines) { out <- out[seq_len(max_lines)]; out[max_lines] <- paste0(substr(out[max_lines], 1, max(1, width - 3)), "...") }
  out
}

.rdr_box <- function(x, y, w, h, style, lines, cex, fonts = 1, badge = FALSE) {
  st <- RDR_NODE_STYLE[[style]]
  graphics::rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = st$fill, border = st$border, lty = st$lty, lwd = st$lwd)
  n <- length(lines); fonts <- rep_len(fonts, n); cexs <- rep_len(cex, n); lh <- 0.13 * max(cexs) / 0.62
  for (i in seq_len(n)) graphics::text(x, y + ((n + 1) / 2 - i) * lh, lines[i], cex = cexs[i], font = fonts[i])
  if (isTRUE(badge)) graphics::rect(x + w / 2 - 0.09, y + h / 2 - 0.09, x + w / 2 - 0.015, y + h / 2 - 0.015, col = RDR_BADGE_COL, border = NA)
}

.rdr_tag <- function(x, y, s, cex, col = "#333333", bg = "#FFFFFF", adj = 0.5, font = 1) {
  w <- graphics::strwidth(s, cex = cex, units = "user", font = font) + 0.04; h <- graphics::strheight("Ag", cex = cex, units = "user") + 0.03
  x0 <- x - adj * w + 0.02 * (adj == 0)
  graphics::rect(x0 - 0.02, y - h / 2, x0 + w - 0.02, y + h / 2, col = bg, border = NA)
  graphics::text(x, y, s, cex = cex, col = col, adj = c(adj, 0.5), font = font)
}

.rdr_exit <- function(cx, cy, hw, hh, dx, dy) {
  t <- min(if (abs(dx) > 1e-9) hw / abs(dx) else Inf, if (abs(dy) > 1e-9) hh / abs(dy) else Inf)
  c(cx + t * dx, cy + t * dy)
}

.rdr_bezier <- function(p0, p1, p2, n = 40) { t <- seq(0, 1, length.out = n); cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
                                                                                    (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2]) }

#' 꺾은 선 꼭짓점 -> 촘촘한 점 (충돌 검사 · 끝 화살표용)
.rdr_dense <- function(p, step = 0.03) {
  out <- p[1, , drop = FALSE]
  for (i in seq_len(nrow(p))[-1]) { d <- sqrt(sum((p[i, ] - p[i - 1, ])^2)); n <- max(2L, ceiling(d / step) + 1L); t <- seq(0, 1, length.out = n)[-1]
    out <- rbind(out, cbind(p[i - 1, 1] + (p[i, 1] - p[i - 1, 1]) * t, p[i - 1, 2] + (p[i, 2] - p[i - 1, 2]) * t)) }
  out
}

.rdr_hit_count <- function(p, boxes, skip, pad = 0.03) {
  q <- .rdr_dense(p); q <- q[-c(1, nrow(q)), , drop = FALSE]; b <- boxes[!boxes$id %in% skip, , drop = FALSE]
  if (!nrow(b) || !nrow(q)) return(0L)
  sum(vapply(seq_len(nrow(b)), function(i) any(abs(q[, 1] - b$x[i]) < b$w[i] / 2 + pad & abs(q[, 2] - b$y[i]) < b$h[i] / 2 + pad), NA))
}

.rdr_label_xy <- function(p) {
  if (nrow(p) > 6) return(p[ceiling(nrow(p) / 2), ])
  len <- sqrt(rowSums((p[-1, , drop = FALSE] - p[-nrow(p), , drop = FALSE])^2)); i <- which.max(len); (p[i, ] + p[i + 1, ]) / 2
}

#' 상자 a -> b 경로. 후보를 순서대로 보고 다른 상자를 지나지 않는 첫 후보를 고른다 (없으면 가장 덜 겹치는 후보). attr(, "kind") = 고른 후보.
#'   main    : 직선 -> 층 사이 틈으로 꺾은 선 -> 작은 곡선 -> 같은 열을 건너뛰는 연결용 옆 통로 (오른쪽 · 왼쪽) -> 큰 곡선
#'   ref_main: 직선 -> 참고 열 옆 통로로 꺾은 선 -> 곡선
#'   ref_ref : 직선 -> 참고 열 안 오른쪽 옆 통로 -> 곡선
#' side_k 는 같은 열에서 옆 통로를 쓰는 연결끼리 겹치지 않게 띄우는 순번이다.
.rdr_route <- function(a, b, kind, boxes, gap_mid, yg_off = 0, lane_x = NULL, side_k = 1L) {
  g <- 0.02; c0 <- c(a$x, a$y); c1 <- c(b$x, b$y); d <- c1 - c0; len <- sqrt(sum(d^2))
  straight <- rbind(.rdr_exit(a$x, a$y, a$w / 2 + g, a$h / 2 + g, d[1], d[2]), .rdr_exit(b$x, b$y, b$w / 2 + g, b$h / 2 + g, -d[1], -d[2]))
  arc <- function(off) { nrm <- c(-d[2], d[1]) / len; ctl <- (c0 + c1) / 2 + off * nrm
    .rdr_bezier(.rdr_exit(a$x, a$y, a$w / 2 + g, a$h / 2 + g, ctl[1] - a$x, ctl[2] - a$y), ctl, .rdr_exit(b$x, b$y, b$w / 2 + g, b$h / 2 + g, ctl[1] - b$x, ctl[2] - b$y)) }
  side <- function(sgn, gap) { dy <- (side_k - 1.5) * 0.07
    xe <- if (sgn > 0) max(a$x + a$w / 2, b$x + b$w / 2) + gap + 0.08 * (side_k - 1) else min(a$x - a$w / 2, b$x - b$w / 2) - gap - 0.08 * (side_k - 1)
    rbind(c(a$x + sgn * (a$w / 2 + g), a$y + dy), c(xe, a$y + dy), c(xe, b$y + dy), c(b$x + sgn * (b$w / 2 + g), b$y + dy)) }
  b_top <- b$y + b$h / 2 + g; yg <- b$y + b$h / 2 + gap_mid + yg_off
  ortho <- if (a$y - a$h / 2 - g > yg) rbind(c(a$x, a$y - a$h / 2 - g), c(a$x, yg), c(b$x, yg), c(b$x, b_top))
  lane <- if (!is.null(lane_x)) rbind(c(a$x + a$w / 2 + g, a$y), c(lane_x, a$y), c(lane_x, yg), c(b$x, yg), c(b$x, b_top))
  cand <- switch(kind,
    main = list(straight = straight, ortho = ortho, arc = arc(0.45), arc = arc(-0.45), side = side(1, 0.12), side = side(-1, 0.12),
                arc = arc(0.8), arc = arc(-0.8), arc = arc(1.2), arc = arc(-1.2)),
    ref_main = list(straight = straight, lane = lane, arc = arc(0.45), arc = arc(-0.45)),
    ref_ref = list(straight = straight, side = side(1, 0.06), arc = arc(0.45), arc = arc(-0.45), arc = arc(0.8), arc = arc(-0.8)))
  cand <- cand[!vapply(cand, is.null, NA)]
  wide <- vapply(cand, .rdr_hit_count, 0L, boxes = boxes, skip = c(a$id, b$id), pad = 0.08)    # 다른 상자 가장자리에 붙어 지나가는 경로도 먼저 피한다
  tight <- vapply(cand, .rdr_hit_count, 0L, boxes = boxes, skip = c(a$id, b$id), pad = 0.03)
  i <- if (any(wide == 0L)) which(wide == 0L)[1] else which(tight == min(tight))[1]; p <- cand[[i]]; attr(p, "kind") <- names(cand)[i]; p
}

.rdr_draw_path <- function(p, style) {
  st <- RDR_EDGE_STYLE[[style]]; q <- .rdr_dense(p, 0.05); m <- nrow(q)
  graphics::lines(q[1:(m - 1), 1], q[1:(m - 1), 2], col = st$col, lty = st$lty, lwd = st$lwd)
  if (st$arrow) graphics::arrows(q[m - 1, 1], q[m - 1, 2], q[m, 1], q[m, 2], length = 0.055, angle = 22, col = st$col, lwd = st$lwd, lty = 1)
  else graphics::lines(q[(m - 1):m, 1], q[(m - 1):m, 2], col = st$col, lty = st$lty, lwd = st$lwd)
}


# ── 배치 ─────────────────────────────────────────────────────────────────────
#' 층 배치: 층 = 묶음 안 최장 경로 깊이 (1 부터), x = 부모 평균 (첫 층은 원본 layout_x 순서), 최소 간격 1 slot.
rdr_layer_layout <- function(ids, from, to, ord) {
  layer <- stats::setNames(rep(1L, length(ids)), ids)
  if (length(from)) {
    g <- igraph::graph_from_data_frame(data.frame(from, to), vertices = data.frame(id = ids))
    for (v in names(igraph::topo_sort(g, mode = "out"))) for (k in to[from == v]) layer[k] <- max(layer[k], layer[v] + 1L)
  }
  x <- stats::setNames(numeric(length(ids)), ids)
  for (L in sort(unique(layer))) {
    v <- names(layer)[layer == L]
    key <- if (L == 1L) rank(ord[v], ties.method = "first") - 1 else vapply(v, function(z) mean(x[from[to == z & from %in% names(layer)[layer < L]]]), 0)
    o <- order(key, ord[v]); v <- v[o]; key <- key[o]; pos <- key
    for (i in seq_along(pos)[-1]) pos[i] <- max(pos[i], pos[i - 1] + 1)
    x[v] <- pos - mean(pos) + mean(key)
  }
  data.frame(id = ids, layer = unname(layer[ids]), x = unname(x[ids] - min(x)), stringsAsFactors = FALSE)
}

.rdr_role <- function(rv) stats::setNames(rv$proposed_nodes$view_role, rv$proposed_nodes$proposed_id)

.rdr_node_style <- function(id, role) if (id %in% RDR_SPLIT_IDS) "split" else role[[id]]

.rdr_edge_style <- function(em) ifelse(em$treatment == "keep", "keep", ifelse(em$treatment == "unresolved", "unresolved", ifelse(em$treatment == "held", "held",
                                ifelse(em$view_panel == "product_reference", "product", "reference"))))

#' 계통 묶음의 노드 · 연결 (그림과 테스트가 같이 쓴다). 반환 list(main, ref, tax, edges)
rdr_lineage_view <- function(rv, src, group) {
  gdef <- RDR_LINEAGE_GROUPS[[group]]; pn <- rv$proposed_nodes; em <- rv$source_edge_map
  main <- pn$proposed_id[pn$module %in% gdef$modules & pn$view_role %in% c("rna_state_candidate", "product_reference", "held_candidate")]
  assay <- pn$proposed_id[pn$view_role == "assay_defined_reference"]
  gm <- src$gate_membership
  ref <- unique(c(intersect(c(em$source[em$target %in% main], em$target[em$source %in% main]), assay), gm$gate_node[gm$subset_node %in% main]))
  repeat { up <- intersect(em$source[em$target %in% ref], assay); if (all(up %in% ref)) break; ref <- unique(c(ref, up)) }
  tx <- src$taxonomy; tax <- unique(tx$parent[tx$child %in% main])
  shown <- c(main, ref, RDR_SPLIT_SOURCE)
  e <- em[em$source %in% shown & em$target %in% shown, ]
  edges <- data.frame(source_edge_id = e$source_edge_id, from = ifelse(e$source == RDR_SPLIT_SOURCE, "origin", e$source), to = e$target,
                      style = .rdr_edge_style(e), treatment = e$treatment, stringsAsFactors = FALSE)
  if (length(tax)) { t2 <- tx[tx$parent %in% tax & tx$child %in% main, ]
    edges <- rbind(edges, data.frame(source_edge_id = paste0("is_a:", t2$child), from = t2$child, to = t2$parent, style = "taxonomy", treatment = "taxonomy")) }
  g2 <- gm[gm$subset_node %in% main & gm$gate_node %in% ref, ]
  if (nrow(g2)) edges <- rbind(edges, data.frame(source_edge_id = g2$membership_id, from = g2$subset_node, to = g2$gate_node, style = "membership", treatment = "membership"))
  list(main = main, ref = ref, tax = tax, edges = edges)
}

#' 그리는 노드의 원본 layout 좌표만 수치로 바꾼다 (taxonomy 의 "NA" 좌표는 읽지 않는다). 유한하지 않으면 멈춘다.
.rdr_group_layout <- function(src, view) {
  n <- src$nodes[match(view$main, src$nodes$node_id), ]
  ord <- as.numeric(n$layout_x) + as.numeric(n$layout_y) / 100
  if (any(!is.finite(ord))) stop("layout_x / layout_y must be finite for drawn nodes: ", paste(view$main[!is.finite(ord)], collapse = ", "))
  names(ord) <- view$main
  e <- view$edges[view$edges$from %in% view$main & view$edges$to %in% view$main & !view$edges$style %in% c("taxonomy", "membership"), ]
  rdr_layer_layout(view$main, e$from, e$to, ord)
}

.rdr_boxes <- function(id, x, y, w, h, style, badge = FALSE) {
  n <- length(id)
  if (!n) return(data.frame(id = character(0), x = numeric(0), y = numeric(0), w = numeric(0), h = numeric(0), style = character(0), badge = logical(0), stringsAsFactors = FALSE))
  data.frame(id = id, x = rep_len(x, n), y = rep_len(y, n), w = rep_len(w, n), h = rep_len(h, n), style = rep_len(style, n), badge = rep_len(badge, n), stringsAsFactors = FALSE)
}

.rdr_legend_bbox <- function(x, ytop) c(x0 = x, x1 = x + sum(RDR_LEGEND_COLW), y0 = ytop - 3 * RDR_LEGEND_ROWH - 0.1, y1 = ytop + 0.1)

.rdr_draw_legend <- function(x, ytop, tx, cex = 0.58) {
  cols <- list(list(c("n", "rna_state_candidate", tx$lg_rna), c("n", "split", tx$lg_split), c("n", "assay_defined_reference", tx$lg_assay), c("n", "product_reference", tx$lg_prodnode)),
               list(c("n", "held_candidate", tx$lg_heldnode), c("n", "taxonomy_reference", tx$lg_taxnode), c("b", "", tx$lg_badge)),
               list(c("e", "keep", tx$lg_keep), c("e", "unresolved", tx$lg_unres), c("e", "held", tx$lg_held)),
               list(c("e", "reference", tx$lg_ref), c("e", "product", tx$lg_product), c("e", "membership", tx$lg_member)))
  xs <- x + c(0, cumsum(RDR_LEGEND_COLW)[-length(RDR_LEGEND_COLW)])
  for (j in seq_along(cols)) for (i in seq_along(cols[[j]])) {
    it <- cols[[j]][[i]]; yy <- ytop - (i - 1) * RDR_LEGEND_ROWH; xx <- xs[j]
    if (it[1] == "n") { .rdr_box(xx + 0.16, yy, 0.3, 0.13, it[2], "", cex); graphics::text(xx + 0.4, yy, it[3], adj = c(0, 0.5), cex = cex) }
    else if (it[1] == "b") { graphics::rect(xx + 0.12, yy - 0.04, xx + 0.2, yy + 0.04, col = RDR_BADGE_COL, border = NA); graphics::text(xx + 0.4, yy, it[3], adj = c(0, 0.5), cex = cex) }
    else { st <- RDR_EDGE_STYLE[[it[2]]]
      if (st$arrow) graphics::arrows(xx, yy, xx + 0.42, yy, length = 0.05, col = st$col, lty = st$lty, lwd = st$lwd) else graphics::segments(xx, yy, xx + 0.42, yy, col = st$col, lty = st$lty, lwd = st$lwd)
      graphics::text(xx + 0.52, yy, it[3], adj = c(0, 0.5), cex = cex) }
  }
}


# ── 핵심안 배치 ──────────────────────────────────────────────────────────────
rdr_core_spec <- function(rv, src) {
  K <- RDR_CORE; pn <- rv$proposed_nodes; em <- rv$source_edge_map; role <- .rdr_role(rv)
  badge <- stats::setNames(grepl("protein_or_function_defined", pn$badges, fixed = TRUE), pn$proposed_id)
  views <- lapply(names(RDR_LINEAGE_GROUPS), function(g) rdr_lineage_view(rv, src, g)); names(views) <- names(RDR_LINEAGE_GROUPS)
  lays <- lapply(views, function(v) .rdr_group_layout(src, v))
  widths <- vapply(lays, function(l) max(l$x) + 1, 0); nl <- max(vapply(lays, function(l) max(l$layer), 0))
  ref_ids <- pn$proposed_id[pn$view_role == "assay_defined_reference"]; tax_ids <- pn$proposed_id[pn$view_role == "taxonomy_reference"]
  n_lines <- vapply(c(ref_ids, tax_ids), function(id) if (id %in% ref_ids) sum(em$source == id | em$target == id) + sum(src$gate_membership$gate_node == id)
                    else 1L + sum(src$taxonomy$parent == id), 0)
  W <- 1.2 + sum(widths) * K$slot + K$gap_g * (length(widths) - 1)
  H_graph <- K$top + 0.55 + (nl - 1) * K$dy + K$bh / 2 + 0.45; ref_h <- 1.4 + max(n_lines) * 0.15; leg_h <- 1.3; H <- H_graph + ref_h + leg_h
  ytop <- H - K$top - 0.55; xoff <- 0.6; boxes <- .rdr_boxes(character(0), numeric(0), numeric(0), numeric(0), numeric(0), character(0), logical(0)); groups <- list()
  for (g in names(lays)) {
    l <- lays[[g]]; cx <- xoff + (l$x + 0.5) * K$slot; cy <- ytop - (l$layer - 1) * K$dy
    groups[[g]] <- c(x0 = xoff - 0.08, x1 = xoff + widths[[g]] * K$slot + 0.08, y0 = ytop - (max(l$layer) - 1) * K$dy - K$bh / 2 - 0.13, y1 = ytop + 0.42,
                     hx = xoff + widths[[g]] * K$slot / 2, hy = ytop + 0.3)
    boxes <- rbind(boxes, .rdr_boxes(l$id, cx, cy, K$bw, K$bh, vapply(l$id, .rdr_node_style, "", role = role, USE.NAMES = FALSE), unname(badge[l$id])))
    xoff <- xoff + widths[[g]] * K$slot + K$gap_g
  }
  ox <- W / 2; oy <- H - K$top + 0.35
  bracket <- c(x0 = ox - 1.75, x1 = ox + 1.75, y0 = oy - 0.27, y1 = oy + 0.33, xc = ox)
  cands <- .rdr_boxes(RDR_SPLIT_IDS, ox + c(-0.85, 0.85), oy + 0.03, 1.5, 0.36, "split")
  routes <- list(); bus <- list(); bus_target <- list(); bus_y <- bracket[["y0"]] - 0.12; n_side <- 0L
  for (g in names(views)) { e <- views[[g]]$edges
    for (i in seq_len(nrow(e))) {
      if (e$style[i] %in% c("taxonomy", "membership") || !e$to[i] %in% boxes$id) next
      b <- boxes[boxes$id == e$to[i], ]
      if (e$from[i] == "origin") {
        bus[[e$source_edge_id[i]]] <- rbind(c(ox, bracket[["y0"]]), c(ox, bus_y), c(b$x, bus_y), c(b$x, b$y + b$h / 2 + 0.02)); bus_target[[e$source_edge_id[i]]] <- b$id; next }
      if (!e$from[i] %in% boxes$id) next
      a <- boxes[boxes$id == e$from[i], ]
      p <- .rdr_route(a, b, "main", boxes, (K$dy - K$bh) / 2, side_k = n_side + 1L); n_side <- n_side + (attr(p, "kind") == "side")
      routes[[length(routes) + 1L]] <- list(id = e$source_edge_id[i], from = a$id, to = b$id, style = e$style[i], pts = p)
    } }
  py <- H - H_graph - 0.25; cw <- (W - 1.2) / (length(ref_ids) + length(tax_ids)); bwr <- min(cw - 0.15, 1.45)
  ref_boxes <- .rdr_boxes(c(ref_ids, tax_ids), 0.6 + (seq_along(c(ref_ids, tax_ids)) - 0.5) * cw, py - 0.62, bwr, 0.34, unname(role[c(ref_ids, tax_ids)]))
  leg_top <- 0.3 + 0.1 + 3 * RDR_LEGEND_ROWH
  list(W = W, H = H, groups = groups, boxes = rbind(boxes, ref_boxes), graph_ids = boxes$id, candidates = cands, bracket = bracket, bus = bus, bus_target = bus_target,
       routes = routes, ref_panel = c(x0 = 0.3, x1 = W - 0.3, y0 = py - ref_h + 0.15, y1 = py + 0.1, title_y = py - 0.08), ref_ids = ref_ids, tax_ids = tax_ids,
       legend_xy = c(0.5, leg_top), legend_bbox = .rdr_legend_bbox(0.5, leg_top))
}


# ── 계통별 배치 ──────────────────────────────────────────────────────────────
rdr_lineage_spec <- function(rv, src, group) {
  K <- RDR_LIN; pn <- rv$proposed_nodes; em <- rv$source_edge_map; role <- .rdr_role(rv)
  badge <- stats::setNames(grepl("protein_or_function_defined", pn$badges, fixed = TRUE), pn$proposed_id)
  v <- rdr_lineage_view(rv, src, group); l <- .rdr_group_layout(src, v); e <- v$edges; ly <- stats::setNames(l$layer, l$id)
  # 참고 열 층: 본체 도착점 층 - 0.5, 참고끼리는 도착점 - 1, gate 포함 관계는 하위 노드 층
  ry <- stats::setNames(rep(NA_real_, length(v$ref)), v$ref)
  for (it in 1:6) for (r in v$ref) {
    tg <- e$to[e$from == r]; s <- e$from[e$to == r & e$style == "membership"]
    val <- c(ly[intersect(tg, names(ly))] - 0.5, ry[intersect(tg, v$ref)] - 1, ly[intersect(s, names(ly))]); val <- val[!is.na(val)]
    if (length(val)) ry[r] <- min(val)
  }
  ry[is.na(ry)] <- 1; if (length(ry)) { o <- order(ry); rr <- ry[o]; for (i in seq_along(rr)[-1]) rr[i] <- max(rr[i], rr[i - 1] + 0.85); ry[names(rr)] <- rr }
  kids <- lapply(stats::setNames(v$tax, v$tax), function(t) e$from[e$to == t & e$style == "taxonomy"])
  ty <- vapply(v$tax, function(t) mean(ly[kids[[t]]]), 0)
  if (length(ty)) { o <- order(ty); tt <- ty[o]; for (i in seq_along(tt)[-1]) tt[i] <- max(tt[i], tt[i - 1] + 1.6); ty[names(tt)] <- tt }
  tax_drop <- vapply(v$tax, function(t) (K$bh / 2 + length(kids[[t]]) * 0.14 + 0.15) / K$dy, 0)
  top_layer <- min(c(1, ry)) - 1.4; bottom_layer <- max(c(ly, ry, ty + tax_drop))
  ref_w <- if (length(v$ref)) K$ref_w else 0.3; tax_w <- if (length(v$tax)) K$tax_w else 0; main_w <- (max(l$x) + 1) * K$slot
  need <- 0.15 + ref_w + 0.25 + main_w + 0.25 + tax_w + 0.3; gw <- max(need, sum(RDR_LEGEND_COLW) + 0.9)
  tab <- em[em$source_edge_id %in% e$source_edge_id, ]; tab <- tab[order(tab$source_edge_id), ]
  extra <- e[e$style %in% c("taxonomy", "membership"), ]
  gh <- 1.5 + (bottom_layer - top_layer) * K$dy + 0.55 + 1.35; th <- 1.25 + (nrow(tab) + nrow(extra) + 1) * 0.19 + 0.8
  H <- max(gh, th); W <- gw + K$tab_w
  Y <- function(layer) H - 1.5 - (layer - top_layer) * K$dy
  x_ref <- 0.15 + 0.3 + K$bw / 2; mx <- 0.15 + ref_w + 0.25 + max(0, gw - need) / 2 + (l$x + 0.5) * K$slot; x_tax <- gw - 0.3 - tax_w / 2
  boxes <- rbind(.rdr_boxes(l$id, mx, Y(l$layer), K$bw, K$bh, vapply(l$id, .rdr_node_style, "", role = role, USE.NAMES = FALSE), unname(badge[l$id])),
                 .rdr_boxes(v$ref, rep(x_ref, length(v$ref)), Y(ry[v$ref]), K$bw, K$bh, rep("assay_defined_reference", length(v$ref))),
                 .rdr_boxes(v$tax, rep(x_tax, length(v$tax)), Y(ty[v$tax]), K$bw, K$bh, rep("taxonomy_reference", length(v$tax))))
  ox <- mean(range(mx)); oy <- Y(top_layer)
  bracket <- c(x0 = ox - 1.6, x1 = ox + 1.6, y0 = oy - 0.3, y1 = oy + 0.3, xc = ox)
  cands <- .rdr_boxes(RDR_SPLIT_IDS, ox + c(-0.78, 0.78), oy, 1.4, 0.4, "split")
  B <- function(id) boxes[boxes$id == id, ]
  rm_edges <- e[e$from %in% v$ref & e$to %in% l$id | e$style == "membership", ]
  lane_lo <- x_ref + K$bw / 2 + 0.22; lane_hi <- min(mx) - K$bw / 2 - 0.1
  lanes <- if (nrow(rm_edges)) stats::setNames(seq(lane_lo, lane_hi, length.out = nrow(rm_edges) + 2)[-c(1, nrow(rm_edges) + 2)], rm_edges$source_edge_id) else numeric(0)
  routes <- list(); n_into <- list(); offs <- c(0, 0.05, -0.05, 0.1, -0.1, 0.03, -0.03); n_side <- 0L
  for (i in seq_len(nrow(e))) {
    if (e$style[i] == "taxonomy" || e$from[i] == "origin") next
    memb <- e$style[i] == "membership"; a <- B(if (memb) e$to[i] else e$from[i]); b <- B(if (memb) e$from[i] else e$to[i])
    kind <- if (a$id %in% v$ref && b$id %in% v$ref) "ref_ref" else if (a$id %in% v$ref) "ref_main" else "main"
    key <- as.character(round(b$y, 3)); n_into[[key]] <- (if (is.null(n_into[[key]])) 0L else n_into[[key]]) + 1L
    p <- .rdr_route(a, b, kind, boxes, (K$dy - K$bh) / 2, offs[(n_into[[key]] - 1L) %% length(offs) + 1L], if (kind == "ref_main") lanes[[e$source_edge_id[i]]] else NULL,
                    side_k = n_side + 1L)
    n_side <- n_side + (attr(p, "kind") == "side")
    routes[[length(routes) + 1L]] <- list(id = e$source_edge_id[i], from = a$id, to = b$id, style = e$style[i], pts = p, label_xy = .rdr_label_xy(p))
  }
  bus <- list(); bus_target <- list(); bus_lab <- list(); bus_y <- bracket[["y0"]] - 0.14; k_ref <- 0L
  for (i in which(e$from == "origin")) {
    b <- B(e$to[i]); id <- e$source_edge_id[i]
    if (b$id %in% v$ref) { k_ref <- k_ref + 1L; xl <- x_ref - K$bw / 2 - 0.02 - 0.07 * k_ref
      bus[[id]] <- rbind(c(ox, bracket[["y0"]]), c(ox, bus_y), c(xl, bus_y), c(xl, b$y), c(x_ref - K$bw / 2 - 0.02, b$y)); bus_lab[[id]] <- c(xl + 0.03, b$y + 0.1)
    } else { bus[[id]] <- rbind(c(ox, bracket[["y0"]]), c(ox, bus_y), c(b$x, bus_y), c(b$x, b$y + b$h / 2 + 0.02)); bus_lab[[id]] <- c(b$x + 0.04, bus_y - 0.16) }   # 이름표는 위쪽 (층 사이 틈의 참고 연결과 겹치지 않게)
    bus_target[[id]] <- b$id
  }
  leg_top <- 0.3 + 0.1 + 3 * RDR_LEGEND_ROWH
  list(W = W, H = H, gw = gw, boxes = boxes, graph_ids = boxes$id, candidates = cands, bracket = bracket, bus = bus, bus_target = bus_target, bus_lab = bus_lab, routes = routes,
       ref_rect = if (length(v$ref)) c(x0 = 0.15, x1 = 0.15 + ref_w, y0 = Y(bottom_layer) - 0.55, y1 = Y(top_layer + 1) + 0.5) else NULL,
       tax_rect = if (length(v$tax)) c(x0 = gw - 0.3 - tax_w, x1 = gw - 0.3, y0 = Y(bottom_layer) - 0.55, y1 = Y(top_layer + 1) + 0.5) else NULL,
       tax_kids = kids, tab = tab, extra = extra, legend_xy = c(0.4, leg_top), legend_bbox = .rdr_legend_bbox(0.4, leg_top))
}


# ── 그림 1: 원본 대 제안 ──────────────────────────────────────────────────────
.rdr_flow_panel <- function(x0, y0, w, h, left, right, flows, title, lt, rt, cols, unit) {
  graphics::text(x0, y0 + h + 0.38, title, adj = c(0, 0), cex = 0.95, font = 2)
  graphics::text(x0, y0 + h + 0.12, lt, adj = c(0, 0), cex = 0.7, col = "#555555"); graphics::text(x0 + w, y0 + h + 0.12, rt, adj = c(1, 0), cex = 0.7, col = "#555555")
  bw <- 0.22; gapy <- 0.09
  seg <- function(v) { tot <- sum(v) * unit + gapy * (length(v) - 1); top <- y0 + h - (h - tot) / 2; s <- list()
    for (k in names(v)) { s[[k]] <- c(top - v[[k]] * unit, top); top <- top - v[[k]] * unit - gapy }; s }
  L <- seg(left); R <- seg(right); lu <- stats::setNames(rep(0, length(left)), names(left)); ru <- stats::setNames(rep(0, length(right)), names(right))
  xl <- x0 + 2.25; xr <- x0 + w - 2.25
  for (i in seq_len(nrow(flows))) {
    a <- flows$from[i]; b <- flows$to[i]; nl <- flows$n_left[i]; nr <- flows$n_right[i]
    la <- L[[a]][2] - lu[[a]] * unit; lb <- la - nl * unit; lu[[a]] <- lu[[a]] + nl
    ra <- R[[b]][2] - ru[[b]] * unit; rb <- ra - nr * unit; ru[[b]] <- ru[[b]] + nr
    t <- seq(0, 1, length.out = 50); s <- 3 * t^2 - 2 * t^3; xx <- xl + bw + (xr - xl - bw) * t
    graphics::polygon(c(xx, rev(xx)), c(la + (ra - la) * s, rev(lb + (rb - lb) * s)), col = grDevices::adjustcolor(cols[[b]], 0.45), border = NA)
  }
  for (k in names(left)) { graphics::rect(xl, L[[k]][1], xl + bw, L[[k]][2], col = "#666666", border = NA)
    graphics::text(xl - 0.06, mean(L[[k]]), sprintf("%s  %d", k, left[[k]]), adj = c(1, 0.5), cex = 0.66) }
  for (k in names(right)) { graphics::rect(xr, R[[k]][1], xr + bw, R[[k]][2], col = cols[[k]], border = NA)
    graphics::text(xr + bw + 0.06, mean(R[[k]]), sprintf("%s  %d", k, right[[k]]), adj = c(0, 0.5), cex = 0.66) }
}

rdr_draw_comparison <- function(rv, src, out_dir, tx, kind, fam) {
  nm <- rv$source_node_map; em <- rv$source_edge_map
  rk <- ifelse(nm$treatment == "split", "split", nm$view_role); rl <- tx$role[rk]
  kinds <- c("stem_state", "assay_gate", "precursor_state", "stage_state", "mature_state", "activation_state", "product", "held_candidate", "taxonomy_category")
  left_n <- stats::setNames(as.integer(table(factor(nm$source_node_kind, kinds))), kinds)
  rorder <- c("split", "rna_state_candidate", "assay_defined_reference", "product_reference", "held_candidate", "taxonomy_reference")
  right_n <- stats::setNames(vapply(rorder, function(r) if (r == "split") length(RDR_SPLIT_IDS) else if (r == "rna_state_candidate") sum(nm$treatment == "keep")
                                    else sum(rv$proposed_nodes$view_role == r), 0L), tx$role[rorder])
  fl <- stats::aggregate(list(n = rep(1L, nrow(nm))), list(from = nm$source_node_kind, to = unname(rl)), sum)
  fl$n_left <- fl$n; fl$n_right <- ifelse(fl$from == "stem_state", length(RDR_SPLIT_IDS), fl$n)
  fl <- fl[order(match(fl$from, kinds), match(fl$to, tx$role[rorder])), ]
  ncols <- stats::setNames(c("#B07D00", "#4D4D4D", "#8A8A8A", "#7A7A7A", "#C0392B", "#2A8C8C"), tx$role[rorder])
  sc <- c("coarse", "conditional", "supported", "held"); tr <- c("keep", "unresolved", "reference_only", "held")
  el <- stats::setNames(as.integer(table(factor(em$support_class, sc))), sc); er <- stats::setNames(as.integer(table(factor(em$treatment, tr))), tx$treat[tr])
  fe <- stats::aggregate(list(n = rep(1L, nrow(em))), list(from = em$support_class, to = unname(tx$treat[em$treatment])), sum)
  fe$n_left <- fe$n; fe$n_right <- fe$n; fe <- fe[order(match(fe$from, sc), match(fe$to, tx$treat[tr])), ]
  ecols <- stats::setNames(c("#333333", "#3B6FB6", "#9C9C9C", "#C0392B"), tx$treat[tr])
  hs <- em[em$treatment == "unresolved", ]; W <- 16.2; H <- 10.2
  draw <- function() {
    graphics::text(0.4, H - 0.35, tx$cmp_title, adj = c(0, 1), cex = 1.25, font = 2)
    .rdr_flow_panel(0.4, H - 5.3, 7.4, 4.0, left_n, right_n, fl, tx$cmp_nodes, tx$cmp_left_nodes, tx$cmp_right_nodes, ncols, 0.045)
    .rdr_flow_panel(8.4, H - 5.3, 7.4, 4.0, el, er, fe, tx$cmp_edges, tx$cmp_left_edges, tx$cmp_right_edges, ecols, 0.04)
    y_top <- 3.75; tgt <- hs$target; nt <- length(tgt); xs <- seq(0.7, 7.3, length.out = nt)
    graphics::text(0.4, y_top + 0.45, tx$hspc_left, adj = c(0, 0), cex = 0.85, font = 2)
    .rdr_box(4, y_top - 0.05, 2.0, 0.42, "split", c("HSPC_continuum", "(v1.1.0)"), c(0.68, 0.55), c(2, 1))
    for (i in seq_len(nt)) { graphics::segments(4, y_top - 0.26, xs[i], 1.62, col = "#555555", lwd = 0.8); graphics::points(xs[i], 1.55, pch = 16, cex = 0.5, col = "#555555")
      graphics::text(xs[i], 1.45, paste0(hs$source_edge_id[i], " ", tgt[i]), srt = 90, adj = c(1, 0.5), cex = 0.48) }
    x2 <- 8.4; xs2 <- seq(x2 + 0.3, x2 + 6.9, length.out = nt)
    graphics::text(x2, y_top + 0.45, tx$hspc_right, adj = c(0, 0), cex = 0.85, font = 2)
    graphics::rect(x2 + 1.95, y_top - 0.38, x2 + 5.25, y_top + 0.28, border = RDR_EDGE_STYLE$unresolved$col, lty = 2, lwd = 1.1)
    .rdr_box(x2 + 2.8, y_top - 0.05, 1.45, 0.42, "split", c("HSC_candidate", tx$split_line), c(0.66, 0.5), c(2, 1))
    .rdr_box(x2 + 4.4, y_top - 0.05, 1.45, 0.42, "split", c("MPP_candidate", tx$split_line), c(0.66, 0.5), c(2, 1))
    for (i in seq_len(nt)) { graphics::segments(x2 + 3.6, y_top - 0.38, xs2[i], 1.62, col = RDR_EDGE_STYLE$unresolved$col, lty = 2, lwd = 0.8)
      graphics::points(xs2[i], 1.55, pch = 16, cex = 0.5, col = RDR_EDGE_STYLE$unresolved$col)
      graphics::text(xs2[i], 1.45, paste0(hs$source_edge_id[i], " ", tgt[i]), srt = 90, adj = c(1, 0.5), cex = 0.48) }
    graphics::text(x2, 0.3, tx$hspc_note, adj = c(0, 0), cex = 0.7, col = RDR_EDGE_STYLE$unresolved$col)
    graphics::text(0.4, 0.3, tx$footer, adj = c(0, 0), cex = 0.6, col = "#555555")
  }
  .rdr_save(out_dir, "rna_dag_original_vs_proposed", W, H, draw, kind, fam)
}


# ── 그림 2: 핵심안 ───────────────────────────────────────────────────────────
rdr_draw_core <- function(rv, src, out_dir, tx, kind, fam) {
  s <- rdr_core_spec(rv, src); K <- RDR_CORE; pn <- rv$proposed_nodes; em <- rv$source_edge_map; W <- s$W; H <- s$H
  r <- table(factor(pn$view_role, c("rna_state_candidate", "assay_defined_reference", "product_reference", "held_candidate", "taxonomy_reference")))
  draw <- function() {
    graphics::text(0.4, H - 0.3, tx$core_title, adj = c(0, 1), cex = 1.25, font = 2)
    graphics::text(0.4, H - 0.62, sprintf(tx$core_sub, r[[1]], r[[2]], r[[3]], r[[4]], r[[5]]), adj = c(0, 1), cex = 0.72, col = "#444444")
    graphics::text(0.4, H - 0.86, tx$layer_note, adj = c(0, 1), cex = 0.66, col = "#444444")
    for (g in names(s$groups)) { q <- s$groups[[g]]; graphics::rect(q[["x0"]], q[["y0"]], q[["x1"]], q[["y1"]], col = "#F7F7F7", border = NA) }
    b <- s$bracket; graphics::rect(b[["x0"]], b[["y0"]], b[["x1"]], b[["y1"]], border = RDR_EDGE_STYLE$unresolved$col, lty = 2, lwd = 1.1)
    graphics::text(b[["xc"]], b[["y1"]] + 0.1, tx$origin, cex = 0.68, font = 2, col = "#7A5A00")
    graphics::text(b[["x1"]] + 0.15, (b[["y0"]] + b[["y1"]]) / 2, tx$origin_note, adj = c(0, 0.5), cex = 0.58, col = RDR_EDGE_STYLE$unresolved$col)
    cd <- s$candidates; for (i in seq_len(nrow(cd))) .rdr_box(cd$x[i], cd$y[i], cd$w[i], cd$h[i], "split", cd$id[i], 0.7, 2)
    for (p in s$bus) .rdr_draw_path(p, "unresolved")
    for (q in s$routes) .rdr_draw_path(q$pts, q$style)
    gb <- s$boxes[s$boxes$id %in% s$graph_ids, ]
    for (i in seq_len(nrow(gb))) .rdr_box(gb$x[i], gb$y[i], gb$w[i], gb$h[i], gb$style[i], .rdr_wrap_id(gb$id[i], 15), 0.6, 1, gb$badge[i])
    for (g in names(s$groups)) { q <- s$groups[[g]]; .rdr_tag(q[["hx"]], q[["hy"]], RDR_LINEAGE_GROUPS[[g]][[tx$lang]], 0.68, "#222222", "#F7F7F7", font = 2) }
    rp <- s$ref_panel; graphics::rect(rp[["x0"]], rp[["y0"]], rp[["x1"]], rp[["y1"]], col = "#FAFAFA", border = "#BBBBBB", lty = 2)
    graphics::text(0.5, rp[["title_y"]], tx$ref_panel, adj = c(0, 1), cex = 0.78, font = 2)
    rb <- s$boxes[!s$boxes$id %in% s$graph_ids, ]; gm <- src$gate_membership; txm <- src$taxonomy
    for (j in seq_len(nrow(rb))) {
      id <- rb$id[j]; .rdr_box(rb$x[j], rb$y[j], rb$w[j], rb$h[j], rb$style[j], .rdr_wrap_id(id, 18), 0.62, 2)
      if (id %in% s$ref_ids) { k <- em[em$source == id | em$target == id, ]
        other <- ifelse(k$source == id, k$target, ifelse(k$source == RDR_SPLIT_SOURCE, "HSC/MPP", k$source))
        ln <- sprintf("%s %s %s (%s)", k$source_edge_id, ifelse(k$source == id, "→", "←"), other, tx$treat[k$treatment])
        m <- gm[gm$gate_node == id, ]; if (nrow(m)) ln <- c(ln, sprintf("%s: %s %s", m$membership_id, m$subset_node, tx$membership))
      } else ln <- c(paste0(tx$tax_children, ":"), txm$child[txm$parent == id])
      for (q in seq_along(ln)) graphics::text(rb$x[j] - rb$w[j] / 2, rb$y[j] - 0.3 - (q - 1) * 0.15, ln[q], adj = c(0, 0.5), cex = 0.5, col = "#444444")
    }
    .rdr_draw_legend(s$legend_xy[1], s$legend_xy[2], tx)
    graphics::text(W - 0.4, 0.3, tx$footer, adj = c(1, 0), cex = 0.58, col = "#555555")
  }
  .rdr_save(out_dir, "rna_dag_core", W, H, draw, kind, fam)
}


# ── 그림 3: 계통별 상세 ──────────────────────────────────────────────────────
rdr_draw_lineage <- function(rv, src, group, out_dir, tx, kind, fam) {
  s <- rdr_lineage_spec(rv, src, group); K <- RDR_LIN; pn <- rv$proposed_nodes; lab <- stats::setNames(pn$label, pn$proposed_id); W <- s$W; H <- s$H; gw <- s$gw
  draw <- function() {
    graphics::text(0.35, H - 0.3, sprintf(tx$lin_title, RDR_LINEAGE_GROUPS[[group]][[tx$lang]]), adj = c(0, 1), cex = 1.15, font = 2)
    graphics::text(0.35, H - 0.6, tx$lin_sub, adj = c(0, 1), cex = 0.62, col = "#444444")
    graphics::text(0.35, H - 0.8, tx$layer_note, adj = c(0, 1), cex = 0.62, col = "#444444")
    if (!is.null(s$ref_rect)) { q <- s$ref_rect; graphics::rect(q[["x0"]], q[["y0"]], q[["x1"]], q[["y1"]], col = "#F3F3F3", border = NA)
      graphics::text((q[["x0"]] + q[["x1"]]) / 2, q[["y0"]] + 0.18, tx$ref_col, cex = 0.62, font = 2, col = "#666666") }
    if (!is.null(s$tax_rect)) { q <- s$tax_rect; graphics::rect(q[["x0"]], q[["y0"]], q[["x1"]], q[["y1"]], col = "#EEF7F7", border = NA)
      graphics::text((q[["x0"]] + q[["x1"]]) / 2, q[["y0"]] + 0.18, tx$tax_col, cex = 0.62, font = 2, col = "#2A8C8C") }
    b <- s$bracket; graphics::rect(b[["x0"]], b[["y0"]], b[["x1"]], b[["y1"]], border = RDR_EDGE_STYLE$unresolved$col, lty = 2, lwd = 1.1)
    graphics::text(b[["xc"]], b[["y1"]] + 0.12, tx$origin, cex = 0.6, font = 2, col = "#7A5A00")
    cd <- s$candidates; for (i in seq_len(nrow(cd))) .rdr_box(cd$x[i], cd$y[i], cd$w[i], cd$h[i], "split", c(cd$id[i], tx$ref_unres), c(0.62, 0.5), c(2, 1))
    for (p in s$bus) .rdr_draw_path(p, "unresolved")
    for (q in s$routes) .rdr_draw_path(q$pts, q$style)
    bx <- s$boxes
    for (i in seq_len(nrow(bx))) .rdr_box(bx$x[i], bx$y[i], bx$w[i], bx$h[i], bx$style[i], c(.rdr_wrap_id(bx$id[i], 19)[1], .rdr_wrap_text(lab[[bx$id[i]]], 25, 2)),
                                          c(0.62, 0.5, 0.5), c(2, 1, 1), bx$badge[i])
    for (q in s$routes) .rdr_tag(q$label_xy[1], q$label_xy[2], q$id, 0.48, RDR_EDGE_STYLE[[q$style]]$col)
    for (k in names(s$bus_lab)) .rdr_tag(s$bus_lab[[k]][1], s$bus_lab[[k]][2], k, 0.48, RDR_EDGE_STYLE$unresolved$col, adj = 0)
    for (t in names(s$tax_kids)) { tb <- bx[bx$id == t, ]; ln <- c(paste0(tx$tax_children, ":"), s$tax_kids[[t]])
      for (q in seq_along(ln)) graphics::text(tb$x - tb$w / 2, tb$y - tb$h / 2 - 0.14 * q, ln[q], adj = c(0, 0.5), cex = 0.5, col = "#2A8C8C") }
    .rdr_draw_legend(s$legend_xy[1], s$legend_xy[2], tx, cex = 0.55)
    x0 <- gw + 0.15; y0 <- H - 1.25; cx <- x0 + c(0, 0.62, 2.95, 4.0, 5.55); tab <- s$tab; extra <- s$extra
    graphics::text(x0, H - 1.0, sprintf(tx$tab_title, nrow(tab) + nrow(extra)), adj = c(0, 0), cex = 0.68, font = 2)
    for (k in seq_along(cx)) graphics::text(cx[k], y0, tx$tab_head[k], adj = c(0, 0.5), cex = 0.58, font = 2)
    graphics::segments(x0, y0 - 0.1, x0 + K$tab_w - 0.3, y0 - 0.1, col = "#999999")
    for (i in seq_len(nrow(tab))) {
      yy <- y0 - 0.1 - i * 0.19; st <- RDR_EDGE_STYLE[[.rdr_edge_style(tab[i, ])]]
      ev <- .rdr_split(tab$source_evidence_ids[i]); evs <- if (length(ev) > 2) paste0(paste(ev[1:2], collapse = ";"), sprintf(" +%d", length(ev) - 2)) else paste(ev, collapse = ";")
      vals <- c(tab$source_edge_id[i], sprintf("%s → %s", tab$source[i], tab$target[i]), tx$treat[[tab$treatment[i]]], sprintf("%s / %s", tab$support_class[i], tab$adult_support[i]), evs)
      for (k in seq_along(cx)) graphics::text(cx[k], yy, vals[k], adj = c(0, 0.5), cex = 0.5, col = if (k == 3) st$col else "#222222")
    }
    for (j in seq_len(nrow(extra))) { yy <- y0 - 0.1 - (nrow(tab) + j) * 0.19
      vals <- c(extra$source_edge_id[j], sprintf("%s → %s", extra$from[j], extra$to[j]), tx$treat[[extra$treatment[j]]])
      for (k in 1:3) graphics::text(cx[k], yy, vals[k], adj = c(0, 0.5), cex = 0.5, col = "#2A8C8C") }
    graphics::text(x0, 0.3, tx$footer, adj = c(0, 0), cex = 0.48, col = "#555555")
  }
  .rdr_save(out_dir, paste0("rna_dag_lineage_", group), W, H, draw, kind, fam)
}


# ── 전체 ─────────────────────────────────────────────────────────────────────
rdr_render_all <- function(rv, src, out_dir, lang = "ko") {
  tx <- RDR_TXT[[lang]]; if (is.null(tx)) stop("unknown --lang: ", lang)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  kind <- .rdr_device_kind(); fam <- .rdr_family(kind)
  if (lang == "ko" && fam == "sans") message("Hangul font not registered; Korean text may not render. Use --lang en if needed.")
  c(rdr_draw_comparison(rv, src, out_dir, tx, kind, fam), rdr_draw_core(rv, src, out_dir, tx, kind, fam),
    unlist(lapply(names(RDR_LINEAGE_GROUPS), function(g) rdr_draw_lineage(rv, src, g, out_dir, tx, kind, fam))))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); arg <- function(k, default) { i <- match(k, a); if (is.na(i)) default else a[i + 1L] }
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  source(file.path(dirname(self), "build_rna_dag_review.R"))
  src <- rdr_read_source(arg("--source", file.path("data", "priors", "rna_dag_review", "source_v1_1_0")))
  rv <- rdr_build_review(src); pd <- arg("--proposal", file.path("data", "priors", "rna_dag_review", "proposal"))
  for (k in RDR_REVIEW_TABLES)
    if (!identical(rdr_read_tsv(file.path(pd, paste0(k, ".tsv"))), rv[[k]])) stop(k, ".tsv in ", pd, " differs from a fresh build; rerun build_rna_dag_review.R")
  files <- rdr_render_all(rv, src, arg("--out", file.path("docs", "rna_dag_review", "figures")), arg("--lang", "ko"))
  cat(paste0("  ", files), sep = "\n")
}
