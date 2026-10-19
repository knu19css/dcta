# =============================================================================
# 34_mbv_umap_overlay.R — 저장 graph(R33 audit)를 원 UMAP 위에 겹쳐 그리기 (보기 전용: 새 계산 · 평가 없음)
# =============================================================================
#
# 무엇을 그리는가 (배치마다 한 장: 왼쪽 경로 소속, 가운데 세포 원주석, 오른쪽 표)
#   - 배경: 원 h5ad 의 obsm/X_umap 전체 세포. UMAP 을 다시 맞추거나 옮기거나 돌리지 않는다. 두 패널은 같은 범위 · 같은 비율(asp = 1).
#   - 세포 순서: NRV structure 세포 ID 와 h5 obs ID 가 정확히 같은 순서여야 한다(섞임 · 빠짐 · 중복이면 멈춘다).
#   - graph: R33 audit 의 graph_nodes.csv · graph_edges.csv 를 그대로 읽는다(그 audit manifest hash 확인). 경로 · 순서 · 방향 · 숫자를 바꾸지 않는다.
#   - node 위치 = 그 군집 세포들의 원 UMAP 좌표 중앙값 (실제 세포가 아니다). 표 umap_node_xy.tsv 에 남긴다. 이름표 때문에 좌표를 옮기지 않는다:
#     node 안에는 군집 ID 만 적고, 최다 fine 라벨 · 비율 · 소속은 오른쪽 표에 적는다.
#   - 왼쪽: 회색 전체 세포 위에 경로 군집 세포(옅게)와 실제 경계 세포(연속 경로 군집 사이에 세포 edge 가 있는 세포)를 경로 색으로.
#   - 가운데: 세포마다 자신의 원저자 주석(AuthorCellType)의 broad family 색 (MBV author 단계의 고정 대응표). 군집 최다 색을 퍼뜨리지 않는다.
#   - edge: 두 node 중앙값을 잇는 점선 직선 화살표(곡선으로 꾸미지 않는다). 굵기 = pooled q (R33 과 같은 눈금). 같은 block edge 0 은 자홍색.
#     직선은 고차원 kNN graph 에 있는 군집 연결의 요약이며, UMAP 빈 공간을 지나는 실제 연속 경로가 아니다.
#   - 화살표 = anchor 조건부 경로 탐색 순서. 생물학적 방향은 검증되지 않았다(supported 수를 그림에 적는다).
# 실행: run-r dcta R/branch_biology/34_mbv_umap_overlay.R --attempt-dir <MBV attempt> --data-root <MAIN> [--graph-audit audit_03_route_graph]
#        [--audit-name audit_04_umap_overlay]
# =============================================================================

MBVU_SOURCES <- c("05_cell_trajectory_io.R", "10_expression_edge_graph.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                  "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "28_view_nrv.R",
                  "29_audit_nrv.R", "30_mbv_reference.R", "31_mbv_branch.R", "32_run_mbv.R", "33_mbv_graph_view.R", "34_mbv_umap_overlay.R")


# ── 세포 · 좌표 · 색 ─────────────────────────────────────────────────────────
#' NRV 세포 ID 와 h5 obs ID 가 정확히 같은 순서인지 (다르면 멈춘다; 순서를 고쳐 맞추지 않는다)
mbvu_align_cells <- function(cell_ids, h5_ids) {
  if (length(cell_ids) != length(h5_ids)) stop("세포 수가 다르다: NRV ", length(cell_ids), ", h5 ", length(h5_ids))
  if (anyDuplicated(cell_ids) || anyDuplicated(h5_ids)) stop("세포 ID 중복이 있다.")
  if (!setequal(cell_ids, h5_ids)) stop("한쪽에 없는 세포 ID 가 있다: ", paste(utils::head(setdiff(union(cell_ids, h5_ids), intersect(cell_ids, h5_ids)), 3), collapse = ", "))
  if (!identical(as.character(cell_ids), as.character(h5_ids))) stop("세포 ID 순서가 다르다 (h5 obs 순서와 NRV 세포 순서).")
  TRUE
}

#' node 위치: 군집 세포들의 원 UMAP 좌표 중앙값 (실제 세포가 아니다)
mbvu_node_xy <- function(U, cluster, clusters) {
  do.call(rbind, lapply(clusters, function(k) { i <- which(cluster == k); if (!length(i)) stop("세포가 없는 군집이다: ", k)
    data.frame(cluster = k, x_median = stats::median(U[i, 1]), y_median = stats::median(U[i, 2]), n_cells = length(i), stringsAsFactors = FALSE) }))
}

#' 한 배치의 edge: 저장 graph_edges 의 행을 순서 그대로 (원래 열은 바꾸지 않는다) + 그릴 좌표 + same-block 0 표시
mbvu_layout_edges <- function(ge, layout, xy) {
  e <- ge[as.character(ge$layout) == as.character(layout), , drop = FALSE]; if (!nrow(e)) return(e)
  i0 <- match(e$from, xy$cluster); i1 <- match(e$to, xy$cluster)
  if (anyNA(c(i0, i1))) stop("UMAP 위치가 없는 군집의 edge 다: ", paste(e$from[is.na(i0) | is.na(i1)], e$to[is.na(i0) | is.na(i1)], sep = "-", collapse = ", "))
  e$x0 <- xy$x_median[i0]; e$y0 <- xy$y_median[i0]; e$x1 <- xy$x_median[i1]; e$y1 <- xy$y_median[i1]
  e$weak <- !is.na(suppressWarnings(as.numeric(e$n_same_block_edges))) & as.numeric(e$n_same_block_edges) == 0
  e
}

#' 세포마다 자기 원주석의 broad family 와 색 (라벨 없음 = no_author_label, 규칙 없음 = unmapped)
mbvu_cell_colors <- function(au, bm) {
  fam <- bm$broad_family[match(au, bm$author_label)]; fam[is.na(fam)] <- "unmapped"; fam[is.na(au)] <- "no_author_label"
  data.frame(family = fam, col = unname(MBVG_FAMILY_COLS[fam]), stringsAsFactors = FALSE)
}
mbvu_color_mapping <- function(au, bm) {
  t <- table(au, useNA = "ifany"); lab <- names(t); lab[is.na(lab)] <- NA_character_
  cc <- mbvu_cell_colors(lab, bm)
  data.frame(author_label = lab, broad_family = cc$family, color = cc$col, n_cells = as.integer(t), stringsAsFactors = FALSE)
}

#' 두 패널 공통 범위 (전체 세포, 2 % 여유)
mbvu_bounds <- function(U) { p <- function(r) r + c(-1, 1) * 0.02 * diff(r); list(xlim = p(range(U[, 1])), ylim = p(range(U[, 2]))) }
mbvu_open_panel <- function(b) { graphics::plot.new(); graphics::plot.window(xlim = b$xlim, ylim = b$ylim, asp = 1); graphics::par("usr") }


#' 화살표 끝점: 양끝 여백은 기본 r, 짧은 edge 는 길이의 30 % 까지만 줄인다 (빠뜨리지 않는다). 길이 0 (두 node 위치가 같음) 만 그릴 수 없다고 적는다.
mbvu_edge_segments <- function(edges, r) {
  dx <- edges$x1 - edges$x0; dy <- edges$y1 - edges$y0; L <- sqrt(dx^2 + dy^2); ins <- pmin(r, 0.3 * L); s <- ifelse(L > 0, ins / L, NA_real_)
  data.frame(ax0 = edges$x0 + dx * s, ay0 = edges$y0 + dy * s, ax1 = edges$x1 - dx * s, ay1 = edges$y1 - dy * s, length = L,
             render_status = ifelse(L > 0, "drawn", "coincident_not_drawable"), stringsAsFactors = FALSE)
}

#' 확대 창: 짧은 edge (0 < 길이 < 3 r) 로 이어진 node 묶음마다 정사각 창 하나. 좌표는 원 UMAP 그대로 (창만 좁힌다).
#'   서로 떨어진 묶음을 한 창에 몰지 않는다. 창 폭이 max_width 이상이면 확대가 아니므로 만들지 않는다. 짧은 edge 가 없으면 빈 목록.
mbvu_zoom_windows <- function(nodes, edges, r, max_width = Inf) {
  L <- sqrt((edges$x1 - edges$x0)^2 + (edges$y1 - edges$y0)^2); s <- which(L > 0 & L < 3 * r); if (!length(s)) return(list())
  grp <- list()                                                                    # 짧은 edge 로 이어진 묶음 (edge 순서대로 합친다)
  for (i in s) { a <- edges$from[i]; b <- edges$to[i]; hit <- which(vapply(grp, function(g) a %in% g || b %in% g, logical(1)))
    if (!length(hit)) grp[[length(grp) + 1L]] <- c(a, b) else { m <- unique(c(unlist(grp[hit]), a, b)); grp <- c(grp[-hit], list(m)) } }
  out <- lapply(grp, function(cl) { i <- match(cl, nodes$cluster); x <- nodes$x_median[i]; y <- nodes$y_median[i]
    h <- max(1.0, 1.5 * max(diff(range(x)), diff(range(y))))
    list(clusters = cl, xlim = mean(range(x)) + c(-h, h), ylim = mean(range(y)) + c(-h, h)) })
  Filter(function(w) diff(w$xlim) < max_width, out)
}

#' 가까운 node 의 이름표: node 좌표는 그대로 두고, 이름표만 떨어뜨려 지시선으로 잇는다 (가까움 = 다른 node 와 거리 < 4 % 범위)
mbvu_label_positions <- function(nodes, span) {
  x <- nodes$x_median; y <- nodes$y_median; n <- length(x); D <- as.matrix(stats::dist(cbind(x, y))); diag(D) <- Inf
  close <- if (n > 1L) apply(D, 1, min) < 0.04 * span else rep(FALSE, n)
  ang <- (seq_len(n) - 1L) * (2 * pi / max(n, 1L)) + pi / 4; off <- 0.06 * span
  data.frame(nodes, label_x = ifelse(close, x + off * cos(ang), x), label_y = ifelse(close, y + off * sin(ang), y), callout = close, stringsAsFactors = FALSE)
}


# ── 그림 ─────────────────────────────────────────────────────────────────────
.mbvu_draw_graph_layer <- function(nodes, edges, view, b) {
  span <- diff(b$xlim); r <- 0.012 * span
  if (nrow(edges)) { sg <- mbvu_edge_segments(edges, r)
    for (i in which(sg$render_status == "drawn")) { e <- edges[i, ]
      col <- if (e$weak) "#C51B7D" else if (view == "route") MBVG_ROUTE_COLS[[e$segment]] else "black"; lw <- .mbvg_lwd(as.numeric(e$q_pooled))
      graphics::segments(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lty = 1, lwd = lw + 3, col = "white")          # 흰 테두리: 세포 구름 위에서도 보이게
      graphics::arrows(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lty = 2, lwd = lw, col = col, length = if (sg$length[i] < 2 * r) 0.05 else 0.1) } }
  fill <- if (view == "route") MBVG_ROUTE_COLS[nodes$membership] else rep("white", nrow(nodes))
  near <- mbvu_label_positions(nodes, span)$callout; ncex <- ifelse(near, 1.0, 2.6)              # 가까운 node 는 작은 기호 (이름표는 지시선으로 뺀다)
  graphics::points(nodes$x_median, nodes$y_median, pch = 21, cex = ncex, bg = fill, col = "black", lwd = 1.2)
  sp <- nodes$is_split %in% c("TRUE", TRUE)
  if (any(sp)) graphics::points(nodes$x_median[sp], nodes$y_median[sp], pch = 1, cex = ncex[sp] + 1.2, col = "#E69F00", lwd = 3)
  if (nrow(edges)) { short <- which(sg$render_status == "drawn" & sg$length < 3 * r)              # 짧은 edge 는 node 기호 위에 한 번 더 (가려지지 않게)
    for (i in short) { e <- edges[i, ]; col <- if (e$weak) "#C51B7D" else if (view == "route") MBVG_ROUTE_COLS[[e$segment]] else "black"
      graphics::segments(e$x0, e$y0, e$x1, e$y1, lty = 1, lwd = max(1.5, .mbvg_lwd(as.numeric(e$q_pooled)) * 0.5), col = col) } }
  tc <- if (view == "route") ifelse(nodes$membership %in% c("trunk", "arm2"), "white", "black") else rep("black", nrow(nodes))
  lp <- mbvu_label_positions(nodes, span); co <- lp$callout
  if (any(!co)) graphics::text(lp$x_median[!co], lp$y_median[!co], lp$cluster[!co], cex = 0.62, font = 2, col = tc[!co])
  if (any(co)) {                                                                  # 가까운 node: 지시선 + 흰 바탕 이름표 (node 는 제자리)
    graphics::segments(lp$x_median[co], lp$y_median[co], lp$label_x[co], lp$label_y[co], col = "grey20", lwd = 1)
    graphics::points(lp$label_x[co], lp$label_y[co], pch = 22, cex = 2.6, bg = "white", col = "grey20")
    graphics::text(lp$label_x[co], lp$label_y[co], lp$cluster[co], cex = 0.62, font = 2)
  }
}

.mbvu_footer <- function(n_sup, n_rows) {
  foot <- c("Node = median original UMAP position of the cluster's cells (not a cell). Straight dashed lines summarise cluster connections in the high-dimensional kNN graph,",
            "not continuous paths through empty UMAP space. Arrows = anchor-conditioned route search order; line width = pooled q; magenta = 0 same-block cell edges.",
            paste0("Biological differentiation direction NOT verified (supported ", n_sup, "/", n_rows, ")."))
  for (j in seq_along(foot)) graphics::mtext(foot[j], side = 1, line = 0.3 + 0.95 * j, cex = 0.6, font = if (j == 3L) 2 else 1)
}

#' 한 배치: 왼쪽 경로 소속, 가운데 세포 원주석 family, 오른쪽 표 · 범례
mbvu_draw_layout <- function(U, cells, cell_col, bnd, nodes, edges, b, title, n_sup, n_rows, fam_tab, zoom = FALSE) {
  graphics::layout(matrix(1:3, 1), widths = c(1, 1, 0.62)); us <- list()
  cp <- if (zoom) list(pch = 16, cex = 0.35, ab = 0.3) else list(pch = ".", cex = 1, ab = 0.75)  # 확대 창: 세포를 조금 크게, 경계 세포는 옅게
  # (A) 경로 소속
  graphics::par(mar = c(4.2, 1, 3.6, 1)); us$route <- mbvu_open_panel(b)
  graphics::points(U, pch = cp$pch, cex = cp$cex, col = grDevices::adjustcolor("grey78", 0.6))
  for (k in seq_len(nrow(nodes))) { i <- which(cells$cluster == nodes$cluster[k])
    graphics::points(U[i, , drop = FALSE], pch = cp$pch, cex = cp$cex, col = grDevices::adjustcolor(MBVG_ROUTE_COLS[[nodes$membership[k]]], 0.18)) }
  for (k in seq_len(nrow(edges))) { ids <- bnd[[k]]; if (length(ids))
    graphics::points(U[ids, , drop = FALSE], pch = 16, cex = if (zoom) 0.45 else 0.3, col = grDevices::adjustcolor(if (edges$weak[k]) "#C51B7D" else MBVG_ROUTE_COLS[[edges$segment[k]]], cp$ab)) }
  .mbvu_draw_graph_layer(nodes, edges, "route", b); graphics::box(col = "grey60")
  graphics::title(main = paste0(title, "\n(A) route membership: route-cluster cells (pale) + actual boundary cells"), cex.main = 0.8); .mbvu_footer(n_sup, n_rows)
  # (B) 세포 원주석
  graphics::par(mar = c(4.2, 1, 3.6, 1)); us$author <- mbvu_open_panel(b)
  graphics::points(U, pch = cp$pch, cex = cp$cex, col = grDevices::adjustcolor(cell_col, 0.7))
  .mbvu_draw_graph_layer(nodes, edges, "author", b); graphics::box(col = "grey60")
  graphics::title(main = paste0(title, "\n(B) each cell coloured by its own author annotation (broad family)"), cex.main = 0.8); .mbvu_footer(n_sup, n_rows)
  # 표 · 범례
  graphics::par(mar = c(1, 0.5, 3.6, 0.5)); graphics::plot.new(); graphics::plot.window(c(0, 1), c(0, 1))
  lab <- paste0(nodes$cluster, "  [", nodes$membership, if (any(nodes$is_split %in% c("TRUE", TRUE))) ifelse(nodes$is_split %in% c("TRUE", TRUE), ", split", "") else "", "]  ",
                ifelse(is.na(nodes$top_author_label) | nodes$top_author_label == "", "(no author label)", nodes$top_author_label), " ",
                sprintf("%.1f%%", 100 * as.numeric(nodes$top_fraction)), ifelse(nodes$top_kind == "plurality", " (plurality)", ""), "  n=", format(as.integer(nodes$cluster_n_cells), big.mark = ","))
  graphics::text(0, 0.99, "node | membership | top fine author label (fraction) | cells", adj = c(0, 1), cex = 0.85, font = 2)
  y <- 0.95 - 0.038 * (seq_along(lab) - 1L); graphics::text(0, y, lab, adj = c(0, 1), cex = 0.8)
  y0 <- min(y) - 0.06
  graphics::legend(0, y0, legend = c("trunk", "arm 1", "arm 2", "shared", "split cluster"), pch = c(21, 21, 21, 21, 1), pt.bg = c(unname(MBVG_ROUTE_COLS), NA),
                   col = c(rep("grey20", 4), "#E69F00"), pt.lwd = c(0.6, 0.6, 0.6, 0.6, 3), pt.cex = 2.2, cex = 0.8, y.intersp = 1.3, bty = "n", title = "(A) route membership")
  graphics::text(0, 0.01, paste0("Close nodes (< 4 % of the map width apart): small symbol, ID in a white box joined by a leader line (node not moved).\n",
                                 "Very short edges (< 3.6 % of the width) are also redrawn as a thin solid line on top of the node symbols."), adj = c(0, 0), cex = 0.62)
  graphics::legend(0.42, y0, legend = paste0(fam_tab$broad_family, " (", format(fam_tab$n_cells, big.mark = ","), ")"), pch = 15, pt.cex = 1.8,
                   col = unname(MBVG_FAMILY_COLS[fam_tab$broad_family]), cex = 0.72, bty = "n", title = "(B) cell author family (cells)")
  invisible(us)
}


# ── 실행기 ────────────────────────────────────────────────────────────────────
mbvu_main <- function(argv, root = NULL, reader = NULL) {
  a <- list(attempt = NULL, data_root = NULL, graph = "audit_03_route_graph", name = "audit_04_umap_overlay"); i <- 1L
  keys <- c("--attempt-dir" = "attempt", "--data-root" = "data_root", "--graph-audit" = "graph", "--audit-name" = "name")
  while (i <= length(argv)) {
    k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (k %in% names(keys)) a[[keys[[k]]]] <- v else { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    i <- i + 2L
  }
  ok_name <- function(x) !is.null(x) && !is.na(x) && grepl("^audit_[A-Za-z0-9_]+$", x)
  if (is.null(a$attempt) || is.na(a$attempt) || !dir.exists(a$attempt) || is.null(a$data_root) || is.na(a$data_root) || !ok_name(a$name) || !ok_name(a$graph)) {
    cat("--attempt-dir(있는 MBV attempt) · --data-root 와 audit_ 로 시작하는 --graph-audit · --audit-name 이 필요하다.\n"); return(2L) }
  out <- file.path(a$attempt, a$name)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), a$name, ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "audit.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", MBVU_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    ld <- mdu_load_snapshot(snap, files = MBVU_SOURCES)
    get("mbvu_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, reader = reader)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

mbvu_build <- function(a, out, snap, log, argv, reader) {
  attempt <- a$attempt
  # 1. 원 attempt 와 앞 graph 감사 확인
  mbv_check_upstream(attempt, "view")
  fz <- jsonlite::fromJSON(nrv_stage_file(attempt, "freeze", "plan_freeze.json")); m <- list(attempt = attempt, nrv_attempt = fz$args$nrv_attempt)
  gdir <- file.path(attempt, a$graph); gmf <- file.path(gdir, "audit_manifest.json"); if (!file.exists(gmf)) stop("graph 감사 manifest 가 없다: ", gmf)
  gm <- jsonlite::fromJSON(gmf)
  gchk <- function(f) { p <- file.path(gdir, f); h <- gm$files[[f]]; if (is.null(h) || !file.exists(p) || !identical(.nrv_sha(p), h)) stop("graph 감사 파일 hash 가 그 manifest 와 다르다: ", f); p }
  inp <- c(graph_nodes = gchk("graph_nodes.csv"), graph_edges = gchk("graph_edges.csv"), graph_manifest = gmf,
           cells = .mbv_nrv_path(m, "structure", "cells_lossless.rds"), edges = .mbv_nrv_path(m, "structure", "edges.rds"),
           broad_map = nrv_stage_file(attempt, "author", "broad_label_mapping.csv"), verdict = nrv_stage_file(attempt, "verdict", "verdict_lossless.rds"))
  if (is.null(reader)) reader <- nrv_h5_reader(a$data_root)
  .mbv_check_inputs(m, reader)                                                     # 원 입력 hash = NRV structure 고정값 (MBV freeze 에 고정)
  # 2. 세포 순서 · 좌표 · 원주석
  cells <- readRDS(inp[["cells"]])$cells; cells$cluster <- as.character(cells$cluster)
  mbvu_align_cells(as.character(cells$cell_id), as.character(reader$obs("cell_id")))
  U <- reader$umap(); if (!is.matrix(U) || nrow(U) != nrow(cells) || ncol(U) != 2L) stop("UMAP 행렬이 세포 수 x 2 가 아니다.")
  au <- as.character(reader$author()); if (length(au) != nrow(cells)) stop("저자 라벨 수가 세포 수와 다르다.")
  bm <- utils::read.csv(inp[["broad_map"]], stringsAsFactors = FALSE); bs <- readRDS(inp[["verdict"]])$branch_summary
  n_sup <- sum(startsWith(bs$biological_status, "supported")); n_rows <- nrow(bs)
  gn <- utils::read.csv(inp[["graph_nodes"]], stringsAsFactors = FALSE, colClasses = "character")
  ge <- utils::read.csv(inp[["graph_edges"]], stringsAsFactors = FALSE, colClasses = "character")
  xy <- mbvu_node_xy(U, cells$cluster, unique(gn$cluster))
  sz <- unique(gn[, c("cluster", "cluster_n_cells")])
  if (any(as.integer(sz$cluster_n_cells) != xy$n_cells[match(sz$cluster, xy$cluster)])) stop("graph node 의 세포 수와 세포 표의 군집 크기가 다르다.")
  cc <- mbvu_cell_colors(au, bm); mp <- mbvu_color_mapping(au, bm); b <- mbvu_bounds(U)
  fam_tab <- stats::aggregate(list(n_cells = mp$n_cells), list(broad_family = mp$broad_family), sum); fam_tab <- fam_tab[order(-fam_tab$n_cells), ]
  e <- readRDS(inp[["edges"]])
  # 3. 배치마다 한 장
  files <- character(0); put <- function(p) { files <<- c(files, p); invisible(p) }
  fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE); oe_all <- list(); zoom_rows <- list()
  for (r in sort(unique(as.integer(gn$layout)))) {
    nd <- gn[as.integer(gn$layout) == r, , drop = FALSE]; nd <- cbind(nd, xy[match(nd$cluster, xy$cluster), c("x_median", "y_median")])
    ed <- mbvu_layout_edges(ge, r, xy)
    cl <- cells$cluster; keep <- cl[e$i] %in% nd$cluster & cl[e$j] %in% nd$cluster
    sd <- if (nrow(ed) && any(keep)) nrv_set_definitions(e[keep, , drop = FALSE], cells) else NULL
    bnd <- lapply(seq_len(nrow(ed)), function(k) if (is.null(sd)) integer(0) else
      unique(c(sd$pool$sets[[paste0("pool|", ed$from[k], ">", ed$to[k])]], sd$pool$sets[[paste0("pool|", ed$to[k], ">", ed$from[k])]])))
    if (nrow(ed)) {
      sg <- mbvu_edge_segments(ed, 0.012 * diff(b$xlim)); ed$render_status <- sg$render_status; ed$umap_length <- sg$length
      ed$n_boundary_cells_drawn <- vapply(bnd, length, 1L); oe_all[[length(oe_all) + 1L]] <- ed
      log("layout", r, "| saved edges:", nrow(ed), "| drawn:", sum(sg$render_status == "drawn"), "| coincident (not drawable):", sum(sg$render_status != "drawn"))
    }
    nk <- if (nrow(ed)) ed$n_layout_keys[1] else "?"
    ttl <- paste0("[", nd$tier[1], "] ", nd$arm1[1], " vs ", nd$arm2[1], " | shown layout: ", nd$shown_structural_status[1], " (", nk, " keys)")
    p <- file.path(fig, paste0("fig_umap_overlay_", r, "_", nd$arm1[1], "_vs_", nd$arm2[1], ".png"))
    grDevices::png(p, width = 26, height = 10.5, units = "in", res = 110)
    tryCatch(mbvu_draw_layout(U, cells, cc$col, bnd, nd, ed, b, ttl, n_sup, n_rows, fam_tab), finally = grDevices::dev.off()); put(p)
    zws <- if (nrow(ed)) mbvu_zoom_windows(nd, ed, 0.012 * diff(b$xlim), max_width = 0.5 * diff(b$xlim)) else list()   # 짧은 edge 묶음마다 같은 좌표의 확대 그림
    for (z in seq_along(zws)) { zw <- zws[[z]]
      zoom_rows[[length(zoom_rows) + 1L]] <- data.frame(layout = r, arm1 = nd$arm1[1], arm2 = nd$arm2[1], zoom = z, clusters = paste(zw$clusters, collapse = ";"),
                                                        xmin = zw$xlim[1], xmax = zw$xlim[2], ymin = zw$ylim[1], ymax = zw$ylim[2], stringsAsFactors = FALSE)
      pz <- file.path(fig, paste0("fig_umap_overlay_", r, "_", nd$arm1[1], "_vs_", nd$arm2[1], "_zoom", z, ".png"))
      grDevices::png(pz, width = 26, height = 10.5, units = "in", res = 110)
      tryCatch(mbvu_draw_layout(U, cells, cc$col, bnd, nd, ed, zw, paste0(ttl, "\nZOOM on clusters ", paste(zw$clusters, collapse = ", "),
                                                                        " (same original coordinates; window in zoom_windows.tsv)"), n_sup, n_rows, fam_tab, zoom = TRUE),
               finally = grDevices::dev.off()); put(pz)
    }
  }
  p <- file.path(out, "zoom_windows.tsv")
  mdu_write(if (length(zoom_rows)) do.call(rbind, zoom_rows) else data.frame(layout = integer(0), arm1 = character(0), arm2 = character(0), zoom = integer(0), clusters = character(0),
                                                                              xmin = numeric(0), xmax = numeric(0), ymin = numeric(0), ymax = numeric(0)), p); put(p)
  # 4. 표 · 출처 · 환경
  xy$position_method <- "median of original obsm/X_umap over the cluster's cells (not a cell; coordinates not moved)"
  p <- file.path(out, "umap_node_xy.tsv"); mdu_write(xy, p); put(p)
  p <- file.path(out, "overlay_edges.tsv"); mdu_write(do.call(rbind, oe_all), p); put(p)
  p <- file.path(out, "cell_color_mapping.tsv"); mdu_write(mp, p); put(p)
  inf <- if (is.function(reader$input_files)) reader$input_files() else character(0)
  pv <- data.frame(input = c(names(inp), if (length(inf)) paste0("original:", names(inf)) else character(0)), path = c(unname(inp), unname(inf)),
                   sha256 = c(vapply(inp, .nrv_sha, "", USE.NAMES = FALSE), if (length(inf)) rep("checked against NRV structure input_sha256.tsv (frozen in MBV freeze)", length(inf)) else character(0)),
                   stringsAsFactors = FALSE)
  p <- file.path(out, "overlay_provenance.tsv"); mdu_write(pv, p); put(p)
  pk <- c("jsonlite", "digest", "Matrix", "hdf5r")                                 # hdf5r = 원 h5ad 의 X_umap · AuthorCellType 읽기
  pv_of <- function(x) tryCatch(as.character(utils::packageVersion(x)), error = function(e) "not_installed")
  env <- data.frame(item = c("R", "platform", "libPaths", paste0("package:", pk), "command", "working_dir", "written_at"),
                    value = c(R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"), vapply(pk, pv_of, ""),
                              paste(c("34_mbv_umap_overlay.R", argv), collapse = " "), getwd(), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), stringsAsFactors = FALSE)
  p <- file.path(out, "environment.tsv"); mdu_write(env, p); put(p)
  p <- file.path(out, "source_snapshot.tsv"); mdu_write(snap, p); put(p)
  rel <- substring(normalizePath(files), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(audit = basename(out), graph_audit = a$graph, written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   files = as.list(setNames(vapply(files, .nrv_sha, ""), rel))), auto_unbox = TRUE, pretty = TRUE), file.path(out, "audit_manifest.json"))
  log("layouts:", length(unique(gn$layout)), "| node positions:", nrow(xy), "| edges drawn:", sum(vapply(oe_all, nrow, 0L)), "| cells:", nrow(cells),
      "| author labels:", nrow(mp), "| supported:", n_sup, "/", n_rows)
}


# 스크립트로 실행할 때: bootstrap 에 필요한 파일만 먼저 읽는다. 계산은 스냅샷 사본에서 한다.
if (!interactive() && sys.nframe() == 0L) {
  .mbvu_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "27_run_nrv.R", "30_mbv_reference.R", "32_run_mbv.R", "33_mbv_graph_view.R")) source(file.path(.mbvu_bb, .f))
  quit(status = mbvu_main(commandArgs(trailingOnly = TRUE)))
}
