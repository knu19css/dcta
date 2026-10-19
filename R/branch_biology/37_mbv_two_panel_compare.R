# =============================================================================
# 37_mbv_two_panel_compare.R — 원저자 주석(A) 대 우리 결과(B) 2패널 비교 (보기 전용: 새 주석 · 모델 · 경로 · 통계 없음)
# =============================================================================
#
# 두 그림 (각각 2패널, 같은 세포 좌표 · 같은 창 · 같은 비율)
#   1. 초기 확대: neutrophil_early vs eosinophil_progenitor 배치, audit_09 의 확대 창 그대로.
#   2. 전체 지도: Mono · Neu · Eos 를 포함한 모든 anchor 후보와 팔마다 가장 흔한 저장 경로.
# A (왼쪽) = 원저자 세포 주석 (AuthorCellType; 사후 비교용)
#   - 초기 확대: 세포마다 자기 fine label 색 (audit_09 와 같다). 전체 지도: 세포마다 자기 라벨의 broad family 색 (audit_08 대응표).
# B (오른쪽) = 우리 결과만. 저자 라벨을 쓰지 않는다.
#   - 저장된 anchor 후보 (MBV select anchors.tsv: 역할 · 군집 · 후보로 고른 key 수), 경로 소속 (줄기 · 팔), 저장 경로 (candidate_routes.tsv).
#   - 중간 군집의 단계는 추정하지 않았으므로 "단계 미확정" 이라고 적는다. anchor 이름은 우리 참고 node panel 이름(HSC, Myeloblast, SN …)이다.
#   - 실행 중 확인: B 이름표가 저자 라벨을 섞어도 같다, 좌우 창이 같다, node 중앙값이 audit_09 와 같다 (checks.tsv).
# 화살표 = anchor 조건부 경로 탐색 순서이며 검증된 분화 방향이 아니다.
# 실행: run-r dcta R/branch_biology/37_mbv_two_panel_compare.R --attempt-dir <MBV attempt> [--stage-audit audit_09_fine_stage_overlay]
#        [--graph-audit audit_03_route_graph] [--overlay-audit audit_08_umap_overlay] [--audit-name audit_13_two_panel_compare]
# =============================================================================

MBVT_SOURCES <- c("05_cell_trajectory_io.R", "10_expression_edge_graph.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                  "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "28_view_nrv.R",
                  "29_audit_nrv.R", "30_mbv_reference.R", "31_mbv_branch.R", "32_run_mbv.R", "33_mbv_graph_view.R", "34_mbv_umap_overlay.R",
                  "35_mbv_fine_stage_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R")
MBVT_ROLE_COLS <- c(S = "#111111", neutrophil_full = "#0072B2", neutrophil_early = "#56B4E9", monocyte = "#009E73",
                    eosinophil_mature = "#CC79A7", eosinophil_progenitor = "#E69F00", shared = "#D55E00")
MBVT_TXT <- list(
  ko = list(s_cand = "S 후보", t_cand = "T 후보", trunk = "줄기", arm1 = "팔 1", arm2 = "팔 2", shared = "공유", end = "끝", off = "경로 밖",
            nd = "단계 미확정", key = "key", panel_cand = "panel 후보 (key/%d)"),
  en = list(s_cand = "S cand.", t_cand = "T cand.", trunk = "trunk", arm1 = "arm 1", arm2 = "arm 2", shared = "shared", end = "end", off = "off route",
            nd = "stage not determined", key = "keys", panel_cand = "panel candidate (keys/%d)"))
.mbvt_roles <- function() c("S", MBV_ARMS$arm)
.mbvt_node_of <- function(role) if (identical(role, "S")) "HSC" else MBV_ARMS$terminal_node[MBV_ARMS$arm == role]   # 우리 참고 node panel 이름


# ── B 의 재료: 저장된 anchor 후보 · 경로 (저자 라벨 없음) ───────────────────────
#' 역할 x 군집마다 후보로 고른 key 수 (anchors.tsv 의 candidate == TRUE 만), 전체 key 수
mbvt_anchor_counts <- function(an) {
  n <- length(unique(paste(an$holdout_scheme, an$fold))); c <- unique(an[an$candidate %in% c("TRUE", TRUE), c("holdout_scheme", "fold", "anchor", "cluster")])
  if (!nrow(c)) return(data.frame(role = character(0), cluster = character(0), n_keys_candidate = integer(0), n_keys = integer(0)))
  t <- stats::aggregate(list(n_keys_candidate = rep(1L, nrow(c))), list(role = c$anchor, cluster = as.character(c$cluster)), sum)
  t$n_keys <- as.integer(n); t[order(match(t$role, .mbvt_roles()), -t$n_keys_candidate, t$cluster, method = "radix"), , drop = FALSE]
}

#' 초기 확대 B 이름표: 경로 node 마다 우리 anchor 역할(있으면) · 경로 소속, 없으면 "단계 미확정". 경로 밖 anchor 후보도 덧붙인다.
#'   입력은 node 의 군집 · 소속 · 역할(start/intermediate/end)과 anchor 후보 수뿐이다(저자 라벨을 받지 않는다).
mbvt_b_labels_early <- function(nodes, counts, arm1, arm2, lang = "ko") {
  tx <- MBVT_TXT[[lang]]; roles <- c("S", arm1, arm2); mt <- c(trunk = tx$trunk, arm1 = tx$arm1, arm2 = tx$arm2, shared = tx$shared)
  one <- function(k, on_route, memb, nrole) {
    cr <- counts[counts$cluster == k & counts$role %in% roles, , drop = FALSE]; cr <- cr[order(match(cr$role, roles)), , drop = FALSE]
    place <- if (!on_route) tx$off else paste0(mt[[memb]], if (identical(nrole, "end")) paste0(" ", tx$end) else "")
    if (!nrow(cr)) return(c(paste0(k, "  ", place, " · ", tx$nd), ""))
    head <- paste0(k, "  ", if (any(cr$role == "S")) tx$s_cand else tx$t_cand, " · ", place)
    c(paste0(head, "\n", paste0(vapply(cr$role, .mbvt_node_of, ""), " panel ", cr$n_keys_candidate, "/", cr$n_keys, " ", tx$key, collapse = "; ")),
      paste(cr$role, collapse = ";"))
  }
  on <- lapply(seq_len(nrow(nodes)), function(i) one(nodes$cluster[i], TRUE, nodes$membership[i], nodes$node_role[i]))
  off_cl <- unique(counts$cluster[counts$role %in% roles & !counts$cluster %in% nodes$cluster])
  off <- lapply(off_cl, function(k) one(k, FALSE, NA, NA))
  data.frame(cluster = c(nodes$cluster, off_cl), text = vapply(c(on, off), `[`, "", 1), roles = vapply(c(on, off), `[`, "", 2),
             on_route = c(rep(TRUE, nrow(nodes)), rep(FALSE, length(off_cl))), stringsAsFactors = FALSE)
}

#' 전체 지도 B 이름표: anchor 후보 군집마다 우리 panel 이름과 후보로 고른 key 수 (여러 역할이면 모두)
mbvt_b_labels_whole <- function(counts, lang = "ko") {
  tx <- MBVT_TXT[[lang]]; cl <- unique(counts$cluster); n <- if (nrow(counts)) counts$n_keys[1] else 0L
  rows <- lapply(cl, function(k) { cr <- counts[counts$cluster == k, , drop = FALSE]; cr <- cr[order(match(cr$role, .mbvt_roles())), , drop = FALSE]
    data.frame(cluster = k, text = paste0(k, "  ", paste0(vapply(cr$role, .mbvt_node_of, ""), " ", cr$n_keys_candidate, collapse = " · "), "\n", sprintf(tx$panel_cand, n)),
               roles = paste(cr$role, collapse = ";"), shared = nrow(cr) > 1L, main_role = if (nrow(cr) > 1L) "shared" else cr$role[1], stringsAsFactors = FALSE) })
  do.call(rbind, rows)
}

#' 전체 지도 경로: 팔마다 가장 흔한 저장 경로(동률은 경로 문자열 순), 그 경로들의 edge (여러 팔이 쓰면 shared)
mbvt_whole_routes <- function(cr) {
  paths <- do.call(rbind, lapply(unique(cr$arm), function(a) { t <- table(cr$primary_path[cr$arm == a]); o <- order(-as.integer(t), names(t), method = "radix")[1]
    data.frame(arm = a, path = names(t)[o], n_keys = as.integer(t[o]), n_keys_total = sum(cr$arm == a), stringsAsFactors = FALSE) }))
  ed <- do.call(rbind, lapply(seq_len(nrow(paths)), function(i) { p <- strsplit(paths$path[i], ">", fixed = TRUE)[[1]]; if (length(p) < 2L) return(NULL)
    data.frame(arm = paths$arm[i], from = p[-length(p)], to = p[-1], stringsAsFactors = FALSE) }))
  key <- paste(ed$from, ed$to, sep = ">"); u <- !duplicated(key)
  edges <- data.frame(from = ed$from[u], to = ed$to[u], stringsAsFactors = FALSE)
  edges$arms <- vapply(key[u], function(k) paste(sort(unique(ed$arm[key == k]), method = "radix"), collapse = ";"), "", USE.NAMES = FALSE)
  edges$shared <- grepl(";", edges$arms, fixed = TRUE); list(paths = paths, edges = edges)
}


# ── 표 저장 (줄바꿈이 있는 이름표 글자) ─────────────────────────────────────────
#' 줄바꿈은 '\n' 두 글자로, 탭은 공백으로 (한 행 = 한 줄). 그림의 글자는 바꾸지 않는다.
mbvt_escape_text <- function(x) gsub("\t", " ", gsub("\n", "\\n", x, fixed = TRUE), fixed = TRUE)
mbvt_unescape_text <- function(x) gsub("\\n", "\n", x, fixed = TRUE)
#' 글자 열을 escape 해 쓰고, 다시 읽어 행 수 · 칸 수 · 글자 열이 같은지 확인한다 (TRUE / FALSE)
mbvt_write_tsv <- function(df, path) {
  ch <- which(vapply(df, is.character, logical(1))); for (j in ch) df[[j]] <- mbvt_escape_text(df[[j]])
  mdu_write(df, path)
  back <- utils::read.delim(path, colClasses = "character", quote = "", check.names = FALSE)
  nf <- utils::count.fields(path, sep = "\t", quote = "", blank.lines.skip = FALSE)
  same <- function(a, b) identical(ifelse(is.na(a), "<NA>", a), ifelse(is.na(b), "<NA>", b))
  nrow(back) == nrow(df) && length(nf) == nrow(df) + 1L && all(nf == ncol(df)) && all(vapply(names(df)[ch], function(k) same(back[[k]], df[[k]]), logical(1)))
}


# ── 그리기 ────────────────────────────────────────────────────────────────────
.mbvt_draw_edges <- function(ed, b, col, lwd = NULL) {                          # 흰 테두리 + 점선 화살표 (R36 과 같은 모양)
  if (!nrow(ed)) return(invisible())
  sg <- mbvu_edge_segments(ed, 0.012 * diff(b$xlim))
  for (i in which(sg$render_status == "drawn")) { lw <- if (is.null(lwd)) .mbvg_lwd(ed$q_pooled[i]) else lwd
    graphics::segments(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lwd = lw + 3, col = "white")
    graphics::arrows(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lty = 2, lwd = lw, col = col[i], length = 0.12) }
}
.mbvt_draw_labels <- function(x, y, text, border, cex, b) {                      # 겹치지 않는 자리 + 지시선 (R36 의 배치 함수)
  usr <- graphics::par("usr"); span <- diff(b$xlim); inset <- 0.01 * span
  w <- graphics::strwidth(text, cex = cex) * 1.08 + 0.02 * span; h <- graphics::strheight(text, cex = cex) * 1.5 + 0.012 * span
  lp <- mbvr_place_labels(x, y, w, h, xlim = usr[1:2] + c(inset, -inset), ylim = usr[3:4] + c(inset, -inset), r_node = graphics::strheight("M", cex = 2.4) * 0.75)
  for (i in seq_len(nrow(lp))) { e <- mbvr_leader_end(lp$x[i], lp$y[i], lp$x0[i], lp$x1[i], lp$y0[i], lp$y1[i]); graphics::segments(lp$x[i], lp$y[i], e[1], e[2], col = "grey15", lwd = 2) }
  lp$text <- text; lp$border <- border; lp
}
.mbvt_draw_boxes <- function(lp, cex) {
  graphics::rect(lp$x0, lp$y0, lp$x1, lp$y1, col = grDevices::adjustcolor("white", 0.92), border = lp$border, lwd = 3)
  graphics::text(lp$lx, lp$ly, lp$text, cex = cex, font = 2)
}
.mbvt_points <- function(U, idx, col, zoom) graphics::points(U[idx, , drop = FALSE], pch = if (zoom) 16 else ".", cex = if (zoom) 0.35 else 1, col = col)


# ── 실행기 ────────────────────────────────────────────────────────────────────
mbvt_main <- function(argv, root = NULL) {
  a <- list(attempt = NULL, stage = "audit_09_fine_stage_overlay", graph = "audit_03_route_graph", overlay = "audit_08_umap_overlay", name = "audit_13_two_panel_compare"); i <- 1L
  keys <- c("--attempt-dir" = "attempt", "--stage-audit" = "stage", "--graph-audit" = "graph", "--overlay-audit" = "overlay", "--audit-name" = "name")
  while (i <= length(argv)) {
    k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (k %in% names(keys)) a[[keys[[k]]]] <- v else { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    i <- i + 2L
  }
  ok_name <- function(x) !is.null(x) && !is.na(x) && grepl("^audit_[A-Za-z0-9_]+$", x)
  if (is.null(a$attempt) || is.na(a$attempt) || !dir.exists(a$attempt) || !all(vapply(a[c("stage", "graph", "overlay", "name")], ok_name, logical(1)))) {
    cat("--attempt-dir(있는 MBV attempt) 와 audit_ 로 시작하는 감사 이름이 필요하다.\n"); return(2L) }
  out <- file.path(a$attempt, a$name)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), a$name, ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "audit.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", MBVT_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    ld <- mdu_load_snapshot(snap, files = MBVT_SOURCES)
    get("mbvt_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

mbvt_build <- function(a, out, snap, log, argv) {
  attempt <- a$attempt; mbv_check_upstream(attempt, "view")
  achk <- function(aud, f) { d <- file.path(attempt, aud); mf <- file.path(d, "audit_manifest.json"); if (!file.exists(mf)) stop("감사 manifest 가 없다: ", mf)
    h <- jsonlite::fromJSON(mf)$files[[f]]; p <- file.path(d, f); if (is.null(h) || !file.exists(p) || !identical(.nrv_sha(p), h)) stop("감사 파일 hash 가 그 manifest 와 다르다: ", aud, "/", f); p }
  inp <- c(cells = achk(a$stage, "cell_fine_stage_colors.tsv.gz"), nodes9 = achk(a$stage, "node_fine_annotation.tsv"), palette = achk(a$stage, "stage_palette.tsv"),
           edges9 = achk(a$stage, "fine_stage_edges.tsv"), zoom = achk(a$stage, "fine_stage_zoom.tsv"), graph_nodes = achk(a$graph, "graph_nodes.csv"),
           family_map = achk(a$overlay, "cell_color_mapping.tsv"), anchors = nrv_stage_file(attempt, "select", "anchors.tsv"),
           routes = nrv_stage_file(attempt, "verdict", "candidate_routes.tsv"), verdict = nrv_stage_file(attempt, "verdict", "verdict_lossless.rds"))
  rel <- c(cells = paste0(a$stage, "/cell_fine_stage_colors.tsv.gz"), nodes9 = paste0(a$stage, "/node_fine_annotation.tsv"), palette = paste0(a$stage, "/stage_palette.tsv"),
           edges9 = paste0(a$stage, "/fine_stage_edges.tsv"), zoom = paste0(a$stage, "/fine_stage_zoom.tsv"), graph_nodes = paste0(a$graph, "/graph_nodes.csv"),
           family_map = paste0(a$overlay, "/cell_color_mapping.tsv"), anchors = "stage_select/anchors.tsv", routes = "stage_verdict/candidate_routes.tsv",
           verdict = "stage_verdict/verdict_lossless.rds")
  # 1. 읽기
  ct <- utils::read.delim(gzfile(inp[["cells"]]), stringsAsFactors = FALSE,
                          colClasses = c(cell_id = "character", cluster = "character", author_label = "character", stage_group = "character", color = "character"))
  U <- cbind(ct$umap_x, ct$umap_y)
  n9 <- utils::read.delim(inp[["nodes9"]], stringsAsFactors = FALSE, colClasses = "character")
  pal <- utils::read.delim(inp[["palette"]], stringsAsFactors = FALSE, colClasses = c(author_label = "character", color = "character")); pal_c <- pal[pal$author_label != "Other", ]
  e9 <- utils::read.delim(inp[["edges9"]], stringsAsFactors = FALSE, colClasses = "character")
  for (k in c("x0", "y0", "x1", "y1", "q_pooled")) e9[[k]] <- as.numeric(e9[[k]]); e9$weak <- e9$weak %in% "TRUE"
  zw <- utils::read.delim(inp[["zoom"]]); bz <- list(xlim = c(zw$xmin, zw$xmax), ylim = c(zw$ymin, zw$ymax)); bw <- mbvu_bounds(U)
  gn <- utils::read.csv(inp[["graph_nodes"]], stringsAsFactors = FALSE, colClasses = "character")
  fm <- utils::read.delim(inp[["family_map"]], stringsAsFactors = FALSE, colClasses = c(author_label = "character"))
  an <- utils::read.delim(inp[["anchors"]], stringsAsFactors = FALSE, colClasses = "character")
  cr <- utils::read.delim(inp[["routes"]], stringsAsFactors = FALSE, colClasses = "character")
  bs <- readRDS(inp[["verdict"]])$branch_summary; n_sup <- sum(startsWith(bs$biological_status, "supported"))
  arm1 <- n9$arm1[1]; arm2 <- n9$arm2[1]; status <- unique(gn$shown_structural_status[gn$arm1 == arm1 & gn$arm2 == arm2])[1]; n_lk <- e9$n_layout_keys[1]
  # 2. B 재료 (저자 라벨 없음)
  fam <- if (mbvr_hangul_font()) MBVR_TITLE_FAMILY else ""; lang <- if (nzchar(fam)) "ko" else "en"; tx <- MBVT_TXT[[lang]]
  counts <- mbvt_anchor_counts(an); nd9 <- data.frame(cluster = n9$cluster, membership = n9$membership, node_role = n9$node_role, stringsAsFactors = FALSE)
  bl9 <- mbvt_b_labels_early(nd9, counts, arm1, arm2, lang); blw <- mbvt_b_labels_whole(counts, lang); wr <- mbvt_whole_routes(cr)
  # 3. 위치: 군집 세포의 원 UMAP 중앙값 (초기 node 는 audit_09 와 같아야 한다)
  need <- unique(c(n9$cluster, bl9$cluster, blw$cluster, wr$edges$from, wr$edges$to)); xy <- mbvu_node_xy(U, ct$cluster, need)
  j <- match(n9$cluster, xy$cluster); med_ok <- all(abs(xy$x_median[j] - as.numeric(n9$x_median)) <= 1e-12 & abs(xy$y_median[j] - as.numeric(n9$y_median)) <= 1e-12)
  if (!med_ok) stop("다시 구한 node 중앙값이 audit_09 와 다르다.")
  P <- function(k) xy[match(k, xy$cluster), c("x_median", "y_median")]
  # 저자 라벨 독립 확인: 저자 라벨을 섞은 세포표에서 B 재료를 다시 만들어도 같다
  set.seed(20261014); ct_sh <- ct; ct_sh$author_label <- sample(ct_sh$author_label)
  b_indep <- identical(bl9, mbvt_b_labels_early(nd9, mbvt_anchor_counts(an), arm1, arm2, lang)) && identical(blw, mbvt_b_labels_whole(mbvt_anchor_counts(an), lang)) &&
    identical(mbvu_node_xy(cbind(ct_sh$umap_x, ct_sh$umap_y), ct_sh$cluster, need), xy)
  # 4. 그림
  files <- character(0); put <- function(p) { files <<- c(files, p); invisible(p) }; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  bounds <- list(); trace <- list(); tr <- function(...) trace[[length(trace) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
  ttlA <- if (lang == "ko") "A  원저자 주석" else "A  Original author labels"          # 짧은 제목 (설명은 범례 · 아래 주의 문구에)
  ttlB <- if (lang == "ko") "B  우리 후보 경로" else "B  Our candidate routes"
  foot_common <- if (lang == "ko") c("A: 원저자 라벨은 사후 비교용으로만 썼다. B: 저장된 anchor 후보(우리 marker panel 점수 선택) · 경로 소속 · 저장 경로만 썼고 단계를 추정하지 않았다.",
                                     paste0("화살표 = anchor 조건부 경로 탐색 순서 (직선은 고차원 kNN 연결의 요약이며 빈 공간을 지나는 실제 경로가 아니다). 분화 방향 미검증 (supported ",
                                            n_sup, "/", nrow(bs), "). 같은 세포 좌표 · 같은 창 · 같은 비율."))
                 else c("A: author labels used only for post hoc comparison. B: stored anchor candidates (our marker-panel scores), route membership and stored routes only; no stage inference.",
                        paste0("Arrows = anchor-conditioned search order (straight lines summarise high-dimensional kNN connections, not paths through empty space); direction NOT verified (supported ",
                               n_sup, "/", nrow(bs), "). Same cell coordinates, window and aspect."))
  open2 <- function(p, main) { grDevices::png(p, width = 2800, height = 1550, units = "px", res = 100, pointsize = 12)
    graphics::layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(5.2, 1.25)); graphics::par(oma = c(0, 0, 3.2, 0), family = fam) }
  # 4a. 초기 확대
  p1 <- file.path(fig, paste0("fig_two_panel_early_zoom_", arm1, "_vs_", arm2, ".png")); open2(p1)
  tryCatch({
    nx <- as.numeric(n9$x_median); ny <- as.numeric(n9$y_median)
    graphics::par(mar = c(0.4, 0.4, 3, 0.4)); bounds$early_A <- mbvu_open_panel(bz); .mbvf_draw_cells(U, ct[, c("stage_group", "color")], pal_c, TRUE)
    .mbvt_draw_edges(e9, bz, ifelse(e9$weak, "#C51B7D", "black"))
    labA <- mbvr_label_text(n9$cluster, n9$top_author_label); bcA <- pal_c$color[match(n9$top_author_label, pal_c$author_label)]; bcA[is.na(bcA)] <- "grey40"
    lpA <- .mbvt_draw_labels(nx, ny, labA, bcA, 1.6, bz); graphics::points(nx, ny, pch = 21, cex = 2.4, bg = "white", lwd = 1.5)
    .mbvt_draw_boxes(lpA, 1.6); graphics::box(col = "grey60"); graphics::title(main = ttlA, cex.main = 1.45)
    for (i in seq_along(labA)) tr(figure = "early_zoom", panel = "A", element = "node_label", id = n9$cluster[i], displayed_text = labA[i],
                                  source_file = rel[["nodes9"]], source_columns = "cluster;top_author_label", source_values = paste(n9$cluster[i], n9$top_author_label[i], sep = ";"))
    graphics::par(mar = c(0.4, 0.4, 3, 0.4)); bounds$early_B <- mbvu_open_panel(bz)
    .mbvt_points(U, seq_len(nrow(ct)), grDevices::adjustcolor("grey85", 0.6), TRUE)
    for (k in seq_len(nrow(n9))) .mbvt_points(U, which(ct$cluster == n9$cluster[k]), grDevices::adjustcolor(MBVG_ROUTE_COLS[[n9$membership[k]]], 0.55), TRUE)
    offc <- bl9$cluster[!bl9$on_route]; for (k in offc) .mbvt_points(U, which(ct$cluster == k), grDevices::adjustcolor("grey40", 0.55), TRUE)
    .mbvt_draw_edges(e9, bz, ifelse(e9$weak, "#C51B7D", "black"))
    pb <- P(bl9$cluster); bcB <- ifelse(grepl("(^|;)S($|;)", bl9$roles), MBVT_ROLE_COLS[["S"]], ifelse(nzchar(bl9$roles), MBVT_ROLE_COLS[[arm1]], "grey55"))
    bcB[nzchar(bl9$roles) & !grepl("(^|;)S($|;)", bl9$roles) & grepl(arm2, bl9$roles, fixed = TRUE)] <- MBVT_ROLE_COLS[[arm2]]
    lpB <- .mbvt_draw_labels(pb$x_median, pb$y_median, bl9$text, bcB, 1.45, bz)
    gl <- ifelse(grepl("(^|;)S($|;)", bl9$roles), 22, ifelse(nzchar(bl9$roles), 23, 21))
    graphics::points(pb$x_median, pb$y_median, pch = gl, cex = 2.5, bg = "white", lwd = 1.8); .mbvt_draw_boxes(lpB, 1.45)
    graphics::box(col = "grey60"); graphics::title(main = ttlB, cex.main = 1.45)
    for (i in seq_len(nrow(bl9))) tr(figure = "early_zoom", panel = "B", element = "node_label", id = bl9$cluster[i], displayed_text = bl9$text[i],
                                     source_file = paste(rel[["anchors"]], rel[["nodes9"]], sep = " + "), source_columns = "anchor;cluster;candidate + cluster;membership;node_role",
                                     source_values = paste0("roles=", bl9$roles[i], "; membership=", ifelse(bl9$on_route[i], n9$membership[match(bl9$cluster[i], n9$cluster)], "off_route")))
    graphics::par(mar = c(0, 0.4, 0, 0.4)); graphics::plot.new()
    graphics::legend("topleft", legend = c(pal_c$author_label, "Other"), pch = 15, pt.cex = 2.2, col = c(pal_c$color, MBVF_OTHER_COL), ncol = 3, cex = 1.15, bty = "n",
                     title = if (lang == "ko") "A 색 = 원저자 fine label" else "A colours = author fine label")
    gtxt <- if (lang == "ko") c("줄기 세포", "팔 1 세포", "팔 2 세포", "경로 밖 anchor 후보 세포", "S anchor 후보 (사각)", "T anchor 후보 (마름모)", "경로 군집 (원)")
            else c("trunk cells", "arm 1 cells", "arm 2 cells", "off-route anchor cells", "S anchor candidate (square)", "T anchor candidate (diamond)", "route cluster (circle)")
    graphics::legend("topright", legend = gtxt, pch = c(15, 15, 15, 15, 22, 23, 21), pt.cex = 2.2, col = c(MBVG_ROUTE_COLS[c("trunk", "arm1", "arm2")], "grey40", rep("black", 3)),
                     pt.bg = "white", ncol = 2, cex = 1.15, bty = "n", title = if (lang == "ko") "B 색 · 모양 = 우리 경로 소속 · anchor 역할" else "B colours/glyphs = our route membership / anchor role")
    graphics::mtext(c(foot_common, paste0(if (lang == "ko") "이 배치: " else "This layout: ", status, " (", n_lk, "/11 keys). ",
                                         if (lang == "ko") "magenta = 같은 block 세포 edge 0. 경로 판정 insufficient_evidence." else "magenta = 0 same-block cell edges. Route verdict insufficient_evidence.")),
                    side = 1, line = c(-3.6, -2.4, -1.2), cex = 1.0)
    graphics::mtext(if (lang == "ko") "골수 초기 후보 경로" else "Early bone marrow candidate routes", side = 3, outer = TRUE, line = 1.0, cex = 1.9, font = 2)
  }, finally = grDevices::dev.off()); put(p1)
  tr(figure = "early_zoom", panel = "A", element = "cell_colour", id = "all cells", displayed_text = "author fine label colours (8 + Other)",
     source_file = rel[["cells"]], source_columns = "author_label;stage_group;color", source_values = paste(nrow(ct), "cells"))
  tr(figure = "early_zoom", panel = "B", element = "cell_colour", id = "all cells", displayed_text = "route membership / off-route anchor / grey",
     source_file = paste(rel[["nodes9"]], rel[["anchors"]], sep = " + "), source_columns = "cluster;membership + anchor;cluster;candidate", source_values = paste(nrow(n9), "route clusters"))
  tr(figure = "early_zoom", panel = "A+B", element = "edges", id = paste(e9$from, e9$to, sep = ">", collapse = ","), displayed_text = "dashed arrows (search order)",
     source_file = rel[["edges9"]], source_columns = "from;to;segment;search_order;q_pooled;weak", source_values = paste(nrow(e9), "edges"))
  # 4b. 전체 지도 (Mono · Neu · Eos 후보)
  p2 <- file.path(fig, "fig_two_panel_whole_bm_candidates.png"); open2(p2)
  tryCatch({
    famc <- fm$color[match(ct$author_label, fm$author_label)]; famc[is.na(famc)] <- MBVG_FAMILY_COLS[["unmapped"]]
    graphics::par(mar = c(0.4, 0.4, 3, 0.4)); bounds$whole_A <- mbvu_open_panel(bw); .mbvt_points(U, seq_len(nrow(ct)), grDevices::adjustcolor(famc, 0.7), FALSE)
    wx <- P(blw$cluster); topA <- vapply(blw$cluster, function(k) { t <- sort(table(ct$author_label[ct$cluster == k]), decreasing = TRUE); names(t)[1] }, "")
    labA2 <- mbvr_label_text(blw$cluster, topA)
    lpA2 <- .mbvt_draw_labels(wx$x_median, wx$y_median, labA2, "grey30", 1.35, bw); graphics::points(wx$x_median, wx$y_median, pch = 21, cex = 2.2, bg = "white", lwd = 1.5)
    .mbvt_draw_boxes(lpA2, 1.35); graphics::box(col = "grey60"); graphics::title(main = ttlA, cex.main = 1.45)
    for (i in seq_along(labA2)) tr(figure = "whole_bm", panel = "A", element = "node_label", id = blw$cluster[i], displayed_text = labA2[i],
                                   source_file = rel[["cells"]], source_columns = "cluster;author_label", source_values = paste0("top=", topA[i]))
    graphics::par(mar = c(0.4, 0.4, 3, 0.4)); bounds$whole_B <- mbvu_open_panel(bw); .mbvt_points(U, seq_len(nrow(ct)), grDevices::adjustcolor("grey85", 0.6), FALSE)
    rc <- unique(c(wr$edges$from, wr$edges$to)); .mbvt_points(U, which(ct$cluster %in% setdiff(rc, blw$cluster)), grDevices::adjustcolor("grey55", 0.6), FALSE)
    for (k in seq_len(nrow(blw))) .mbvt_points(U, which(ct$cluster == blw$cluster[k]), grDevices::adjustcolor(MBVT_ROLE_COLS[[blw$main_role[k]]], 0.75), FALSE)
    ew <- wr$edges; pf <- P(ew$from); pt <- P(ew$to); ew$x0 <- pf$x_median; ew$y0 <- pf$y_median; ew$x1 <- pt$x_median; ew$y1 <- pt$y_median
    ecol <- vapply(seq_len(nrow(ew)), function(i) if (ew$shared[i]) "grey25" else MBVT_ROLE_COLS[[ew$arms[i]]], "")   # 여러 팔 공유 = 회색
    .mbvt_draw_edges(ew, bw, ecol, lwd = 2.5)
    pr <- P(setdiff(rc, blw$cluster)); graphics::points(pr$x_median, pr$y_median, pch = 21, cex = 1.3, bg = "white", lwd = 1)
    lpB2 <- .mbvt_draw_labels(wx$x_median, wx$y_median, blw$text, vapply(blw$main_role, function(x) MBVT_ROLE_COLS[[x]], ""), 1.35, bw)
    graphics::points(wx$x_median, wx$y_median, pch = ifelse(grepl("(^|;)S($|;)", blw$roles), 22, 23), cex = 2.4, bg = vapply(blw$main_role, function(x) MBVT_ROLE_COLS[[x]], ""), lwd = 1.5)
    .mbvt_draw_boxes(lpB2, 1.35); graphics::box(col = "grey60"); graphics::title(main = ttlB, cex.main = 1.45)
    for (i in seq_len(nrow(blw))) tr(figure = "whole_bm", panel = "B", element = "node_label", id = blw$cluster[i], displayed_text = blw$text[i],
                                     source_file = rel[["anchors"]], source_columns = "anchor;cluster;candidate", source_values = paste0("roles=", blw$roles[i]))
    for (i in seq_len(nrow(wr$paths))) tr(figure = "whole_bm", panel = "B", element = "route", id = wr$paths$arm[i], displayed_text = wr$paths$path[i],
                                          source_file = rel[["routes"]], source_columns = "arm;primary_path", source_values = paste0(wr$paths$n_keys[i], "/", wr$paths$n_keys_total[i], " keys"))
    graphics::par(mar = c(0, 0.4, 0, 0.4)); graphics::plot.new()
    fl <- unique(fm$broad_family[fm$author_label %in% ct$author_label]); fl <- fl[order(match(fl, names(MBVG_FAMILY_COLS)))]
    graphics::legend("topleft", legend = fl, pch = 15, pt.cex = 2, col = MBVG_FAMILY_COLS[fl], ncol = 4, cex = 1.0, bty = "n", title = if (lang == "ko") "A 색 = 원저자 라벨의 broad family" else "A colours = author broad family")
    rl <- c(.mbvt_roles(), "shared"); rtxt <- c(paste0(vapply(.mbvt_roles(), .mbvt_node_of, ""), " panel (", .mbvt_roles(), ")"), if (lang == "ko") "여러 역할 공유" else "shared by several roles")
    graphics::legend("topright", legend = rtxt, pch = 15, pt.cex = 2, col = MBVT_ROLE_COLS[rl], ncol = 2, cex = 1.0, bty = "n",
                     title = if (lang == "ko") "B 색 = 우리 anchor 후보 역할 (경로 선: 팔 색, 회색 = 여러 팔 공유)" else "B colours = our anchor role (route lines: arm colour, grey = shared)")
    mis <- if (lang == "ko") "성숙 SN · Eos_IV panel 후보가 A 에서는 단핵구 · B 세포 군집에 놓인다 (불일치, 숨기지 않음). 모든 경로 판정 insufficient_evidence."
           else "Mature SN / Eos_IV panel candidates sit on monocyte / B-cell clusters in A (mismatch, not hidden). All route verdicts insufficient_evidence."
    graphics::mtext(c(foot_common, mis), side = 1, line = c(-3.6, -2.4, -1.2), cex = 1.0)
    graphics::mtext(if (lang == "ko") "골수 후보 경로" else "Bone marrow candidate routes", side = 3, outer = TRUE, line = 1.0, cex = 1.9, font = 2)
  }, finally = grDevices::dev.off()); put(p2)
  tr(figure = "whole_bm", panel = "A", element = "cell_colour", id = "all cells", displayed_text = "author broad family colours", source_file = paste(rel[["cells"]], rel[["family_map"]], sep = " + "),
     source_columns = "author_label + author_label;broad_family;color", source_values = paste(nrow(ct), "cells"))
  tr(figure = "whole_bm", panel = "B", element = "cell_colour", id = "all cells", displayed_text = "our anchor role colours / route clusters / grey", source_file = paste(rel[["anchors"]], rel[["routes"]], sep = " + "),
     source_columns = "anchor;cluster;candidate + arm;primary_path", source_values = paste(nrow(blw), "anchor clusters"))
  # 5. 표 · 확인 · 출처 · 환경
  bd <- do.call(rbind, lapply(names(bounds), function(k) { u <- bounds[[k]]; s <- strsplit(k, "_")[[1]]
    data.frame(figure = paste(s[-length(s)], collapse = "_"), panel = s[length(s)], x0 = sprintf("%.17g", u[1]), x1 = sprintf("%.17g", u[2]), y0 = sprintf("%.17g", u[3]), y1 = sprintf("%.17g", u[4]),
               stringsAsFactors = FALSE) }))
  td <- do.call(rbind, trace)
  lay <- function(lp, figure, panel, cl) data.frame(figure = figure, panel = panel, cluster = cl, text = lp$text, placed_ok = lp$ok,
                                                    x0 = lp$x0, x1 = lp$x1, y0 = lp$y0, y1 = lp$y1, stringsAsFactors = FALSE)
  ll <- rbind(lay(lpA, "early_zoom", "A", n9$cluster), lay(lpB, "early_zoom", "B", bl9$cluster), lay(lpA2, "whole_bm", "A", blw$cluster), lay(lpB2, "whole_bm", "B", blw$cluster))
  ck <- data.frame(check = c("panel_bounds_equal_early", "panel_bounds_equal_whole", "b_labels_author_independent", "early_node_medians_equal_stage_audit", "b_sources_no_author_columns"),
                   pass = c(identical(bounds$early_A, bounds$early_B), identical(bounds$whole_A, bounds$whole_B), b_indep, med_ok, !any(td$panel == "B" & grepl("author", td$source_columns))),
                   detail = c("par('usr') of A == B", "par('usr') of A == B", "B labels / node positions recomputed with shuffled author labels are identical",
                              paste0("max |diff| <= 1e-12 vs ", rel[["nodes9"]]), "display_trace B rows list no author columns"), stringsAsFactors = FALSE)
  # 표: 글자 열의 줄바꿈은 '\n' 두 글자로 저장하고, 다시 읽어 행 수 · 칸 수 · 글자를 확인한다 (Codex 검증: 줄바꿈이 TSV 행을 깨뜨렸다)
  rt <- vapply(list(list(bd, "panel_bounds.tsv"), list(td, "display_trace.tsv"), list(bl9, "b_labels_early.tsv"), list(blw, "b_labels_whole.tsv"),
                    list(ll, "label_layout.tsv"), list(wr$paths, "whole_routes.tsv"), list(wr$edges, "whole_route_edges.tsv")),
               function(x) { p <- file.path(out, x[[2]]); ok <- mbvt_write_tsv(x[[1]], p); put(p); ok }, logical(1))
  ck <- rbind(ck, data.frame(check = "tables_roundtrip_ok", pass = all(rt), detail = "each saved TSV re-read: same rows, same field count per line, same text columns (newlines stored as \\n)",
                             stringsAsFactors = FALSE))
  if (!all(ck$pass)) stop("확인 실패: ", paste(ck$check[!ck$pass], collapse = ", "))
  p <- file.path(out, "checks.tsv"); if (!mbvt_write_tsv(ck, p)) stop("checks.tsv 다시 읽기 확인 실패"); put(p)
  rd <- c("# 원저자 주석(A) 대 우리 결과(B) 2패널 비교 (보기 전용)", "",
          "- 그림 1 (초기 확대): `figures/fig_two_panel_early_zoom_*.png`. audit_09 의 확대 창 · 세포 색 · node · edge 를 그대로 썼다.",
          "- 그림 2 (전체 지도): `figures/fig_two_panel_whole_bm_candidates.png`. Mono · Neu · Eos 를 포함한 모든 anchor 후보와 팔마다 가장 흔한 저장 경로.",
          "- A = 원저자 세포 주석(사후 비교). B = 우리 결과만: 저장 anchor 후보(우리 marker panel 점수 선택, 후보로 고른 key 수) · 경로 소속 · 저장 경로. 중간 군집 단계는 추정하지 않았다(단계 미확정).",
          "- B 의 이름 HSC · Myeloblast · SN · Mono_classical · Eos_IV · EoP 는 우리 참고 DAG 의 node panel 이름이다(원저자 라벨이 아니다).",
          "- 확인(checks.tsv): 좌우 창이 같다, B 이름표 · node 위치는 저자 라벨을 섞어도 같다, 초기 node 중앙값 = audit_09, B 출처에 저자 라벨 열이 없다.",
          "- 표시 하나하나의 원본 파일 · 열 · 값: display_trace.tsv. 이름표 배치(겹침 없이 놓였는지): label_layout.tsv.",
          "- 표의 글자 열에서 그림의 줄바꿈은 '\\n' 두 글자로 저장했다(한 행 = 한 줄; 다시 읽어 확인, checks.tsv 의 tables_roundtrip_ok).",
          "- 화살표는 탐색 순서이며 분화 방향은 검증되지 않았다. 직선은 고차원 kNN graph 의 군집 연결 요약이며 UMAP 빈 공간을 지나는 실제 경로가 아니다.",
          paste0("- 성숙 SN · Eos_IV panel 후보의 불일치와 경로 판정(insufficient_evidence, supported ", n_sup, "/", nrow(bs), ")을 그림에 그대로 적었다."))
  p <- file.path(out, "README_ko.md"); writeLines(rd, p); put(p)
  pv <- data.frame(input = unname(rel), path = unname(inp[names(rel)]), sha256 = vapply(inp[names(rel)], .nrv_sha, "", USE.NAMES = FALSE), stringsAsFactors = FALSE)
  p <- file.path(out, "two_panel_provenance.tsv"); mdu_write(pv, p); put(p)
  pk <- c("jsonlite", "digest", "stats", "graphics", "grDevices", "utils"); pv_of <- function(x) tryCatch(as.character(utils::packageVersion(x)), error = function(e) "not_installed")
  env <- data.frame(item = c("R_executable", "R", "platform", "libPaths", paste0("package:", pk), "bitmapType", "font_family", "command", "working_dir", "png_px", "written_at"),
                    value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"), vapply(pk, pv_of, ""), getOption("bitmapType"),
                              if (nzchar(fam)) "AppleSDGothicNeo (installed macOS system font via quartzFonts)" else "default (Korean font unavailable; English text)",
                              paste(c("37_mbv_two_panel_compare.R", argv), collapse = " "), getwd(), "2800x1550 @ 100 dpi, pointsize 12", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                    stringsAsFactors = FALSE)
  p <- file.path(out, "environment.tsv"); mdu_write(env, p); put(p)
  p <- file.path(out, "source_snapshot.tsv"); mdu_write(snap, p); put(p)
  relf <- substring(normalizePath(files), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(audit = basename(out), stage_audit = a$stage, graph_audit = a$graph, overlay_audit = a$overlay, written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   files = as.list(setNames(vapply(files, .nrv_sha, ""), relf))), auto_unbox = TRUE, pretty = TRUE), file.path(out, "audit_manifest.json"))
  log("early labels A/B:", nrow(lpA), "/", nrow(lpB), "placed ok", sum(lpA$ok), "/", sum(lpB$ok), "| whole labels A/B placed ok", sum(lpA2$ok), "/", sum(lpB2$ok),
      "| anchor clusters:", nrow(blw), "| arms:", nrow(wr$paths), "| checks:", paste(ck$check, ck$pass, collapse = ", "), "| lang:", lang)
}


# 스크립트로 실행할 때: bootstrap 에 필요한 파일만 먼저 읽는다. 계산은 스냅샷 사본에서 한다.
if (!interactive() && sys.nframe() == 0L) {
  .mbvt_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "27_run_nrv.R", "30_mbv_reference.R", "32_run_mbv.R", "33_mbv_graph_view.R",
               "34_mbv_umap_overlay.R", "35_mbv_fine_stage_overlay.R", "36_mbv_readable_stage_zoom.R")) source(file.path(.mbvt_bb, .f))
  quit(status = mbvt_main(commandArgs(trailingOnly = TRUE)))
}
