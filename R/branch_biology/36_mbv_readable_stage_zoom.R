# =============================================================================
# 36_mbv_readable_stage_zoom.R — 읽기 쉬운 단일 패널 progenitor 확대 그림 (보기 전용: 글자 크기 · 배치만 바꾼다)
# =============================================================================
#
# 무엇이 같은가 (audit_09 저장 표를 hash 확인 후 그대로 읽는다; h5ad 를 다시 읽지 않는다)
#   - 세포 · 원 UMAP 좌표 · 세포 원주석 색 (cell_fine_stage_colors.tsv.gz), 색표 (stage_palette.tsv)
#   - node 위치(군집 중앙값) · 최다 라벨 (node_fine_annotation.tsv), 경로 edge (fine_stage_edges.tsv), 확대 창 (fine_stage_zoom.tsv)
# 무엇이 다른가
#   - 1400 x 1400 px 한 장. 이름표 = "군집 ID  원저자 최다 fine label" (약 30 px 글자), 비율은 표(readable_labels.tsv)에만.
#   - 이름표는 서로 · node 와 겹치지 않는 자리에 두고 지시선으로 잇는다. node 와 세포 좌표는 옮기지 않는다.
#   - 아래에 짧은 색 범례와 주의 문구.
# 실행: run-r dcta R/branch_biology/36_mbv_readable_stage_zoom.R --attempt-dir <MBV attempt> [--stage-audit audit_09_fine_stage_overlay]
#        [--audit-name audit_10_readable_stage_overlay]
# =============================================================================

MBVR_SOURCES <- c("05_cell_trajectory_io.R", "10_expression_edge_graph.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                  "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "28_view_nrv.R",
                  "29_audit_nrv.R", "30_mbv_reference.R", "31_mbv_branch.R", "32_run_mbv.R", "33_mbv_graph_view.R", "34_mbv_umap_overlay.R",
                  "35_mbv_fine_stage_overlay.R", "36_mbv_readable_stage_zoom.R")
MBVR_LABEL_CEX <- 1.8                                                              # 100 dpi, 12 pt 기준 약 30 px 글자
MBVR_TITLE_FAMILY <- "mbvkorean"                                                   # 제목 · 부제에만 쓰는 한글 글꼴 이름 (이름표 글꼴은 그대로)
MBVR_TITLE_FONT_FILE <- "/System/Library/Fonts/AppleSDGothicNeo.ttc"               # 이미 설치된 macOS 시스템 글꼴 (새로 설치하지 않는다)


# ── 제목 · 부제 (선택) ─────────────────────────────────────────────────────────
#' 제목을 주면: 한글 글꼴을 쓸 수 있으면 요청 문구 그대로, 아니면 영어 대체 문구. 주지 않으면 NULL (기존 제목).
mbvr_title_spec <- function(title, subtitle, title_fb, subtitle_fb, hangul_ok) {
  if (is.null(title) || is.na(title)) return(NULL)
  if (isTRUE(hangul_ok)) return(list(title = title, subtitle = subtitle, family = MBVR_TITLE_FAMILY, fallback_used = FALSE))
  list(title = if (is.na(title_fb)) title else title_fb, subtitle = if (is.na(subtitle_fb)) subtitle else subtitle_fb, family = "", fallback_used = TRUE)
}
#' 한글 제목 글꼴 등록: quartz 비트맵 장치이고 시스템 글꼴 파일이 있으면 quartzFonts 에 이름을 붙인다(설치 없음). 성공 여부를 돌려준다.
mbvr_hangul_font <- function() {
  if (!identical(getOption("bitmapType"), "quartz") || !file.exists(MBVR_TITLE_FONT_FILE)) return(FALSE)
  tryCatch({ do.call(grDevices::quartzFonts, stats::setNames(list(grDevices::quartzFont(c("AppleSDGothicNeo-Regular", "AppleSDGothicNeo-Bold",
                                                                                          "AppleSDGothicNeo-Regular", "AppleSDGothicNeo-Bold"))), MBVR_TITLE_FAMILY)); TRUE },
           error = function(e) FALSE)
}


# ── 이름표 ────────────────────────────────────────────────────────────────────
mbvr_label_text <- function(cluster, label) paste0(cluster, "  ", ifelse(is.na(label) | !nzchar(label), "(no author label)", label))

#' 이름표 자리: node 둘레 16 방향 x 여러 거리 후보 중, 창 안이고 이미 놓은 이름표 · 모든 node 원과 겹치지 않는 첫 자리.
#'   붐비는 node 부터 놓는다. node 좌표(x, y)는 그대로 돌려준다. 자리를 못 찾으면 ok = FALSE (마지막 후보, 창 안으로).
mbvr_place_labels <- function(x, y, w, h, xlim, ylim, r_node, gap = NULL) {
  n <- length(x); if (is.null(gap)) gap <- 0.25 * max(h)
  crowd <- vapply(seq_len(n), function(i) sum(sqrt((x - x[i])^2 + (y - y[i])^2) < 3 * max(w)), 1)
  ord <- order(-crowd, seq_len(n)); ang <- seq(0, 2 * pi, length.out = 17)[-17] + pi / 16
  R <- matrix(NA_real_, n, 4, dimnames = list(NULL, c("x0", "x1", "y0", "y1"))); ok <- logical(n)
  clear <- function(r, i) {
    if (r[1] < xlim[1] || r[2] > xlim[2] || r[3] < ylim[1] || r[4] > ylim[2]) return(FALSE)
    for (j in which(!is.na(R[, 1]))) if (!(r[2] + gap <= R[j, 1] || R[j, 2] + gap <= r[1] || r[4] + gap <= R[j, 3] || R[j, 4] + gap <= r[3])) return(FALSE)
    cx <- pmin(pmax(x, r[1]), r[2]); cy <- pmin(pmax(y, r[3]), r[4])                        # 상자에서 각 node 까지 가장 가까운 점
    all(sqrt((cx - x)^2 + (cy - y)^2) >= r_node + 0.25 * gap)
  }
  for (i in ord) { last <- NULL
    for (d in c(1.5, 2.5, 3.5, 5, 7, 10) * max(h[i], r_node)) { for (a in ang) {
      cx <- x[i] + cos(a) * (d + w[i] / 2); cy <- y[i] + sin(a) * (d + h[i] / 2); r <- c(cx - w[i] / 2, cx + w[i] / 2, cy - h[i] / 2, cy + h[i] / 2); last <- r
      if (clear(r, i)) { R[i, ] <- r; ok[i] <- TRUE; break } }
      if (ok[i]) break }
    if (!ok[i]) { dx <- max(0, xlim[1] - last[1]) - max(0, last[2] - xlim[2]); dy <- max(0, ylim[1] - last[3]) - max(0, last[4] - ylim[2]); R[i, ] <- last + c(dx, dx, dy, dy) }
  }
  data.frame(x = x, y = y, lx = (R[, 1] + R[, 2]) / 2, ly = (R[, 3] + R[, 4]) / 2, x0 = R[, 1], x1 = R[, 2], y0 = R[, 3], y1 = R[, 4], ok = ok)
}

#' 지시선 끝: node 중심에서 상자 중심으로 가는 선이 상자 가장자리에 닿는 점
mbvr_leader_end <- function(px, py, x0, x1, y0, y1) {
  dx <- (x0 + x1) / 2 - px; dy <- (y0 + y1) / 2 - py
  tx <- if (dx > 0) (x0 - px) / dx else if (dx < 0) (x1 - px) / dx else -Inf
  ty <- if (dy > 0) (y0 - py) / dy else if (dy < 0) (y1 - py) / dy else -Inf
  t <- max(tx, ty, 0); c(px + t * dx, py + t * dy)
}


# ── 그림 ─────────────────────────────────────────────────────────────────────
mbvr_draw <- function(U, cs, pal, nd, ed, b, title, foot, legend_tab, ts = NULL) {
  graphics::layout(matrix(1:2, 2), heights = c(5, 1.45))
  graphics::par(mar = c(0.4, 0.4, 3.4, 0.4)); usr <- mbvu_open_panel(b); .mbvf_draw_cells(U, cs, pal, TRUE)
  span <- diff(b$xlim); r <- 0.012 * span
  if (nrow(ed)) { sg <- mbvu_edge_segments(ed, r)                                  # edge: R34 와 같은 선 (굵기 = pooled q, 같은 block 0 = 자홍색)
    for (i in which(sg$render_status == "drawn")) { lw <- .mbvg_lwd(ed$q_pooled[i]); col <- if (ed$weak[i]) "#C51B7D" else "black"
      graphics::segments(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lwd = lw + 3, col = "white")
      graphics::arrows(sg$ax0[i], sg$ay0[i], sg$ax1[i], sg$ay1[i], lty = 2, lwd = lw, col = col, length = 0.12) } }
  lab <- mbvr_label_text(nd$cluster, nd$top_author_label)
  w <- graphics::strwidth(lab, cex = MBVR_LABEL_CEX) * 1.08 + 0.02 * span; h <- graphics::strheight(lab, cex = MBVR_LABEL_CEX) * 1.9
  r_node <- graphics::strheight("M", cex = 2.4) * 0.75; inset <- 0.01 * span
  lp <- mbvr_place_labels(nd$x, nd$y, w, h, xlim = usr[1:2] + c(inset, -inset), ylim = usr[3:4] + c(inset, -inset), r_node = r_node)
  for (i in seq_len(nrow(lp))) { e <- mbvr_leader_end(lp$x[i], lp$y[i], lp$x0[i], lp$x1[i], lp$y0[i], lp$y1[i])
    graphics::segments(lp$x[i], lp$y[i], e[1], e[2], col = "grey15", lwd = 2) }
  graphics::points(nd$x, nd$y, pch = 21, cex = 2.4, bg = "white", col = "black", lwd = 1.5)
  sp <- nd$is_split %in% c("TRUE", TRUE); if (any(sp)) graphics::points(nd$x[sp], nd$y[sp], pch = 1, cex = 3.6, col = "#E69F00", lwd = 3)
  bc <- pal$color[match(nd$top_author_label, pal$author_label)]; bc[is.na(bc)] <- "grey40"
  graphics::rect(lp$x0, lp$y0, lp$x1, lp$y1, col = grDevices::adjustcolor("white", 0.92), border = bc, lwd = 3)
  graphics::text(lp$lx, lp$ly, lab, cex = MBVR_LABEL_CEX, font = 2)
  graphics::box(col = "grey60")
  if (is.null(ts)) graphics::title(main = title, cex.main = 1.35) else {          # 제목 · 부제만 바꾼 판 (나머지 그림은 같다)
    graphics::title(main = ts$title, cex.main = 1.75, line = 1.75, family = ts$family)
    graphics::mtext(ts$subtitle, side = 3, line = 0.35, cex = 1.1, family = ts$family) }
  graphics::par(mar = c(0.2, 0.4, 0.2, 0.4)); graphics::plot.new()
  graphics::legend("top", legend = legend_tab$text, pch = 15, pt.cex = 2.4, col = legend_tab$color, ncol = 3, cex = 1.2, bty = "n", x.intersp = 0.8)
  for (j in seq_along(foot)) graphics::mtext(foot[j], side = 1, line = -5.4 + 1.3 * j, cex = 0.95, font = if (j == 2L) 2 else 1)   # 아래 끝에 잘리지 않게
  invisible(lp)
}


# ── 실행기 ────────────────────────────────────────────────────────────────────
mbvr_main <- function(argv, root = NULL) {
  a <- list(attempt = NULL, stage = "audit_09_fine_stage_overlay", name = "audit_10_readable_stage_overlay",
            title = NA_character_, subtitle = NA_character_, title_fb = NA_character_, subtitle_fb = NA_character_); i <- 1L
  keys <- c("--attempt-dir" = "attempt", "--stage-audit" = "stage", "--audit-name" = "name", "--title" = "title", "--subtitle" = "subtitle",
            "--title-fallback" = "title_fb", "--subtitle-fallback" = "subtitle_fb")
  while (i <= length(argv)) {
    k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (k %in% names(keys)) a[[keys[[k]]]] <- v else { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    i <- i + 2L
  }
  ok_name <- function(x) !is.null(x) && !is.na(x) && grepl("^audit_[A-Za-z0-9_]+$", x)
  if (is.null(a$attempt) || is.na(a$attempt) || !dir.exists(a$attempt) || !ok_name(a$stage) || !ok_name(a$name)) {
    cat("--attempt-dir(있는 MBV attempt) 와 audit_ 로 시작하는 --stage-audit · --audit-name 이 필요하다.\n"); return(2L) }
  out <- file.path(a$attempt, a$name)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), a$name, ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "audit.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", MBVR_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    ld <- mdu_load_snapshot(snap, files = MBVR_SOURCES)
    get("mbvr_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

mbvr_build <- function(a, out, snap, log, argv) {
  attempt <- a$attempt; mbv_check_upstream(attempt, "view")
  d9 <- file.path(attempt, a$stage); mf <- file.path(d9, "audit_manifest.json"); if (!file.exists(mf)) stop("단계 채색 감사 manifest 가 없다: ", mf)
  man <- jsonlite::fromJSON(mf)
  chk <- function(f) { p <- file.path(d9, f); h <- man$files[[f]]; if (is.null(h) || !file.exists(p) || !identical(.nrv_sha(p), h)) stop("감사 파일 hash 가 그 manifest 와 다르다: ", f); p }
  inp <- c(cells = chk("cell_fine_stage_colors.tsv.gz"), nodes = chk("node_fine_annotation.tsv"), palette = chk("stage_palette.tsv"),
           edges = chk("fine_stage_edges.tsv"), zoom = chk("fine_stage_zoom.tsv"), stage_manifest = mf, verdict = nrv_stage_file(attempt, "verdict", "verdict_lossless.rds"))
  pr <- strsplit(man$pair, ":", fixed = TRUE)[[1]]
  ct <- utils::read.delim(gzfile(inp[["cells"]]), stringsAsFactors = FALSE,
                          colClasses = c(cell_id = "character", cluster = "character", author_label = "character", stage_group = "character", color = "character"))
  U <- cbind(ct$umap_x, ct$umap_y); cs <- ct[, c("stage_group", "color")]
  pal <- utils::read.delim(inp[["palette"]], stringsAsFactors = FALSE, colClasses = c(author_label = "character", color = "character"))
  pal_c <- pal[pal$author_label != "Other", , drop = FALSE]
  na <- utils::read.delim(inp[["nodes"]], stringsAsFactors = FALSE, colClasses = "character")
  nd <- data.frame(cluster = na$cluster, x = as.numeric(na$x_median), y = as.numeric(na$y_median), top_author_label = na$top_author_label, is_split = na$is_split,
                   stringsAsFactors = FALSE)
  ed <- utils::read.delim(inp[["edges"]], stringsAsFactors = FALSE, colClasses = "character")
  for (k in c("x0", "y0", "x1", "y1", "q_pooled")) ed[[k]] <- as.numeric(ed[[k]]); ed$weak <- ed$weak %in% "TRUE"
  zw <- utils::read.delim(inp[["zoom"]]); b <- list(xlim = c(zw$xmin, zw$xmax), ylim = c(zw$ymin, zw$ymax))
  bs <- readRDS(inp[["verdict"]])$branch_summary; n_sup <- sum(startsWith(bs$biological_status, "supported"))
  cmp <- mbvf_cmp_status(ct$author_label)
  title <- paste0(basename(attempt), "  |  ", pr[1], " vs ", pr[2], "\nprogenitor zoom (original UMAP coordinates; same cells, colours and graph as ", a$stage, ")")
  foot <- c("Colour = each cell's own original author fine label (post hoc; not an automatic annotation). Node = cluster median position (not moved).",
            paste0("Dashed arrows = anchor-conditioned route search order; differentiation direction NOT verified (supported ", n_sup, "/", nrow(bs), ")."),
            paste0("Magenta = 0 same-block cell edges.  ", sub(" \\(reference.*$", "", cmp), ".  Fractions: readable_labels.tsv"))
  leg <- rbind(data.frame(text = pal_c$author_label, color = pal_c$color, stringsAsFactors = FALSE), data.frame(text = "Other labels", color = MBVF_OTHER_COL))
  files <- character(0); put <- function(p) { files <<- c(files, p); invisible(p) }
  fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  p <- file.path(fig, paste0("fig_readable_zoom_", pr[1], "_vs_", pr[2], ".png"))
  grDevices::png(p, width = 1400, height = 1400, units = "px", res = 100, pointsize = 12)
  ts <- mbvr_title_spec(a$title, a$subtitle, a$title_fb, a$subtitle_fb, if (is.na(a$title)) FALSE else mbvr_hangul_font())
  lp <- tryCatch(mbvr_draw(U, cs, pal_c, nd, ed, b, title, foot, leg, ts), finally = grDevices::dev.off()); put(p)
  title_note <- if (is.null(ts)) "기본 제목 (배치 · 분할 이름)" else paste0("제목 \"", ts$title, "\", 부제 \"", ts$subtitle, "\"; 글꼴 ",
    if (ts$fallback_used) "기본 (한글 글꼴을 쓸 수 없어 영어 대체 문구)" else "AppleSDGothicNeo (설치된 macOS 시스템 글꼴, 제목 · 부제에만)")
  lb <- data.frame(cluster = na$cluster, label_text = mbvr_label_text(na$cluster, na$top_author_label), x_median = na$x_median, y_median = na$y_median,
                   top_fraction = na$top_fraction, top_kind = na$top_kind, label_x0 = lp$x0, label_x1 = lp$x1, label_y0 = lp$y0, label_y1 = lp$y1,
                   placed_ok = lp$ok, stringsAsFactors = FALSE)
  pth <- file.path(out, "readable_labels.tsv"); mdu_write(lb, pth); put(pth)
  rd <- c("# 읽기 쉬운 단일 패널 progenitor 확대 그림 (보기 전용)", "",
          paste0("- `", a$stage, "` 의 저장 표를 hash 확인 후 그대로 다시 그렸다. 세포 · 원 UMAP 좌표 · 세포 원주석 색 · node 위치 · 최다 라벨 · 경로 edge · 확대 창은 같다."),
          "- 바뀐 것은 배치와 글자뿐이다: 1400 x 1400 px 한 장, 이름표 = 군집 ID + 원저자 최다 fine label (약 30 px), 비율은 readable_labels.tsv 에.",
          "- 이름표는 서로 · node 와 겹치지 않는 자리에 두고 지시선으로 이었다. node 와 세포는 옮기지 않았다.",
          paste0("- ", cmp, "."), "- 화살표는 anchor 조건부 탐색 순서이며 생물학적 분화 방향은 검증되지 않았다. 색은 사후 원저자 표시이며 자동 주석 결과가 아니다.",
          paste0("- 그림 제목: ", title_note, ". 분할 · 배치 이름: ", basename(attempt), ", ", man$pair, "."))
  pth <- file.path(out, "README_ko.md"); writeLines(rd, pth); put(pth)
  pv <- data.frame(input = names(inp), path = unname(inp), sha256 = vapply(inp, .nrv_sha, "", USE.NAMES = FALSE), stringsAsFactors = FALSE)
  pth <- file.path(out, "readable_provenance.tsv"); mdu_write(pv, pth); put(pth)
  pk <- c("jsonlite", "digest"); pv_of <- function(x) tryCatch(as.character(utils::packageVersion(x)), error = function(e) "not_installed")
  env <- data.frame(item = c("R", "platform", "libPaths", paste0("package:", pk), "command", "working_dir", "png_px", "title_font", "written_at"),
                    value = c(R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"), vapply(pk, pv_of, ""), paste(c("36_mbv_readable_stage_zoom.R", argv), collapse = " "),
                              getwd(), paste0("1400x1400 @ 100 dpi, pointsize 12, label cex 1.8, bitmapType ", getOption("bitmapType")), title_note,
                              format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), stringsAsFactors = FALSE)
  pth <- file.path(out, "environment.tsv"); mdu_write(env, pth); put(pth)
  pth <- file.path(out, "source_snapshot.tsv"); mdu_write(snap, pth); put(pth)
  rel <- substring(normalizePath(files), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(audit = basename(out), stage_audit = a$stage, pair = man$pair, written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   files = as.list(setNames(vapply(files, .nrv_sha, ""), rel))), auto_unbox = TRUE, pretty = TRUE), file.path(out, "audit_manifest.json"))
  log("pair:", man$pair, "| nodes:", nrow(nd), "| edges:", nrow(ed), "| labels placed without overlap:", sum(lp$ok), "/", nrow(lp), "| cells:", nrow(ct))
}


# 스크립트로 실행할 때: bootstrap 에 필요한 파일만 먼저 읽는다. 계산은 스냅샷 사본에서 한다.
if (!interactive() && sys.nframe() == 0L) {
  .mbvr_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "27_run_nrv.R", "30_mbv_reference.R", "32_run_mbv.R", "33_mbv_graph_view.R",
               "34_mbv_umap_overlay.R", "35_mbv_fine_stage_overlay.R")) source(file.path(.mbvr_bb, .f))
  quit(status = mbvr_main(commandArgs(trailingOnly = TRUE)))
}
