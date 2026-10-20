# =============================================================================
# 76_gse253355_manual_panel.R — GSE253355: 최초 수기 마커 패널 (원본 XLS, 35 유형) 대 현재 slide2 + T (30 유형) — 같은 R UMAP · 군집 · 공통 대조 유전자로 점수 비교
# =============================================================================
#
# 계약: Codex 인계서 codex_manual_panel_instruction.md (2026-10-20, 사용자 "마커패널을 완전 처음에 손으로 만든거 써보면 어떤지"). Ruling: claude_manual_panel_progress.md
#   - 원본: main data/priors/final_unique_gene.xls (sha256 고정). 형식 추출은 scripts/r76_extract_manual_panel_xls.py (bonemarrow_scrna python · xlrd) 가 먼저 하고,
#     이 실행기는 그 추출 CSV · 출처 JSON 을 원본 hash 와 대조한 뒤 원본 XLS · 추출 CSV · 출처 · 이름 교정표 사본을 결과에 보존한다.
#   - 이름: main data/priors/marker_symbol_corrections.tsv 의 status alias · typo 행만 적용한다 (원문 기호 정확 일치). unresolved · 표에 없는 이름은 원문 그대로 두고
#     측정 여부를 보고한다. 그 뒤 유형 · 유전자 · 부호 중복만 지운다. 같은 유형의 양 · 음 충돌은 원본대로 양쪽에 두고 따로 보고한다 (R75 의 충돌 보류를 쓰지 않는다).
#     Mono · CD14 Mono · CD16 Mono 는 그대로 둔다 (비독립 계층 이름). 머리글 원문을 유형 이름으로 쓰고, 그림 표시명 (공백 제거 · ï → i) 과 색 대응은 표로 밝힌다.
#   - 두 arm: manual_original (35 유형) · current_slide2 (R75 slide2 + T, 30 유형). 두 arm 의 active 유전자 합집합을 배경에서 빼고 공통 대조 유전자를 한 번 만든다
#     (R75 와 대조군이 다르다: current_slide2 도 이 실험에서만 다시 점수화). 점수 · 순위는 R75 함수 그대로 (합의 판정 없음).
#   - 자료 · 26 군집 · 82,742 세포 · gse253355_umap_03 좌표 · 18,892 세포 부분집합은 R75 와 같다. 저자 라벨은 점수 고정 (calc_manifest.json) 뒤 사후 대조 · 표시만.
#   - 그림 2 장 (전체 · 계통 부분집합): 저자 / 현재판 / 최초 수기판 3 열. 잠정 최고 점수 후보일 뿐 분화 입증 · 정답이 아니다. 결과 · 공개 폴더가 비어 있지 않으면 쓰지 않는다.
#
# 실행 (W 에서): 1) 형식 추출 (python, xlrd)  2) run-r dcta R/branch_biology/76_gse253355_manual_panel.R --out-dir results/neutrophil_atlas_rna_2026-10-19/gse253355_manual_panel_01
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

GSP_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                 "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                 "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R", "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R",
                 "73_gse253355_umap.R", "74_gse253355_marker_umap.R", "75_gse253355_hsc_annotation.R", "76_gse253355_manual_panel.R")
GSP_XLS <- "/Users/edemhw/Honglab/CS/04_Claude/14_dcta/data/priors/final_unique_gene.xls"                     # main repo, 읽기 전용
GSP_CORRECTIONS <- "/Users/edemhw/Honglab/CS/04_Claude/14_dcta/data/priors/marker_symbol_corrections.tsv"      # main repo, 읽기 전용
GSP_MANUAL_INPUT <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_manual_panel_input_01"                  # 형식 추출 결과 (python xlrd)
GSP_BASELINE <- "results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_manual_panel_baseline.json"
GSP_R75 <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_hsc_annotation_02"                                # R75 결과 (대조군 차이 사후 대조만)
GSP_PUBLIC <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_manual_panel_01"
GSP_PINNED <- list(xls = "3ae23f8eb1b5a19f511ac544a478feb8c992957c24f12047cb97bbe9a6d33f34",
                   markers = "55d8f2056ffcdfb53f6dad9cf84af68d53e328ba8b82405d95b838524447f4c0", umap = "1f91bc928f06b736980bd56e18a595346ff1ddc430a96a021e7537149ab3001f",
                   lossless = "c423ecfe33be76951acfc3e71ea0173e949c908d289839a19828e48b3d6bb215", object = "83fd17576c48e710e5a03a60be485beedbc088f82c6a97fa79d7b80d5367186a",
                   calc_manifest = "5842fd04cf6e0fa943538ac671c11bbc778e4439db71662f38e4fd46871c6955")
GSP_SHEETS <- c(final_unique_gene = 1L, `Negative markers` = -1L)
GSP_ARMS <- c("current_slide2", "manual_original")
GSP_ARM_LABEL <- c(current_slide2 = "현재판 slide2 + T", manual_original = "최초 수기판 (원본 XLS)")
GSP_CORE <- list(current_slide2 = c("HSC", "MPP", "CMP", "GMP", "Neutrophil", "Classical_Monocyte", "Nonclassical_Monocyte"),
                 manual_original = c("HSC", "MPP", "CMP", "GMP", "Neutrophils", "Mono", "CD14 Mono", "CD16 Mono"))
GSP_COLOUR_MAP <- c(Macrophages = "Macrophage", `CD14 Mono` = "Classical Monocyte", `CD16 Mono` = "Nonclassical Monocyte", Neutrophils = "Neutrophil", DC = "cDC",
                    `Plasma B` = "Plasma cell", Basophils = "Basophil", NK = "NK cell", NK_T = "NKT cell", Mast = "Mast cell", `CD4 T` = "CD4 T cell",
                    `Cytotoxic CD8` = "CD8 T cell", Platelates = "Platelet")                                   # 그림 색 · 비교 열 전용 표시명 대응 (표의 유형 이름은 원문)
GSP_VERIFY <- list(c(arm = "manual_original", type = "Neutrophils"), c(arm = "manual_original", type = "MPP"), c(arm = "current_slide2", type = "HSC"))


# ── 1. 원본 패널 (추출 CSV → 이름 교정 · 중복 · 충돌) ──────────────────────────
#' 원본 추출 행 + 이름 교정표 → 출처 표 (모든 칸) · 사용 행 (유형 · 유전자 · 부호 중복만 제거) · 양 · 음 충돌 표 · 표시명 대응
gsp_manual_panel <- function(ext, cr) {
  need <- c("sheet", "sign", "excel_row", "excel_col", "col_letter", "header_raw", "symbol_raw", "symbol"); miss <- setdiff(need, names(ext))
  if (length(miss)) stop("추출 CSV 에 열이 없다: ", paste(miss, collapse = ", "))
  sg <- as.integer(ext$sign); okc <- cr[cr$status %in% c("alias", "typo"), , drop = FALSE]; hit <- match(ext$symbol, cr$source_symbol); okh <- match(ext$symbol, okc$source_symbol)
  gene <- ifelse(!is.na(okh), okc$target_symbol[okh], ext$symbol)
  corr <- ifelse(!is.na(okh), paste0("corrected_", okc$status[okh]), ifelse(!is.na(hit), paste0("kept_", cr$status[hit]), "kept_not_in_table"))
  disp <- trimws(gsub("ï", "i", ext$header_raw, fixed = TRUE))
  dup <- duplicated(paste(ext$header_raw, gene, sg, sep = "\r")); cf <- stats::ave(sg, ext$header_raw, gene, FUN = function(s) length(unique(s))) > 1
  prov <- data.frame(ext[, need], gene = gene, correction = corr, correction_category = ifelse(!is.na(hit), cr$category[hit], ""), duplicate_dropped = dup, sign_conflict_kept = cf,
                     display_label = disp, stringsAsFactors = FALSE)
  used <- data.frame(node_id = ext$header_raw[!dup], display_label = disp[!dup], gene = gene[!dup], sign = sg[!dup], sheet = ext$sheet[!dup],
                     excel_cell = paste0(ext$col_letter, ext$excel_row)[!dup], symbol_raw = ext$symbol_raw[!dup], correction = corr[!dup], stringsAsFactors = FALSE)
  cfl <- unique(prov[cf, c("header_raw", "gene")]); conflicts <- do.call(rbind, lapply(seq_len(nrow(cfl)), function(i) { z <- prov[prov$header_raw == cfl$header_raw[i] & prov$gene == cfl$gene[i], ]
    data.frame(header_raw = cfl$header_raw[i], gene = cfl$gene[i], positive_cells = paste(paste0(z$col_letter, z$excel_row, ":", z$symbol_raw)[z$sign == "1"], collapse = "|"),
               negative_cells = paste(paste0(z$col_letter, z$excel_row, ":", z$symbol_raw)[z$sign == "-1"], collapse = "|"), handling = "kept in both positive and negative panels (original)",
               stringsAsFactors = FALSE) }))
  hdr <- unique(ext$header_raw); dm <- data.frame(header_raw = hdr, display_label = trimws(gsub("ï", "i", hdr, fixed = TRUE)), stringsAsFactors = FALSE)
  dm$display_reason <- ifelse(dm$display_label == hdr, "same as header", "figure display only: surrounding spaces removed and/or i-diaeresis written as i (missing in the figure font)")
  dm$colour_key <- ifelse(dm$display_label %in% names(GSP_COLOUR_MAP), GSP_COLOUR_MAP[dm$display_label], dm$display_label)
  dm$hierarchy_note <- ifelse(dm$display_label %in% c("Mono", "CD14 Mono", "CD16 Mono"), "non-independent hierarchy names (Mono contains CD14 / CD16 Mono); kept as in the original", "")
  list(prov = prov, used = used, conflicts = if (is.null(conflicts)) prov[0, c("header_raw", "gene")] else conflicts, display_map = dm)
}


# ── 2. 그림: 저자 / 현재판 / 최초 수기판 3 열 ───────────────────────────────────
#' tops = arm 별 data.frame(cluster, key, label, status). kcol = 색 키별 색. mask = NULL 이면 전체 세포 (번호 = 모든 군집), 아니면 부분집합 (번호 = label_min_frac 이상 군집에 "번호 · 이름")
gsp_draw_three <- function(path, U, author, acols, cluster, tops, kcol, ntypes, mask = NULL, subset_labels = NULL, label_min_frac = 0, sub = "", ko = TRUE, fam = "") {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; cluster <- as.character(cluster); full <- is.null(mask); keep <- if (full) rep(TRUE, nrow(U)) else as.logical(mask)
  Us <- U[keep, , drop = FALSE]; au <- as.character(author[keep]); cl <- cluster[keep]; ns <- nrow(Us); fmt <- .nso_fmt_n; pad <- function(r) r + c(-1, 1) * 0.03 * diff(r)
  b <- if (full) mbvu_bounds(U) else list(xlim = pad(range(Us[, 1])), ylim = pad(range(Us[, 2])))
  cnt <- table(factor(cl, levels = tops[[1]]$cluster)); shown <- names(cnt)[cnt > 0]; numbered <- if (full) shown else shown[as.integer(cnt[shown]) >= label_min_frac * ns]
  undet <- if (ko) "미판정" else "undetermined"; pcex <- if (full) 0.26 else 0.4; out_cols <- list()
  grDevices::png(path, width = 5400, height = 2550, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:6, 2, byrow = TRUE), heights = c(3.15, 1.25)); graphics::par(oma = c(4.6, 0, 7.2, 0), family = fam)
    # 지도 3 개
    graphics::par(mar = c(0.4, 0.6, 5.4, 0.6)); mbvu_open_panel(b); o <- .gsu_order(au)
    graphics::points(Us[o, 1], Us[o, 2], pch = 16, cex = pcex, col = acols[au[o]]); graphics::box(col = "grey60")
    if (!full) { rp <- atlas_cl_rep_cells(Us, au); lp <- .acl_label_layout(rp, 1.6, b, rp$n_cells); .acl_draw_labels(lp, lp$text, 1.6); for (k in lp$text) note("author", "label", k) }
    t1 <- if (ko) "저자 주석 (사후 대조, 정답 아님)" else "Author annotation (post hoc, not truth)"; graphics::title(main = t1, cex.main = 2.1, line = 2.7); note("author", "title", t1)
    s1 <- sprintf(if (ko) "라벨 %d 개" else "%d labels", length(unique(au))); graphics::mtext(s1, side = 3, line = 0.7, cex = 1.35, col = "grey20"); note("author", "sub", s1)
    for (arm in names(tops)) { tp <- tops[[arm]]; key <- stats::setNames(tp$key, tp$cluster); pcol <- unname(kcol[key[cl]])
      graphics::par(mar = c(0.4, 0.6, 5.4, 0.6)); mbvu_open_panel(b); o2 <- .gsu_order(cl)
      graphics::points(Us[o2, 1], Us[o2, 2], pch = 16, cex = pcex, col = pcol[o2]); graphics::box(col = "grey60")
      rc <- atlas_cl_rep_cells(Us, cl); rc <- rc[rc$cluster %in% numbered, , drop = FALSE]
      if (!full) rc$cluster <- vapply(rc$cluster, function(k) { i <- match(k, tp$cluster); sprintf("%s · %s", k, if (tp$status[i] == "undetermined") undet else tp$label[i]) }, "")
      cx <- if (full) 1.55 else 1.7; lc <- .acl_label_layout(rc, cx, b, rc$n_cells); .acl_draw_labels(lc, lc$text, cx)
      for (k in lc$text) note(arm, "label", k); if (any(!lc$free)) note(arm, "labels_overlapping", paste(lc$text[!lc$free], collapse = " | "))
      t2 <- sprintf(if (ko) "%s (%d 유형)" else "%s (%d types)", if (ko) GSP_ARM_LABEL[[arm]] else arm, as.integer(ntypes[[arm]]))
      s2 <- sprintf(if (ko) "군집 1 위 후보 (잠정) · 미판정 %d · 동률 %d" else "cluster top candidate (provisional) · undetermined %d · tied %d",
                    sum(tp$status[match(shown, tp$cluster)] == "undetermined"), sum(tp$status[match(shown, tp$cluster)] == "tied"))
      graphics::title(main = t2, cex.main = 2.1, line = 2.7); graphics::mtext(s2, side = 3, line = 0.7, cex = 1.35, col = "grey20"); note(arm, "title", t2); note(arm, "sub", s2)
      out_cols[[arm]] <- data.frame(cluster = shown, key = unname(key[shown]), colour = vapply(shown, function(k) pcol[cl == k][1], "", USE.NAMES = FALSE),
                                    n_colours = vapply(shown, function(k) length(unique(pcol[cl == k])), 0L, USE.NAMES = FALSE), stringsAsFactors = FALSE) }
    # 범례 3 개
    graphics::par(mar = c(0.2, 1, 0.6, 1)); graphics::plot.new(); labs <- if (full) names(acols)[names(acols) %in% au] else subset_labels; tb <- table(au)
    lg <- sprintf("%s (%s)", labs, fmt(as.integer(tb[labs])))
    graphics::legend("top", legend = lg, pch = 16, col = acols[labs], pt.cex = 2.6, cex = if (full) 1.55 else 1.8, bty = "n", ncol = if (full) 3L else 2L, y.intersp = 1.0, x.intersp = 0.8,
                     title = if (ko) "저자 라벨 (세포 수)" else "author label (cells)", title.adj = 0)
    for (k in lg) note("author_legend", "item", k)
    for (arm in names(tops)) { tp <- tops[[arm]]; kk <- unique(tp$key[match(shown, tp$cluster)]); kk <- c(setdiff(kk, "undetermined"), intersect("undetermined", kk))
      it <- vapply(kk, function(k) { cs <- shown[tp$key[match(shown, tp$cluster)] == k]; nm <- if (k == "undetermined") undet else paste(unique(tp$label[tp$cluster %in% cs]), collapse = " / ")
        ids <- if (full) cs else intersect(cs, numbered); idt <- if (!length(ids)) (if (ko) "번호 없음" else "unnumbered") else if (length(ids) <= 5L) paste(ids, collapse = ",") else sprintf(if (ko) "%d 군집" else "%d clusters", length(ids))
        sprintf(if (ko) "%s  [%s] %s 세포" else "%s  [%s] %s cells", nm, idt, fmt(sum(cnt[cs]))) }, "")
      graphics::par(mar = c(0.2, 1, 0.6, 1)); graphics::plot.new()
      graphics::legend("top", legend = it, pch = 15, col = kcol[kk], pt.cex = 3.2, cex = 1.65, bty = "n", ncol = if (length(it) > 8L) 2L else 1L, y.intersp = 1.05, x.intersp = 0.8,
                       title = if (ko) "1 위 후보 [군집] 세포" else "top candidate [clusters] cells", title.adj = 0)
      for (k in it) note(paste0(arm, "_legend"), "item", k) }
    main <- if (full) sprintf(if (ko) "GSE253355 전체 %s 세포 — 저자 주석 / 현재 slide2 (30 유형) / 최초 수기 패널 (35 유형), 같은 R UMAP · 같은 26 군집" else "GSE253355, all %s cells - author / current slide2 (30 types) / original manual panel (35 types), same R UMAP and 26 clusters", fmt(ns))
            else sprintf(if (ko) "HSC 부터 Neutrophil · Monocyte 까지 저자 라벨 %d 개 (%s 세포) — 저자 / 현재판 / 최초 수기판, 같은 좌표" else "Author labels from HSC to Neutrophil / Monocyte (%d labels, %s cells) - author / current / original manual, same coordinates", length(subset_labels), fmt(ns))
    graphics::mtext(main, side = 3, outer = TRUE, line = 4.3, cex = 2.4, font = 2); graphics::mtext(sub, side = 3, outer = TRUE, line = 1.9, cex = 1.35, col = "grey20")
    note("figure", "title", main); note("figure", "sub", sub)
    foot <- if (ko) c("잠정 최고 점수 후보일 뿐 분화 입증 · 정답이 아니다. 유형 수 (35 vs 30) 와 세분화가 달라 단순 정답률이나 마커 유전자만의 효과로 해석하지 않는다.",
                      "색 = 표시명 대응 (display_name_map.csv) 기준 같은 이름 계열은 같은 색. 표에는 원본 머리글 · 기호를 그대로 둔다. 점수 = 군집 평균 양성 - max(0, 음성), 0 이하 = 미판정.")
            else c("Provisional top-score candidates only; not proof of differentiation or truth. Type numbers (35 vs 30) and granularity differ, so this is not an accuracy or gene-only comparison.",
                   "Colour = same colour for corresponding display names (display_name_map.csv); tables keep the original headers and symbols. Score = cluster mean positive - max(0, negative); <= 0 = undetermined.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.8 + (k - 1) * 1.6, cex = 1.25, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  list(text = nf$table(), coords = Us, numbered = numbered, colours = out_cols)
}


# ── 3. 실행기 ────────────────────────────────────────────────────────────────
gsp_main <- function(argv, root = NULL, pins = GSP_PINNED) {
  keys <- c("--xls" = "xls", "--corrections" = "corrections", "--manual-input" = "manual_input", "--markers" = "markers", "--umap" = "umap", "--scores-dir" = "scores_dir",
            "--object" = "object", "--baseline" = "baseline", "--r75" = "r75", "--out-dir" = "out", "--public-dir" = "public", "--assay" = "assay"); a <- list(); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (is.null(a$out) || is.na(a$out) || !nzchar(a$out)) { cat("--out-dir 가 필요하다.\n"); return(2L) }
  if (is.null(root)) root <- mdu_find_root()
  for (k in c("xls", "corrections")) if (is.null(a[[k]])) a[[k]] <- if (k == "xls") GSP_XLS else GSP_CORRECTIONS
  if (is.null(a$markers)) a$markers <- "/Users/edemhw/Honglab/CS/04_Claude/14_dcta/docs/cytoscape/rna_marker_review_2026-10-18/sign_revision_2026-10-18/node_markers_sign_review_2026-10-18.csv"
  dflt <- list(manual_input = GSP_MANUAL_INPUT, umap = "results/neutrophil_atlas_rna_2026-10-19/gse253355_umap_03/umap_coordinates.rds", scores_dir = "results/neutrophil_atlas_rna_2026-10-19/gse253355_02",
               object = "data/external/gse253355/S2_steady_GSE253355.local_copy.rds", baseline = GSP_BASELINE, r75 = GSP_R75, public = GSP_PUBLIC)
  for (k in names(dflt)) if (is.null(a[[k]])) a[[k]] <- file.path(root, dflt[[k]])
  if (is.null(a$assay)) a$assay <- "RNA"
  for (d in c(a$out, a$public)) if (dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", d, "\n"); return(6L) }
  out <- a$out; dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", GSP_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(GSP_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = GSP_SOURCES)
    get("gsp_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins, root = root)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

gsp_build <- function(a, out, snap, log, argv, pins, root) {
  ck <- .mbvq_checks(); t0 <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"); rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  files <- character(0); crt <- logical(0); fmt <- .nso_fmt_n
  wc <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  rdc <- function(p) utils::read.csv(p, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8")
  sd <- a$scores_dir; f_ls <- file.path(sd, "lossless.rds"); f_cm <- file.path(sd, "calc_manifest.json"); f_cc <- file.path(sd, "cluster_candidates.csv")
  f_ext <- file.path(a$manual_input, "manual_panel_extracted.csv"); f_prov <- file.path(a$manual_input, "manual_panel_extract_provenance.json")
  f_r75a <- file.path(a$r75, "cluster_provisional_annotation.csv"); f_r75c <- file.path(a$r75, "control_genes.csv")
  # ── (1) 입력 고정 · 기준선 ──
  inp <- c(xls = a$xls, corrections = a$corrections, extracted = f_ext, extract_provenance = f_prov, markers = a$markers, umap = a$umap, lossless = f_ls, calc_manifest = f_cm,
           cluster_candidates = f_cc, object = a$object, baseline = a$baseline, r75_annotation = f_r75a, r75_controls = f_r75c)
  for (p in inp) if (!file.exists(p)) stop("입력 파일이 없다: ", p)
  if (grepl("\\.part$", a$object)) stop("다운로드 중인 .part 입력은 쓰지 않는다: ", a$object)
  for (k in names(GSP_PINNED)) if (is.null(pins[[k]])) stop("hash 고정값이 없다: ", k)
  h0 <- vapply(inp, .nrv_sha, ""); for (k in names(GSP_PINNED)) if (!identical(h0[[k]], pins[[k]])) stop(k, " hash 가 고정값과 다르다: ", h0[[k]])
  if (!identical(h0[["cluster_candidates"]], jsonlite::fromJSON(f_cm)$files[["cluster_candidates.csv"]])) stop("cluster_candidates.csv hash 가 calc_manifest.json 기록과 다르다.")
  pv <- jsonlite::fromJSON(f_prov, simplifyVector = FALSE)
  if (!identical(pv$source_sha256, h0[["xls"]]) || !identical(pv$outputs[["manual_panel_extracted.csv"]], h0[["extracted"]])) stop("추출 출처 JSON 이 원본 XLS · 추출 CSV hash 와 맞지 않는다.")
  bl <- jsonlite::fromJSON(a$baseline)                                                                  # 상대경로 = W 기준, 절대경로 (main repo 3 개) = 그대로 (Codex 초안 검토 1)
  blf <- function() vapply(names(bl), function(f) { p <- if (startsWith(f, "/")) f else file.path(root, f); if (file.exists(p)) .nrv_sha(p) else "missing" }, "")
  bl0 <- blf(); if (!identical(unname(bl0), unname(unlist(bl)))) stop("기준선 파일 hash 가 codex_manual_panel_baseline.json 과 다르다: ", paste(names(bl)[bl0 != unlist(bl)], collapse = ", "))
  ck$add("inputs_pinned", TRUE, sprintf("xls %s, markers %s, umap %s, lossless %s, object %s, calc_manifest %s == pinned; extraction provenance matches the xls and extracted CSV hashes; %d baseline files unchanged",
                                        substr(h0[["xls"]], 1, 12), substr(h0[["markers"]], 1, 12), substr(h0[["umap"]], 1, 12), substr(h0[["lossless"]], 1, 12), substr(h0[["object"]], 1, 12),
                                        substr(h0[["calc_manifest"]], 1, 12), length(bl)))
  # ── (2) 두 arm 패널 ──
  ext <- rdc(f_ext); cr <- utils::read.delim(a$corrections, colClasses = "character", quote = "", na.strings = character(0), fileEncoding = "UTF-8")
  mp <- gsp_manual_panel(ext, cr); mu <- mp$used
  cur <- gsh_panels(gsh_read_markers(a$markers))$used; cur <- cur[cur$model == "slide2", , drop = FALSE]
  arms <- list(current_slide2 = data.frame(node_id = cur$node_id, display_label = cur$display_label, gene = cur$gene, sign = cur$sign, stringsAsFactors = FALSE),
               manual_original = data.frame(node_id = mu$node_id, display_label = mu$display_label, gene = mu$gene, sign = mu$sign, stringsAsFactors = FALSE))
  sh <- pv$sheets; nn <- vapply(sh, function(s) as.integer(s$n_nonempty_cells), 0L); names(nn) <- vapply(sh, function(s) s$sheet, "")
  hd <- lapply(sh, function(s) unlist(s$header)); ok_ext <- identical(sort(names(nn)), sort(names(GSP_SHEETS))) && identical(hd[[1]], hd[[2]]) && length(hd[[1]]) == 35L &&
    all(vapply(names(GSP_SHEETS), function(s) sum(ext$sheet == s) == nn[[s]] && all(as.integer(ext$sign[ext$sheet == s]) == GSP_SHEETS[[s]]), TRUE)) &&
    all(vapply(sh, function(s) as.integer(s$n_nontext_cells) == 0L, TRUE)) && setequal(unique(ext$header_raw), hd[[1]])
  ck$add("extraction_complete", ok_ext, sprintf("%d cells extracted (%s); headers identical in both sheets (%d types); all cells text; no cell dropped by the extraction",
                                                nrow(ext), paste(names(nn), nn, sep = " ", collapse = " · "), length(hd[[1]])))
  okc <- cr[cr$status %in% c("alias", "typo"), ]; chg <- mp$prov$symbol != mp$prov$gene
  ok_nm <- all(mp$prov$symbol[chg] %in% okc$source_symbol) && all(mp$prov$gene[chg] == okc$target_symbol[match(mp$prov$symbol[chg], okc$source_symbol)]) &&
    !any(mp$prov$symbol[!chg] %in% okc$source_symbol)
  ck$add("names_corrected_allowed_only", ok_nm, sprintf("%d cells renamed only by alias / typo rows (%s); unresolved kept as written (%s); %d cells with names not in the table kept",
                                                        sum(chg), paste(sort(unique(paste0(mp$prov$symbol[chg], ">", mp$prov$gene[chg]))), collapse = ", "),
                                                        paste(sort(unique(mp$prov$symbol[grepl("^kept_unresolved", mp$prov$correction)])), collapse = ", "), sum(mp$prov$correction == "kept_not_in_table")))
  need_cf <- c("MPP SLAMF1", "CMP KLF1", "pDC ITGAX", "Pro B cell MS4A1", "Cytotoxic CD8 CD28")                    # 인계서가 이름을 든 충돌 (최소)
  ok_dc <- !anyDuplicated(paste(mu$node_id, mu$gene, mu$sign)) && all(need_cf %in% paste(trimws(mp$conflicts$header_raw), mp$conflicts$gene)) &&
    all(vapply(seq_len(nrow(mp$conflicts)), function(i) all(c(1L, -1L) %in% mu$sign[mu$node_id == mp$conflicts$header_raw[i] & mu$gene == mp$conflicts$gene[i]]), TRUE)) &&
    all(c("Mono", "CD14 Mono", "CD16 Mono") %in% mu$node_id)
  ck$add("dedup_and_conflicts_preserved", ok_dc, sprintf("%d cells -> %d entries after removing %d same type / gene / sign duplicates; %d +/- conflicts kept on both sides (%s); Mono, CD14 Mono, CD16 Mono kept",
                                                         nrow(mp$prov), nrow(mu), sum(mp$prov$duplicate_dropped), nrow(mp$conflicts), paste(sprintf("%s %s", trimws(mp$conflicts$header_raw), mp$conflicts$gene), collapse = ", ")))
  mdir <- file.path(out, "manual_input"); dir.create(mdir, showWarnings = FALSE)
  for (p in c(a$xls, f_ext, f_prov, a$corrections)) { q <- file.path(mdir, basename(p)); if (!file.copy(p, q, overwrite = FALSE) || !identical(.nrv_sha(q), .nrv_sha(p))) stop("입력 사본이 원본과 다르다: ", basename(p))
    Sys.chmod(q, "0444"); files <- c(files, q) }
  # ── (3) 자료 (R75 와 같다) ──
  um <- readRDS(a$umap); ls <- readRDS(f_ls); cells <- ls$cells; cl <- as.character(ls$clusters); donor <- as.character(ls$donor); n <- length(cells)
  U <- um$U; mask <- um$subset_mask; sl <- as.character(um$subset_labels); clusters <- unique(rdc(f_cc)$cluster)
  ro <- gsx_open_rds(a$object); obj <- ro$object; md <- gsx_meta(obj); cn <- gsx_counts(obj, a$assay); X <- cn$counts; rm(obj, ro); invisible(gc())
  ok_ids <- identical(um$cells, cells) && identical(as.character(um$new_cluster), cl) && identical(colnames(X), cells) && identical(rownames(md), cells) && length(clusters) == 26L &&
            setequal(clusters, unique(cl)) && identical(dim(U), c(n, 2L)) && identical(rownames(U), cells) && all(is.finite(U)) && length(mask) == n && is.logical(mask)
  if (!ok_ids) stop("세포 ID · 순서 · 군집 · 좌표가 입력끼리 맞지 않는다.")
  lib <- Matrix::colSums(X); if (!identical(unname(as.numeric(lib)), ls$library)) stop("library 가 lossless 와 다르다.")
  ck$add("cells_clusters_coordinates_library", TRUE, sprintf("%d cell IDs in the same order (UMAP RDS, lossless, raw counts, meta.data); %d clusters; U %d x 2 finite and drawn as stored; subset %d cells; all-gene library identical to lossless",
                                                             n, length(clusters), n, sum(mask)))
  rm(md); invisible(gc()); log("inputs fixed:", n, "cells,", length(clusters), "clusters | manual", nrow(mu), "entries in", length(unique(mu$node_id)), "types | current", nrow(cur), "rows in", length(unique(cur$node_id)), "types")
  # ── (4) 집계 · 공통 대조 · 점수 (저자 라벨 없음) ──
  agg <- gsh_aggregate(X, lib, cl, donor, clusters); bins <- gsh_bins(agg$gene_mean)
  active <- unique(c(arms$current_slide2$gene, arms$manual_original$gene)); meas <- intersect(active, names(bins)); ct <- gsh_controls(bins, meas, exclude = active)
  sc <- do.call(rbind, lapply(GSP_ARMS, function(m) data.frame(arm = m, gsh_type_scores(agg, arms[[m]], ct$controls), stringsAsFactors = FALSE)))
  rk <- lapply(stats::setNames(GSP_ARMS, GSP_ARMS), function(m) gsh_rank(sc[sc$arm == m, , drop = FALSE], "score"))
  rd <- lapply(stats::setNames(GSP_ARMS, GSP_ARMS), function(m) gsh_rank(sc[sc$arm == m, , drop = FALSE], "donor_equal_score"))
  sc2 <- sc; names(sc2)[names(sc2) == "arm"] <- "model"; pk <- gsh_type_peaks(sc2); names(pk)[names(pk) == "model"] <- "arm"
  vr <- do.call(rbind, lapply(GSP_VERIFY, function(v) { z <- arms[[v[["arm"]]]]; z <- z[z$node_id == v[["type"]], , drop = FALSE]
    gp <- intersect(z$gene[z$sign == 1L], meas); gn <- intersect(z$gene[z$sign == -1L], meas); c0 <- gsh_cell_level(X, lib, cl, clusters, gp, gn, ct$controls)
    s0 <- sc[sc$arm == v[["arm"]] & sc$type == v[["type"]], , drop = FALSE]; s0 <- s0[match(clusters, s0$cluster), ]
    data.frame(arm = v[["arm"]], type = v[["type"]], cluster = clusters, pos_aggregated = s0$pos_enrichment, pos_cell_level = c0$pos_cell_mean, neg_aggregated = s0$neg_enrichment,
               neg_cell_level = c0$neg_cell_mean, stringsAsFactors = FALSE) }))
  vd <- max(abs(c(vr$pos_aggregated - vr$pos_cell_level, vr$neg_aggregated - vr$neg_cell_level)), na.rm = TRUE); rm(X); invisible(gc())
  cg <- unlist(ct$controls); bsz <- table(bins)
  ok_ct <- !any(cg %in% active) && all(vapply(names(ct$controls), function(g) all(bins[ct$controls[[g]]] == bins[[g]]), TRUE)) && all(ct$table$n_controls == pmin(GSH_N_CTRL, ct$table$n_candidates_in_bin)) &&
           length(bsz) == GSH_N_BINS && (max(bsz) - min(bsz)) <= 1L && setequal(names(ct$controls), meas)
  ck$add("controls_common_exclude_union", ok_ct, sprintf("one control set (seed %d, %d bins) for %d measured marker genes of both arms; no control is any of the %d active genes of either arm; shortfall %d",
                                                         GSH_SEED, length(bsz), length(meas), length(active), sum(ct$table$shortfall)))
  npos <- sc$n_pos_measured > 0; nty <- vapply(GSP_ARMS, function(m) length(unique(arms[[m]]$node_id)), 0L)
  ok_tab <- all(vapply(GSP_ARMS, function(m) sum(sc$arm == m) == length(clusters) * nty[[m]], TRUE)) && all(is.finite(sc$score[npos])) && all(is.na(sc$score[!npos]))
  ck$add("score_tables_complete", ok_tab, sprintf("%s rows = 26 clusters x types; scores finite where >= 1 positive marker measured, NA otherwise (%d rows)",
                                                  paste(sprintf("%s %d (%d types)", GSP_ARMS, as.integer(table(sc$arm)[GSP_ARMS]), nty), collapse = ", "), sum(!npos)))
  ck$add("cell_level_recompute", is.finite(vd) && vd <= 1e-10, sprintf("cluster means of per-cell enrichments for %s equal the aggregated values in all 26 clusters; max |diff| %.3g",
                                                                       paste(vapply(GSP_VERIFY, function(v) paste(v, collapse = ":"), ""), collapse = ", "), vd))
  # 계산 표 고정
  mp_used <- mu; mp_used$measured <- mp_used$gene %in% meas; cur_used <- cur; cur_used$measured <- cur_used$gene %in% meas
  ctab <- ct$table; ctab$mean_log1p_cp10k <- unname(agg$gene_mean[ctab$gene]); ctab$arms <- vapply(ctab$gene, function(g) paste(GSP_ARMS[vapply(GSP_ARMS, function(m) g %in% arms[[m]]$gene, TRUE)], collapse = "|"), "")
  unm <- do.call(rbind, lapply(GSP_ARMS, function(m) { z <- arms[[m]]; u <- unique(z[!z$gene %in% meas, c("node_id", "gene", "sign")])
    if (!nrow(u)) return(NULL); data.frame(arm = m, u, reason = ifelse(u$gene %in% names(bins), "no controls", "symbol not in the raw-count genes"), stringsAsFactors = FALSE) }))
  if (is.null(unm)) unm <- data.frame(arm = character(0), node_id = character(0), gene = character(0), sign = character(0), reason = character(0), stringsAsFactors = FALSE)
  wc(sc, "cluster_type_scores.csv"); wc(pk, "type_top_clusters.csv"); wc(mp$prov, "manual_panel_provenance.csv"); wc(mp_used, "manual_panel_used.csv"); wc(cur_used, "current_slide2_panel_used.csv")
  wc(mp$conflicts, "manual_panel_sign_conflicts.csv"); wc(mp$display_map, "display_name_map.csv"); wc(unm, "markers_unmeasured.csv"); wc(ctab, "control_genes.csv"); wc(vr, "check_cell_level_recompute.csv")
  cmf <- file.path(out, "calc_manifest.json"); calc_files <- files
  writeLines(jsonlite::toJSON(list(phase = "panels, common controls, type scores and ranks fixed before author labels are read", written_at = stamp(), inputs = as.list(h0),
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  log("calculation fixed:", length(calc_files), "files | top1 status:", paste(vapply(GSP_ARMS, function(m) paste(m, paste(names(table(rk[[m]]$status)), table(rk[[m]]$status), collapse = " ")), ""), collapse = " ; "))
  # ── (5) 사후: 저자 구성 · 나란히 표 · 핵심 유형 · R75 대조군 차이 ──
  au <- as.character(um$author_label); comp <- gsm_author_composition(cl, au, clusters)
  a1 <- vapply(clusters, function(k) comp$author_label[comp$cluster == k][1], ""); f1 <- vapply(clusters, function(k) comp$frac_of_cluster[comp$cluster == k][1], 0)
  a2 <- vapply(clusters, function(k) { z <- comp[comp$cluster == k, ]; if (nrow(z) > 1L) z$author_label[2] else "" }, ""); f2 <- vapply(clusters, function(k) { z <- comp[comp$cluster == k, ]; if (nrow(z) > 1L) z$frac_of_cluster[2] else NA_real_ }, 0)
  ckey <- function(m, lab) { if (m == "manual_original") vapply(lab, function(x) { y <- strsplit(x, " / ", fixed = TRUE)[[1]]; paste(ifelse(y %in% names(GSP_COLOUR_MAP), GSP_COLOUR_MAP[y], y), collapse = " / ") }, "", USE.NAMES = FALSE) else lab }
  side <- data.frame(cluster = clusters, n_cells = as.integer(agg$n_cells[clusters]), n_subset_cells = as.integer(table(factor(cl[mask], levels = clusters))),
                     author_top_label = unname(a1), author_top_frac = unname(f1), author_second_label = unname(a2), author_second_frac = unname(f2), stringsAsFactors = FALSE)
  for (m in GSP_ARMS) { r <- rk[[m]]; e <- rd[[m]]; p <- if (m == "manual_original") "manual" else "current"
    side[[paste0(p, "_status")]] <- r$status; side[[paste0(p, "_top1")]] <- r$top1_label; side[[paste0(p, "_top1_score")]] <- r$top1_score; side[[paste0(p, "_top2")]] <- r$top2_label
    side[[paste0(p, "_top2_score")]] <- r$top2_score; side[[paste0(p, "_gap")]] <- r$gap; side[[paste0(p, "_donor_equal_top1")]] <- e$top1_label; side[[paste0(p, "_donor_equal_status")]] <- e$status
    side[[paste0(p, "_display_key")]] <- ifelse(r$status == "undetermined", "undetermined", ckey(m, r$top1_label)) }
  side$same_display_key <- side$current_display_key == side$manual_display_key
  side$note <- rep("provisional top-score candidates; author labels post hoc; display keys use display_name_map.csv (figure colours and this column only)", nrow(side))
  wc(side, "cluster_side_by_side.csv"); comp$note <- rep("post hoc only", nrow(comp)); wc(comp, "posthoc_author_composition.csv")
  core <- do.call(rbind, lapply(GSP_ARMS, function(m) { z <- pk[pk$arm == m & pk$type %in% GSP_CORE[[m]], , drop = FALSE]; z <- z[order(match(z$type, GSP_CORE[[m]])), , drop = FALSE]
    z$peak_cluster_author_top <- vapply(z$peak_clusters, function(s) if (!nzchar(s)) "" else paste(a1[strsplit(s, "|", fixed = TRUE)[[1]]], collapse = "|"), "")
    z$peak_cluster_top1_in_arm <- vapply(z$peak_clusters, function(s) if (!nzchar(s)) "" else paste(rk[[m]]$top1_label[match(strsplit(s, "|", fixed = TRUE)[[1]], rk[[m]]$cluster)], collapse = "|"), "")
    z$clusters_with_this_top1 <- vapply(z$type, function(ty) paste(rk[[m]]$cluster[rk[[m]]$status == "candidate" & rk[[m]]$top1 == ty], collapse = ","), ""); z }))
  wc(core, "core_type_top_clusters.csv")
  r75a <- rdc(f_r75a); r75c <- rdc(f_r75c); g75 <- intersect(r75c$gene, arms$current_slide2$gene)
  same_ctrl <- vapply(g75, function(g) identical(sort(strsplit(r75c$controls[r75c$gene == g], "|", fixed = TRUE)[[1]]), sort(ct$controls[[g]])), TRUE)
  cvr <- data.frame(cluster = clusters, r75_slide2_top1 = r75a$slide2_top1[match(clusters, r75a$cluster)], this_current_slide2_top1 = rk$current_slide2$top1_label, stringsAsFactors = FALSE)
  cvr$same_top1 <- cvr$r75_slide2_top1 == cvr$this_current_slide2_top1
  cvr$note <- rep(sprintf("R75 controls excluded all listed markers of three networks; here the union of both arms' active genes. Identical control sets for %d / %d current_slide2 genes", sum(same_ctrl), length(g75)), nrow(cvr))
  wc(cvr, "current_vs_r75_controls.csv")
  ck$add("author_labels_after_calculation", identical(unname(vapply(calc_files, .nrv_sha, "")), unname(unlist(jsonlite::fromJSON(cmf)$files[rel(calc_files)]))),
         "author labels read only after calc_manifest.json; no cutoff or feature from author labels; fixed calculation files unchanged afterwards")
  # ── (6) 그림 ──
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; acols <- .gsu_cols(au, sl); fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  tops <- lapply(stats::setNames(GSP_ARMS, GSP_ARMS), function(m) { r <- rk[[m]]; data.frame(cluster = r$cluster, status = r$status, label = r$top1_label,
    key = ifelse(r$status == "undetermined", "undetermined", ckey(m, r$top1_label)), stringsAsFactors = FALSE) })
  ukeys <- unique(c(tops$current_slide2$key, tops$manual_original$key)); ukeys <- setdiff(ukeys, "undetermined")
  ord <- unique(c(ckey("current_slide2", unique(arms$current_slide2$display_label)), ckey("manual_original", unique(arms$manual_original$display_label)))); ukeys <- ukeys[order(match(ukeys, ord), ukeys)]
  pal <- c(GSH_CAT_COLS, grDevices::hcl.colors(max(1L, length(ukeys) - length(GSH_CAT_COLS)), "Dark 3")); kcol <- c(stats::setNames(pal[seq_along(ukeys)], ukeys), undetermined = GSH_UNDET_COL)
  sub <- sprintf(if (ko) "R %s base graphics · 좌표 = %s (재계산 없음) · 점수 = 배경 보정 module score, 두 판 공통 대조 유전자 (합집합 제외, seed %d)" else "R %s base graphics · coordinates = %s (not recomputed) · score = background-corrected module score, common controls for both arms (union excluded, seed %d)",
                 as.character(getRversion()), basename(dirname(a$umap)), GSH_SEED)
  draw <- function(path, mk = NULL) gsp_draw_three(path, U, au, acols, cl, tops, kcol, ntypes = nty, mask = mk, subset_labels = sl, label_min_frac = GSH_LABEL_MIN_FRAC, sub = sub, ko = ko, fam = fam)
  fp1 <- file.path(fig, "fig1_umap_full_author_current_manual.png"); r1 <- draw(fp1)
  fp2 <- file.path(fig, "fig2_umap_lineage_author_current_manual.png"); r2 <- draw(fp2, mask)
  files <- c(files, fp1, fp2); ft <- rbind(r1$text, r2$text); wc(ft, "figure_text.csv"); log("figures written")
  same_xy <- function(A, B) identical(as.vector(A), as.vector(B)) && identical(dim(A), dim(B)) && identical(dimnames(A), dimnames(B))
  colok <- all(vapply(c(r1$colours, r2$colours), function(z) all(z$n_colours == 1L) && identical(unname(kcol[z$key]), z$colour), TRUE))
  ck$add("figures_same_coordinates", same_xy(r1$coords, U) && same_xy(r2$coords, U[mask, , drop = FALSE]) && colok && !any(ft$text %in% c("NA", "NaN", "Inf")) &&
           !any(grepl("→|->|−|ï", ft$text)),
         sprintf("figure 1 all %d cells, figure 2 all %d subset cells at the stored coordinates (values and IDs identical); one colour per cluster per arm; no NA text, arrows, U+2212 or i-diaeresis glyphs", n, sum(mask)))
  # ── (7) 기록 · 요약 · 공개 ──
  h1 <- vapply(inp, .nrv_sha, ""); bl1 <- blf()
  ck$add("existing_inputs_unchanged", identical(unname(h1), unname(h0)) && identical(unname(bl1), unname(unlist(bl))),
         sprintf("%d inputs re-hashed after the run (original xls, corrections, extraction, markers, coordinates, lossless, object, R75 tables, baseline json) and %d baseline files unchanged", length(inp), length(bl)))
  wc(data.frame(input = names(inp), path = unname(inp), sha256 = unname(h1), stringsAsFactors = FALSE), "input_sha256.csv"); wc(snap, "source_snapshot.csv")
  pkg <- c("Matrix", "digest", "jsonlite")
  wc(data.frame(item = c("R_executable", "R", "platform", "libPaths", "direct_packages", "base_packages", "rng", "seed", "extraction_python", "extraction_xlrd", "extraction_command", "graphics",
                         "title_font", "command", "working_dir", "started_at", "calc_fixed_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pkg, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), "stats, utils, methods, grDevices, graphics",
                          "Mersenne-Twister / Inversion / Rejection", as.character(GSH_SEED), paste(pv$python, pv$python_version), pv$xlrd_version, paste(unlist(pv$command), collapse = " "),
                          "R base graphics (grDevices::png, res 100)", if (ko) MBVR_TITLE_FAMILY else "default", paste(c("76_gse253355_manual_panel.R", argv), collapse = " "), getwd(),
                          format(t0, "%Y-%m-%dT%H:%M:%S%z"), jsonlite::fromJSON(cmf)$written_at, stamp()), stringsAsFactors = FALSE), "environment.csv")
  tbs <- function(m) { r <- rk[[m]]; sprintf("%s: 1 위 후보 %d · 동률 %d · 미판정 %d", m, sum(r$status == "candidate"), sum(r$status == "tied"), sum(r$status == "undetermined")) }
  keyc <- c("6", "8", "17", "7", "12", "15", "16")
  kc <- paste(vapply(keyc, function(k) { i <- match(k, side$cluster); sprintf("%s (%s) %s / %s", k, side$author_top_label[i], side$current_top1[i], side$manual_top1[i]) }, ""), collapse = " · ")
  cor1 <- function(ty_c, ty_m) { x1 <- core[core$arm == "current_slide2" & core$type == ty_c, ]; x2 <- core[core$arm == "manual_original" & core$type == ty_m, ]
    sprintf("%s 점수 최고 군집 현재 %s (%s) · 수기 %s (%s); 1 위인 군집 현재 [%s] · 수기 [%s]", ty_m, x1$peak_clusters, x1$peak_cluster_author_top, x2$peak_clusters, x2$peak_cluster_author_top,
            x1$clusters_with_this_top1, x2$clusters_with_this_top1) }
  sm <- c("# GSE253355 — 최초 수기 마커 패널 (35 유형) 대 현재 slide2 + T (30 유형) 비교 (요약)", "",
          "- 사용: 1) python (xlrd) scripts/r76_extract_manual_panel_xls.py 로 원본 XLS 형식 추출  2) run-r dcta R/branch_biology/76_gse253355_manual_panel.R --out-dir <새 폴더>.",
          sprintf("- 원본: final_unique_gene.xls (sha256 %s…) 두 sheet %d 칸, alias · typo 교정 %d 칸, 중복 %d 제거 → %d 항목, 양 · 음 충돌 %d 쌍은 양쪽 보존 (%s).", substr(h0[["xls"]], 1, 12), nrow(ext),
                  sum(chg), sum(mp$prov$duplicate_dropped), nrow(mu), nrow(mp$conflicts), paste(sprintf("%s %s", trimws(mp$conflicts$header_raw), mp$conflicts$gene), collapse = ", ")),
          sprintf("- 미측정: %s. unresolved 이름 (%s) 은 원문 유지.", paste(vapply(GSP_ARMS, function(m) sprintf("%s %d", m, sum(unm$arm == m)), ""), collapse = " · "),
                  paste(sort(unique(mp$prov$symbol[grepl("^kept_unresolved", mp$prov$correction)])), collapse = ", ")),
          sprintf("- 점수: R75 와 같은 배경 보정 module score. 대조 유전자는 두 판 합집합 %d 개를 빼고 한 번 만들었다 (R75 와 다름: 현재판 유전자 %d / %d 만 같은 대조군; 현재판 1 위가 R75 와 같은 군집 %d / 26).",
                  length(active), sum(same_ctrl), length(g75), sum(cvr$same_top1)),
          paste0("- ", tbs("current_slide2"), " | ", tbs("manual_original"), sprintf(". 표시명 대응 기준 두 판 1 위가 같은 군집 %d / 26.", sum(side$same_display_key))),
          paste0("- 주요 군집 (저자 최다) 현재 / 수기: ", kc, "."),
          paste0("- ", cor1("HSC", "HSC"), "."), paste0("- ", cor1("MPP", "MPP"), "."), paste0("- ", cor1("Neutrophil", "Neutrophils"), "."),
          paste0("- ", cor1("Classical_Monocyte", "CD14 Mono"), "."),
          "- 해석: 학습 · 독립 검증이 없어 어느 판이 더 정확하다고 말하지 않는다. 유형 수 (35 vs 30) · 세분화 · 계층 이름 (Mono 와 CD14 / CD16 Mono) 이 달라 단순 비교가 아니다.",
          "- 한계: MSC · 내피 등 패널 없는 계통도 혈액 유형 후보가 붙는다. PM 등 성숙 단계는 다루지 않았다. 저자 라벨은 계산 고정 뒤 사후 대조에만 썼다.",
          "- 표: cluster_side_by_side · cluster_type_scores · core_type_top_clusters · manual_panel_provenance / used / sign_conflicts · markers_unmeasured · control_genes · display_name_map.")
  pub_list <- c("summary_ko.md", "cluster_side_by_side.csv", "cluster_type_scores.csv", "type_top_clusters.csv", "core_type_top_clusters.csv", "posthoc_author_composition.csv",
                "manual_panel_provenance.csv", "manual_panel_used.csv", "current_slide2_panel_used.csv", "manual_panel_sign_conflicts.csv", "display_name_map.csv", "markers_unmeasured.csv",
                "control_genes.csv", "current_vs_r75_controls.csv", "check_cell_level_recompute.csv", "checks.csv", "figure_text.csv", "environment.csv", "input_sha256.csv",
                file.path("manual_input", basename(c(a$xls, f_ext, f_prov, a$corrections))), file.path("figures", basename(c(fp1, fp2))))
  ck$add("tables_roundtrip_ok", all(crt), sprintf("%d CSVs re-read with identical dimensions, names and cell text", length(crt)))
  ck$add("public_copy_aggregate_only", !any(grepl("lossless|coordinates|aggregates|\\.rds$|\\.gz$", pub_list)), sprintf("%d public files: aggregate tables, PNG, original XLS / extraction / corrections copies, checks; no cell-level data", length(pub_list)))
  cks <- ck$table(); sm <- c(sm, sprintf("- 검사 %d / %d 통과 (checks.csv).", sum(cks$pass), nrow(cks)))
  writeLines(sm, file.path(out, "summary_ko.md")); files <- c(files, file.path(out, "summary_ko.md")); wc(cks, "checks.csv"); allf <- c(files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), inputs = as.list(h0), files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  pub <- a$public; for (dd in unique(dirname(file.path(pub, pub_list)))) dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(file.path(out, pub_list), file.path(pub, pub_list), overwrite = FALSE))) stop("공개 사본 복사 실패")
  writeLines(jsonlite::toJSON(list(source_output = normalizePath(out), inputs = as.list(h0), files = as.list(stats::setNames(vapply(file.path(pub, pub_list), .nrv_sha, ""), pub_list))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(pub, "public_manifest.json"))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .gsp_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.gsp_bb, .f))
  quit(status = gsp_main(commandArgs(trailingOnly = TRUE)))
}
