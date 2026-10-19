# =============================================================================
# 62_nso_connectivity_audit.R — 호중구 단계 marker 검출 세포가 원 세포 그래프로 이어지는가 (한정 감사: 계산 · 실행 · 작은 그림)
# =============================================================================
#
# 계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_connectivity_instruction.md (2026-10-19). 입력 = NSO attempt_02 (manifest 고정) + NRV edges.rds.
#   - 검출 집합: 단계 marker 중 log > 0 (= raw count > 0) 인 유전자가 1 개 이상(그림 1 색칠 기준) · 2 개 이상(민감도 참고). 원 RNA 를 다시 추정하지 않고
#     저장된 log 열로 정의하고 저장된 검출 수와 대조한다. 단계 겹침은 그대로 두고, 배타 주석 · argmax · 새 cutoff 를 만들지 않는다.
#   - 연결 성분: 전체 그래프와 단계별 유도 부분그래프(그 단계 세포끼리의 원 edge 만). 같은 block 조건 = block 이 같은 원 edge 만 (단계 안 · 사이 모두).
#   - 인접 쌍 A → B: A–B 원 edge 수 · 양끝 고유 세포 · 같은 block edge · 겹침 세포 · A-only ↔ B-only edge. 방향은 기준 순서일 뿐 시간 검증이 아니다.
#   - 순서 도달: 단계별 유도 성분을 단위로, 연속 단계 사이 실제 원 edge 로만 다음 층으로 간다. 무관한 세포를 거친 우회 · 가상 edge ·
#     같은 세포라는 이유만의 zero-hop 전진은 없다. 단계가 겹치므로 같은 원 세포가 실제 edge 를 거쳐 다음 층에 다시 나올 수 있다(n_reappearing_cells).
#     6 개의 배타 세포 상태나 분화 경로가 아니다. 어떤 결과도 성숙 · RNA 연속성 · 확률 · 정확도 입증이 아니다.
# 실행: run-r dcta R/branch_biology/62_nso_connectivity_audit.R --in-dir <NSO attempt_02> --o-root <O> --out-dir <N>/results/neut_stage_observation_2026-10-19/connectivity_audit_01
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph) })

NCA_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "27_run_nrv.R",
                 "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R",
                 "58_nso_inputs.R", "62_nso_connectivity_audit.R")
NCA_PINNED_MANIFEST <- "4c31d681b0da2fc11b747dd865ce68431a10d161896f03515557b2d54d3c36a7"   # NSO attempt_02 manifest.json (지시서 고정값)
NCA_NRV_P92 <- "results/neutrophil_route_validation_2026-10-13/attempt_01_p92"                # O 기준: 원 세포 그래프의 양수 상삼각 edge (edges.rds)
NCA_N_CELLS <- 263159L
NCA_CRITERIA <- c(ge1 = 1L, ge2 = 2L)                                                         # 1 개 이상 = 그림 1 색칠 기준, 2 개 이상 = 민감도 참고 (분류 기준 아님)
NCA_GRAPHS <- c(all_edges = FALSE, same_block = TRUE)


# ── 1. 계산 (순수) ───────────────────────────────────────────────────────────
#' 검출 집합: 세포 x 단계 논리 행렬. markers = 이름 붙은 list (단계 → marker 유전자), L = 세포 x 유전자 log 희소행렬
nca_detection_sets <- function(L, markers, min_markers) {
  miss <- setdiff(unlist(markers), colnames(L)); if (length(miss)) stop("marker 유전자가 log 행렬에 없다: ", paste(unique(miss), collapse = ", "))
  S <- matrix(FALSE, nrow(L), length(markers), dimnames = list(NULL, names(markers)))
  for (s in names(markers)) { g <- markers[[s]]; if (length(g)) S[, s] <- as.numeric(Matrix::rowSums(L[, g, drop = FALSE] > 0)) >= min_markers }
  S
}

.nca_edges <- function(edges, block = NULL, same_block = FALSE) {
  if (!same_block) return(edges)
  if (is.null(block)) stop("same_block 조건에는 block 이 필요하다.")
  edges[block[edges$i] == block[edges$j], , drop = FALSE]
}

#' 연결 성분: keep 이 없으면 전체 그래프, 있으면 keep 세포끼리의 원 edge 만 쓴 유도 부분그래프. membership 은 keep 밖 세포에서 NA.
nca_components <- function(n, edges, keep = NULL, block = NULL, same_block = FALSE) {
  e <- .nca_edges(edges, block, same_block); if (is.null(keep)) keep <- rep(TRUE, n)
  if (length(keep) != n) stop("keep 길이가 세포 수와 다르다.")
  idx <- which(keep); loc <- integer(n); loc[idx] <- seq_along(idx); sel <- keep[e$i] & keep[e$j]
  g <- igraph::make_empty_graph(n = length(idx), directed = FALSE)
  if (any(sel)) g <- igraph::add_edges(g, as.vector(rbind(loc[e$i[sel]], loc[e$j[sel]])))
  m <- if (length(idx)) igraph::components(g)$membership else integer(0); memb <- rep(NA_integer_, n); memb[idx] <- m
  list(membership = memb, n_cells = length(idx), n_edges = as.integer(sum(sel)), n_components = length(unique(m)),
       largest = if (length(m)) as.integer(max(tabulate(m))) else 0L, n_isolated = as.integer(sum(igraph::degree(g) == 0)))
}

#' 인접 쌍 A → B 의 원 edge 집계 (한 edge 는 한 번만 센다)
nca_pair_edges <- function(edges, A, B, block, same_block = FALSE) {
  e <- .nca_edges(edges, block, same_block); i <- e$i; j <- e$j
  ab <- A[i] & B[j]; ba <- B[i] & A[j]; hit <- ab | ba; Ao <- A & !B; Bo <- B & !A; oo <- (Ao[i] & Bo[j]) | (Bo[i] & Ao[j])
  data.frame(n_cells_A = as.integer(sum(A)), n_cells_B = as.integer(sum(B)), n_overlap_cells = as.integer(sum(A & B)), n_edges = as.integer(sum(hit)),
             n_endpoints_A = length(unique(c(i[ab], j[ba]))), n_endpoints_B = length(unique(c(j[ab], i[ba]))), n_same_block_edges = as.integer(sum(hit & block[i] == block[j])),
             n_edges_Aonly_Bonly = as.integer(sum(oo)), n_edges_with_overlap_endpoint = as.integer(sum(hit) - sum(oo)))
}

#' 순서 도달: 첫 단계의 모든 검출 세포에서 시작해, 직전 층 도달 세포와 실제 원 edge 로 닿은 다음 단계 세포의 유도 성분 전체로 전진한다.
#' 반환: table(단계별 세포 · 성분 · 도달 세포 · 도달 성분 · 다시 나온 세포), first_break ("" = 끝까지 도달)
nca_ordered_reach <- function(n, edges, sets, block = NULL, same_block = FALSE) {
  e <- .nca_edges(edges, block, same_block); st <- colnames(sets)
  comps <- lapply(seq_along(st), function(k) nca_components(n, e, keep = sets[, k]))
  R <- sets[, 1]; fb <- ""; rows <- vector("list", length(st))
  for (k in seq_along(st)) { reap <- 0L
    if (k > 1L) { S1 <- sets[, k]
      hit <- unique(c(e$j[R[e$i] & S1[e$j]], e$i[R[e$j] & S1[e$i]]))                       # 실제 edge 로 닿은 이번 단계 세포 (자기 자신은 edge 가 아니다)
      Rn <- S1 & comps[[k]]$membership %in% unique(comps[[k]]$membership[hit])
      reap <- as.integer(sum(Rn & R)); if (!any(Rn) && !nzchar(fb)) fb <- paste(st[k - 1L], st[k], sep = ">"); R <- Rn }
    rows[[k]] <- data.frame(stage = st[k], n_cells = as.integer(sum(sets[, k])), n_components = comps[[k]]$n_components, n_reachable_cells = as.integer(sum(R)),
                            n_reachable_components = length(unique(comps[[k]]$membership[R])), n_reappearing_cells = reap, reached = any(R), stringsAsFactors = FALSE) }
  list(table = do.call(rbind, rows), first_break = fb)
}


# ── 2. 그림 (작은 요약 1 장) ─────────────────────────────────────────────────
#' log10(n + 1) 축의 눈금: 원 세포 수 0, 1, 10, 100, … 을 실제 좌표 log10(값 + 1) 에 둔다 (Codex 검토: 0 이 1 자리에 보이지 않게)
.nca_log_ticks <- function(ymax) {
  v <- c(0, 10^(0:max(0, floor(ymax)))); v <- v[log10(v + 1) <= ymax + 1e-12]
  data.frame(value = v, at = log10(v + 1), label = .nca_fmt_n(v), stringsAsFactors = FALSE)
}
#' 6 단계 축: 기준 x 그래프 4 패널. 회색 = 검출 세포, 색 = 순서 도달 세포 (log10 축, 막대 위 실제 수), 첫 끊김 표시
nca_draw_summary <- function(path, rt, ko, fam) {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; conds <- unique(rt[, c("criterion", "graph")])
  grDevices::png(path, width = 2600, height = 1900, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::par(mfrow = c(2, 2), oma = c(4.6, 0, 3.6, 0), mar = c(4.2, 6, 5, 1.5), family = fam)
    for (r in seq_len(nrow(conds))) { z <- rt[rt$criterion == conds$criterion[r] & rt$graph == conds$graph[r], ]; st <- z$stage
      ymax <- log10(max(z$n_cells, 1) + 1) * 1.18; x <- seq_along(st)
      graphics::plot(NA, xlim = c(0.5, length(st) + 0.5), ylim = c(0, ymax), xaxt = "n", yaxt = "n", xlab = "",
                     ylab = if (ko) "세포 수 (축 = log10(n + 1), 눈금 = 원 세포 수)" else "cells (axis = log10(n + 1), ticks = raw counts)")
      graphics::axis(1, at = x, labels = st, cex.axis = 1.25); tk <- .nca_log_ticks(ymax); graphics::axis(2, at = tk$at, labels = tk$label, las = 1)
      graphics::rect(x - 0.38, 0, x + 0.38, log10(z$n_cells + 1), col = "grey85", border = "grey55")
      yr <- log10(z$n_reachable_cells + 1); graphics::rect(x - 0.24, 0, x + 0.24, yr, col = "#3B7DBF", border = NA)
      graphics::text(x, log10(z$n_cells + 1), .nca_fmt_n(z$n_cells), pos = 3, cex = 1.05)
      inb <- yr >= 0.45; lab <- .nca_fmt_n(z$n_reachable_cells); show <- z$n_cells > 0                     # 짧은 막대 · 0 은 막대 위 검정 글자 (검출 0 이면 회색 글자 하나만)
      if (any(inb & show)) graphics::text(x[inb & show], yr[inb & show] / 2, lab[inb & show], col = "white", font = 2, cex = 1.0)
      if (any(!inb & show)) graphics::text(x[!inb & show], yr[!inb & show], lab[!inb & show], col = "black", font = 2, cex = 1.0, pos = 3)
      fb <- z$first_break[1]
      if (nzchar(fb)) { k <- match(strsplit(fb, ">", fixed = TRUE)[[1]][2], st); graphics::text(k - 0.5, ymax * 0.93, "×", col = "#CB181D", cex = 3.2, font = 2)
        note(paste(conds$criterion[r], conds$graph[r]), "first_break", fb) } else note(paste(conds$criterion[r], conds$graph[r]), "first_break", if (ko) "끊김 없음" else "no break")
      ttl <- sprintf("%s · %s", if (conds$criterion[r] == "ge1") (if (ko) "marker 1 개 이상" else ">= 1 marker") else (if (ko) "marker 2 개 이상" else ">= 2 markers"),
                     if (conds$graph[r] == "all_edges") (if (ko) "모든 원 edge" else "all original edges") else (if (ko) "같은 block edge 만" else "same-block edges only"))
      sub <- if (nzchar(fb)) sprintf(if (ko) "첫 끊김 %s" else "first break %s", fb) else if (ko) "6 단계 끝까지 도달" else "reached all 6 stages"
      graphics::title(main = ttl, cex.main = 1.6, line = 2.4); graphics::mtext(sub, side = 3, line = 0.6, cex = 1.15, col = if (nzchar(fb)) "#CB181D" else "grey20")
      note(paste(conds$criterion[r], conds$graph[r]), "title", ttl); note(paste(conds$criterion[r], conds$graph[r]), "sub", sub) }
    main <- if (ko) "단계 marker 검출 세포의 원 그래프 순서 도달 (생물학적 DAG · 실제 세포 경로 아님)" else "Ordered reachability of stage-marker cells on the original graph (not a DAG or a cell path)"
    graphics::mtext(main, side = 3, outer = TRUE, line = 1.2, cex = 1.85, font = 2); note("figure", "title", main)
    foot <- if (ko) c("회색 = 그 단계 검출 세포, 파랑 = 첫 단계 모든 검출 세포에서 연속 단계 사이 실제 원 edge 로만 전진해 닿은 세포(단계 유도 성분 단위). × = 첫 끊김. 세로축 = log10(n + 1), 0 은 바닥선.",
                      "단계는 겹칠 수 있어 같은 세포가 실제 edge 를 거쳐 다음 층에 다시 나올 수 있다. 배타 세포 상태 · 분화 방향 · 확률 · 정확도가 아니다.")
            else c("Grey = detected cells; blue = cells reached from all first-stage cells using only real edges between consecutive stages (stage-induced components). x = first break. y = log10(n + 1), 0 at baseline.",
                   "Stages overlap, so a cell can reappear in the next layer through a real edge. Not exclusive states, directions, probabilities or accuracy.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.8 + (k - 1) * 1.3, cex = 1.05, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  nf$table()
}
.nca_fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)


# ── 3. 실행기 ────────────────────────────────────────────────────────────────
nca_main <- function(argv, root = NULL, pins = NCA_PINNED_MANIFEST) {
  a <- list(in_dir = NULL, o_root = NULL, out = NULL); keys <- c("--in-dir" = "in_dir", "--o-root" = "o_root", "--out-dir" = "out"); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (any(vapply(a, function(x) is.null(x) || is.na(x) || !nzchar(x), TRUE))) { cat("--in-dir, --o-root, --out-dir 가 필요하다.\n"); return(2L) }
  out <- a$out
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", NCA_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(NCA_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = NCA_SOURCES)
    get("nca_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

nca_build <- function(a, out, snap, log, argv, pins) {
  ck <- .mbvq_checks()
  # ── 입력 고정 ──
  mf <- file.path(a$in_dir, "manifest.json"); hm <- .nrv_sha(mf)
  if (!is.null(pins) && !identical(hm, pins)) stop("NSO 결과 manifest hash 가 지시서 고정값과 다르다: ", hm)
  inp <- c(nso_manifest = mf, nso_lossless = .mbvq_manifest_file(a$in_dir, "manifest.json", "lossless.rds"),
           nso_reference_stages = .mbvq_manifest_file(a$in_dir, "manifest.json", "reference_stages.tsv"),
           nso_input_sha256 = .mbvq_manifest_file(a$in_dir, "manifest.json", "input_sha256.tsv"),
           nrv_p92_edges = nrv_stage_file(file.path(a$o_root, NCA_NRV_P92), "structure", "edges.rds"))
  rec <- utils::read.delim(inp[["nso_input_sha256"]], stringsAsFactors = FALSE, colClasses = "character")
  if (!identical(.nrv_sha(inp[["nrv_p92_edges"]]), rec$sha256[rec$input == "nrv_p92_edges"])) stop("edges.rds hash 가 NSO 실행 기록(input_sha256.tsv)과 다르다.")
  ck$add("inputs_pinned", TRUE, sprintf("NSO manifest sha256 %s == instruction%s; lossless, reference_stages, input_sha256 checked against that manifest; NRV edges.rds against its stage manifest and the NSO input record",
                                        substr(hm, 1, 12), if (is.null(pins)) " (pin disabled)" else ""))
  ls <- readRDS(inp[["nso_lossless"]]); md <- ls$metadata; L <- ls$gene_evidence$cell$log; nd <- ls$programs$cells$n_detected; mk <- ls$reference$stage_markers; rm(ls)
  n <- nrow(md)
  ok_n <- n == NCA_N_CELLS && nrow(L) == n && nrow(nd) == n && !anyDuplicated(md$cell_id)
  ck$add("denominator_263159", ok_n, sprintf("%d cells in metadata, log matrix and saved detection; unique cell_id", n)); if (!ok_n) stop("세포 분모가 263,159 가 아니다: ", n)
  e <- readRDS(inp[["nrv_p92_edges"]]); if (nrow(e) && any(e$i < 1L | e$j > n | e$i >= e$j)) stop("edge 가 1..n 상삼각 세포 번호가 아니다.")
  block <- as.character(md$block); log("inputs pinned:", n, "cells,", nrow(e), "edges")
  # ── 검출 집합 (저장 검출 수와 대조) ──
  markers <- split(mk$gene, factor(mk$stage_id, levels = NSO_STAGES)); markers <- lapply(markers, as.character)
  sets <- lapply(NCA_CRITERIA, function(k) nca_detection_sets(L, markers, k))
  cons <- all(vapply(names(NCA_CRITERIA), function(cr) identical(unname(sets[[cr]]), unname(!is.na(nd) & nd >= NCA_CRITERIA[[cr]])), TRUE))
  ck$add("sets_match_saved_detection", cons, "sets from saved log columns > 0 equal the saved per-cell detected-marker counts (>= 1 and >= 2) for all 6 stages")
  if (!cons) stop("log 열로 만든 검출 집합이 저장된 검출 수와 다르다.")
  det <- data.frame(stage = NSO_STAGES, n_markers = lengths(markers), markers = vapply(markers, paste, "", collapse = "|"),
                    n_cells_ge1 = colSums(sets$ge1), n_cells_ge2 = colSums(sets$ge2), row.names = NULL, stringsAsFactors = FALSE)
  # ── 전체 그래프 성분 ──
  wg <- do.call(rbind, lapply(names(NCA_GRAPHS), function(gr) { z <- nca_components(n, e, block = block, same_block = NCA_GRAPHS[[gr]])
    data.frame(graph = gr, n_cells = z$n_cells, n_edges = z$n_edges, n_components = z$n_components, largest_component = z$largest, n_isolated = z$n_isolated, stringsAsFactors = FALSE) }))
  log("whole graph:", paste(wg$graph, wg$n_components, wg$largest_component, sep = ":", collapse = " "))
  # ── 단계 성분 · 인접 쌍 · 순서 도달 (기준 x 그래프) ──
  sc <- list(); pe <- list(); rt <- list(); fb <- list()
  for (cr in names(NCA_CRITERIA)) for (gr in names(NCA_GRAPHS)) { S <- sets[[cr]]; sb <- NCA_GRAPHS[[gr]]
    for (k in seq_along(NSO_STAGES)) { z <- nca_components(n, e, keep = S[, k], block = block, same_block = sb)
      sc[[length(sc) + 1L]] <- data.frame(criterion = cr, graph = gr, stage = NSO_STAGES[k], n_cells = z$n_cells, n_within_edges = z$n_edges, n_components = z$n_components,
                                          largest_component = z$largest, n_isolated = z$n_isolated, stringsAsFactors = FALSE) }
    for (k in seq_len(length(NSO_STAGES) - 1L)) pe[[length(pe) + 1L]] <- data.frame(criterion = cr, graph = gr, from = NSO_STAGES[k], to = NSO_STAGES[k + 1L],
                                                                                   nca_pair_edges(e, S[, k], S[, k + 1L], block, same_block = sb), stringsAsFactors = FALSE)
    r <- nca_ordered_reach(n, e, S, block = block, same_block = sb)
    rt[[length(rt) + 1L]] <- data.frame(criterion = cr, graph = gr, r$table, first_break = r$first_break, stringsAsFactors = FALSE)
    log(cr, gr, "first break:", if (nzchar(r$first_break)) r$first_break else "none", "| reachable:", paste(r$table$n_reachable_cells, collapse = ",")) }
  bind <- function(x) { d <- do.call(rbind, x); rownames(d) <- NULL; d }; sc <- bind(sc); pe <- bind(pe); rt <- bind(rt)
  ok_rt <- all(tapply(rt$n_reachable_cells <= rt$n_cells, paste(rt$criterion, rt$graph), all)) &&
    all(vapply(split(rt, paste(rt$criterion, rt$graph)), function(z) { b <- which(!z$reached); !length(b) || all(!z$reached[min(b):nrow(z)]) }, TRUE))
  ck$add("reach_tables_consistent", ok_rt, "reachable <= detected in every stage; after the first unreached stage no later stage is reached")
  # ── 표 · marker 표 · 그림 ──
  files <- character(0); rtp <- logical(0)
  wt <- function(df, rel) { p0 <- file.path(out, rel); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); rtp[[rel]] <<- mbvt_write_tsv(df, p0); files <<- c(files, p0) }
  p0 <- file.path(out, "stage_marker_table.tsv"); file.copy(inp[["nso_reference_stages"]], p0); files <- c(files, p0)
  ck$add("marker_table_copied_as_is", identical(.nrv_sha(p0), .nrv_sha(inp[["nso_reference_stages"]])), "stage_marker_table.tsv is a byte copy of NSO reference_stages.tsv")
  wt(det, "detection_sets.tsv"); wt(wg, "whole_graph_components.tsv"); wt(sc, "stage_components.tsv"); wt(pe, "adjacent_pair_edges.tsv"); wt(rt, "ordered_reachability.tsv")
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; dir.create(file.path(out, "figures"), showWarnings = FALSE)
  p0 <- file.path(out, "figures", "fig_connectivity_summary.png"); ft <- nca_draw_summary(p0, rt, ko, fam); files <- c(files, p0); wt(ft, "figure_text.tsv")
  # ── provenance ──
  p0 <- file.path(out, "lossless.rds")
  saveRDS(list(contract = "neut_stage_observation_2026-10-19/codex_connectivity_instruction.md", criteria = NCA_CRITERIA, graphs = NCA_GRAPHS, markers = markers,
               detection = det, sets = sets, whole_graph = wg, stage_components = sc, pair_edges = pe, reach = rt), p0, version = 3); files <- c(files, p0)
  wt(data.frame(input = names(inp), path = unname(inp), sha256 = vapply(inp, .nrv_sha, "", USE.NAMES = FALSE), stringsAsFactors = FALSE), "input_sha256.tsv")
  wt(snap, "source_snapshot.tsv")
  wt(data.frame(item = c("R_executable", "R", "libPaths", "packages", "title_font", "random_seed", "command", "working_dir", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, paste(.libPaths(), collapse = ";"),
                          paste(vapply(c("Matrix", "igraph", "digest", "jsonlite"), function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "),
                          if (ko) MBVR_TITLE_FAMILY else "default", "none (no random numbers are drawn)", paste(c("62_nso_connectivity_audit.R", argv), collapse = " "), getwd(),
                          format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), stringsAsFactors = FALSE), "environment.tsv")
  ck$add("tables_roundtrip_ok", all(rtp), sprintf("%d TSVs re-read with the same rows, field counts and text columns", length(rtp)))
  cks <- ck$table(); wt(cks, "checks.tsv")
  rel <- substring(normalizePath(files), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .nca_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.nca_bb, .f))
  quit(status = nca_main(commandArgs(trailingOnly = TRUE)))
}
