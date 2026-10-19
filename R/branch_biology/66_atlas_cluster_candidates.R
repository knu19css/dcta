# =============================================================================
# 66_atlas_cluster_candidates.R — atlas 6 단계 패널의 군집 평균 점수 → 군집 내부 1 위 후보(score_top_candidate) · 단계별 최고 군집 (계산 · 실행기)
# =============================================================================
#
# 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_cluster_instruction.md (2026-10-19). Ruling: 같은 폴더 claude_cluster_progress.md
#   - 입력: atlas attempt_02 (manifest hash 고정)의 lossless.rds · cluster_programs.tsv. raw RNA 를 다시 읽지 않는다. 6 marker 프로그램만 쓴다(패널 10 개씩, 기존 정규화).
#   - 주 점수 = mean_log_pooled (군집 안 모든 세포의 패널 유전자 logRNA 평균). 보조 = mean_log_equal_block (가용 block 이 없으면 결측, pooled 로 바꾸지 않는다).
#   - 군집마다 6 점수의 1 위 · 2 위 · 차이. 이름은 score_top_candidate 다 (확률 · 확정 cell type · 분화 판정이 아니다).
#     모두 정확히 0 = unassigned_zero, 모두 결측 = unavailable, 1 위 동률(양수 점수끼리 |차이| <= ATLAS_CL_TIE_TOL = R23 tie_tol) = tied (후보 모두),
#     일부 결측 = incomplete_panel_comparison. 허용오차는 동률 비교에만 쓰고 0 판정에는 쓰지 않는다 (Codex 진행 중 검토).
#   - 작은 군집(n < NSO_MIN_CELLS)은 경고만 붙이고 숨기거나 합치지 않는다. 한 군집이 여러 단계의 최고일 수 있고, 어떤 단계도 1 위가 아닐 수 있다(유일 대응 없음).
#   - 단계별 최고 군집(모든 군집 중 그 단계 점수 최고)은 군집 내부 1 위와 다른 비교다 (stage_peak_clusters.csv).
#   - 계산 표를 calc_manifest.json 으로 고정한 뒤에만 원 UMAP(세포 ID 순서 확인, 표시 전용)과 사후 원저자 맥락 표(attempt_02 posthoc 표)를 읽는다.
#   R18 · R23 · R46 · R59 · R64 를 수정 없이 재사용한다. 결과 폴더가 비어 있지 않으면 쓰지 않는다(exit 6).
#
# 실행 (W 에서, dcta 환경만)
#   run-r dcta R/branch_biology/66_atlas_cluster_candidates.R --data-root <M> --out-dir <W>/results/neutrophil_atlas_rna_2026-10-19/cluster_view_01 \
#     [--attempt-dir <W>/results/neutrophil_atlas_rna_2026-10-19/attempt_02]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

ATLAS_CL_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                      "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                      "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                      "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                      "67_atlas_cluster_view.R")
ATLAS_CL_ATTEMPT <- "results/neutrophil_atlas_rna_2026-10-19/attempt_02"                       # W 기준 (읽기 전용)
ATLAS_CL_PINNED <- list(attempt = "3296af42a68671a01f0f87bdabaf78281e7aa47b349e27f822b7c607df31e635",   # attempt_02 manifest.json (지시서)
                        h5ad = "4149b357b7d0027027988b8c78a92ae8894f36f3eae667762500dbc976f3c454", n_cells = 263159L)
ATLAS_CL_TIE_TOL <- 1e-12                                                                     # = R23 NRV_THRESHOLDS tie_tol (실행 때 대조)
ATLAS_CL_SCORE_COLS <- c("n_cells", "n_blocks", "n_blocks_min_cells", "mean_log_pooled", "mean_log_equal_block", "mean_det_rate_pooled")
ATLAS_CL_NOTE <- "score_top_candidate = highest cluster-mean panel score; not a probability, a confirmed cell type or a differentiation call; panel genes are also expressed outside the neutrophil lineage"


# ── 1. 계산 (순수) ───────────────────────────────────────────────────────────
#' 군집 x 6 marker 프로그램 점수 (긴 표). clusters = R59 nso_program_evidence()$clusters (partition · cluster · program_id · n_cells · … · status)
atlas_cl_scores <- function(clusters) {
  need <- c("partition", "cluster", "program_id", ATLAS_CL_SCORE_COLS, "status"); miss <- setdiff(need, names(clusters))
  if (length(miss)) stop("군집 프로그램 표에 필요한 열이 없다: ", paste(miss, collapse = ", "))
  d <- clusters[clusters$program_id %in% paste0("marker_", NSO_STAGES), need, drop = FALSE]
  if (!nrow(d)) stop("군집 프로그램 표에 6 단계 marker 프로그램이 없다.")
  d$stage <- sub("^marker_", "", d$program_id); d$cluster <- as.character(d$cluster)
  parts <- unique(as.character(d$partition))
  out <- do.call(rbind, lapply(parts, function(p) { z <- d[d$partition == p, , drop = FALSE]; u <- fdt_cluster_order_safe(z$cluster)
    k <- table(factor(z$cluster, levels = u)); if (any(k != length(NSO_STAGES))) stop("군집마다 6 단계 marker 프로그램이 하나씩이어야 한다: ", p)
    z[order(match(z$cluster, u), match(z$stage, NSO_STAGES)), , drop = FALSE] }))
  data.frame(partition = out$partition, cluster = out$cluster, stage = out$stage, out[ATLAS_CL_SCORE_COLS], program_status = out$status, row.names = NULL, stringsAsFactors = FALSE)
}

#' 긴 표 → 군집 x 단계 행렬 (군집 순서 = fdt_cluster_order_safe)
atlas_cl_matrix <- function(sc, partition, col) {
  z <- sc[sc$partition == partition, , drop = FALSE]; u <- fdt_cluster_order_safe(z$cluster)
  M <- matrix(NA_real_, length(u), length(NSO_STAGES), dimnames = list(u, NSO_STAGES)); M[cbind(match(z$cluster, u), match(z$stage, NSO_STAGES))] <- z[[col]]; M
}

#' 행마다 1 위 후보 · 2 위 · 차이. 동률 = 1 위와의 차이가 tol 이하 (후보는 단계 순서로 모두). 모두 0 = unassigned_zero, 모두 결측 = unavailable.
atlas_cl_rank <- function(S, tol) {
  st <- if (is.null(colnames(S))) NSO_STAGES[seq_len(ncol(S))] else colnames(S)
  rows <- lapply(seq_len(nrow(S)), function(i) { v <- S[i, ]; ok <- !is.na(v); nav <- sum(ok)
    scope <- if (!nav) "none" else if (nav < length(v)) "incomplete_panel_comparison" else "complete"
    if (!nav) return(data.frame(score_top_status = "unavailable", score_top_candidate = "", score_top_value = NA_real_, score_second_candidate = "", score_second_value = NA_real_,
                                score_top_minus_second = NA_real_, n_top_tied = 0L, comparison_scope = scope, n_stages_available = 0L, stringsAsFactors = FALSE))
    o <- order(-v[ok], match(st[ok], st)); vs <- unname(v[ok][o]); ss <- st[ok][o]; top <- vs[1]
    sv <- if (nav > 1L) vs[2] else NA_real_; sc <- if (nav > 1L) ss[2] else ""
    if (top == 0) return(data.frame(score_top_status = "unassigned_zero", score_top_candidate = "", score_top_value = top, score_second_candidate = "", score_second_value = sv,
                                    score_top_minus_second = top - sv, n_top_tied = 0L, comparison_scope = scope, n_stages_available = nav, stringsAsFactors = FALSE))
    tied <- st[st %in% ss[vs > 0 & vs >= top - tol]]                                           # 0(미검출)은 작은 양수와 동률 후보로 묶지 않는다 (Codex 검토)
    data.frame(score_top_status = if (length(tied) > 1L) "tied" else "score_top_candidate", score_top_candidate = paste(tied, collapse = "|"), score_top_value = top,
               score_second_candidate = sc, score_second_value = sv, score_top_minus_second = top - sv, n_top_tied = length(tied), comparison_scope = scope,
               n_stages_available = nav, stringsAsFactors = FALSE) })
  d <- do.call(rbind, rows); d <- data.frame(cluster = if (is.null(rownames(S))) as.character(seq_len(nrow(S))) else rownames(S), d, stringsAsFactors = FALSE); rownames(d) <- NULL; d
}

#' 단계별 최고 군집: 모든 군집 중 그 단계 점수가 가장 높은 군집 (동률 모두). 군집 내부 1 위와 다른 비교 — 최고 군집 자신의 1 위 후보를 같이 적는다.
atlas_cl_stage_peaks <- function(S, rk, tol, n_cells = NULL) {
  st <- colnames(S); cl <- rownames(S); tops <- strsplit(rk$score_top_candidate[match(cl, rk$cluster)], "|", fixed = TRUE)
  nc <- if (is.null(n_cells)) rep(NA_integer_, length(cl)) else as.integer(n_cells[cl])
  rows <- lapply(st, function(s) { v <- S[, s]; ok <- !is.na(v); ntop <- sum(vapply(tops, function(z) s %in% z, TRUE))
    base <- data.frame(stage = s, peak_status = "", peak_clusters = "", peak_score = NA_real_, n_tied = 0L, peak_n_cells = "", second_cluster = "", second_score = NA_real_,
                       peak_minus_second = NA_real_, peak_cluster_top_candidate = "", peak_is_cluster_top = NA, n_clusters_with_stage_top = ntop, stringsAsFactors = FALSE)
    if (!any(ok)) { base$peak_status <- "unavailable"; return(base) }
    o <- order(-v[ok], match(cl[ok], cl)); vs <- unname(v[ok][o]); cs <- cl[ok][o]; top <- vs[1]
    base$peak_score <- top; if (length(vs) > 1L) { base$second_score <- vs[2]; base$second_cluster <- cs[2]; base$peak_minus_second <- top - vs[2] }
    if (top == 0) { base$peak_status <- "all_zero"; return(base) }
    pc <- cl[cl %in% cs[vs > 0 & vs >= top - tol]]; j <- match(pc, cl)
    base$peak_status <- if (length(pc) > 1L) "tied" else "peak"; base$peak_clusters <- paste(pc, collapse = "|"); base$n_tied <- length(pc)
    base$peak_n_cells <- if (is.null(n_cells)) "" else paste(nc[j], collapse = "|")
    base$peak_cluster_top_candidate <- paste(vapply(tops[j], paste, "", collapse = "|"), collapse = ";"); base$peak_is_cluster_top <- any(vapply(tops[j], function(z) s %in% z, TRUE))
    base })
  d <- do.call(rbind, rows); rownames(d) <- NULL; d
}

#' 단계별 검출 세포 비율 (패널 유전자 >= 1, >= 2 검출), 군집마다. 단계 프로그램이 판정불가면 결측
atlas_cl_detection <- function(md, nd) {
  do.call(rbind, lapply(c("p92", "p51"), function(p) { k <- if (p == "p92") md$cluster92 else md$cluster51; u <- fdt_cluster_order_safe(k); f <- factor(k, levels = u)
    do.call(rbind, lapply(NSO_STAGES, function(s) { v <- nd[, paste0("marker_", s)]; av <- !all(is.na(v))
      data.frame(partition = p, cluster = u, stage = s, frac_cells_ge1 = if (av) as.numeric(tapply(v >= 1, f, mean)) else NA_real_,
                 frac_cells_ge2 = if (av) as.numeric(tapply(v >= 2, f, mean)) else NA_real_, stringsAsFactors = FALSE) })) }))
}

#' 세포 단위 패널 평균(R59 cells$mean_log)의 군집 평균 = mean_log_pooled 의 독립 재계산 (선형이라 같아야 한다)
atlas_cl_recompute_pooled <- function(md, cm) {
  do.call(rbind, lapply(c("p92", "p51"), function(p) { k <- if (p == "p92") md$cluster92 else md$cluster51; u <- fdt_cluster_order_safe(k); f <- factor(k, levels = u)
    do.call(rbind, lapply(NSO_STAGES, function(s) { v <- cm[, paste0("marker_", s)]
      data.frame(partition = p, cluster = u, stage = s, mean_cells = if (all(is.na(v))) NA_real_ else as.numeric(tapply(v, f, mean)), stringsAsFactors = FALSE) })) }))
}

#' 군집 대표 표시 위치: 군집 세포 중 원 UMAP 좌표 중앙값에 가장 가까운 실제 세포 (표시 전용, 계산에 쓰지 않는다)
atlas_cl_rep_cells <- function(U, cluster) {
  cluster <- as.character(cluster); u <- fdt_cluster_order_safe(cluster)
  do.call(rbind, lapply(u, function(k) { i <- which(cluster == k); mx <- stats::median(U[i, 1]); my <- stats::median(U[i, 2])
    j <- i[which.min((U[i, 1] - mx)^2 + (U[i, 2] - my)^2)]
    data.frame(cluster = k, cell = j, x = U[j, 1], y = U[j, 2], n_cells = length(i), stringsAsFactors = FALSE) }))
}

#' 군집 후보 표 (넓은 표): pooled 주 점수 · equal-block 보조 점수 · 1 위 / 2 위 · 작은 군집 경고 · 두 가중의 1 위 일치 여부
atlas_cl_candidates <- function(sc, tol, min_cells) {
  do.call(rbind, lapply(unique(sc$partition), function(p) { S <- atlas_cl_matrix(sc, p, "mean_log_pooled"); E <- atlas_cl_matrix(sc, p, "mean_log_equal_block")
    z <- sc[sc$partition == p & sc$stage == NSO_STAGES[1], , drop = FALSE]; n <- stats::setNames(z$n_cells, z$cluster)[rownames(S)]
    rp <- atlas_cl_rank(S, tol); re <- atlas_cl_rank(E, tol)
    d <- data.frame(partition = p, cluster = rownames(S), n_cells = as.integer(n), small_cluster = as.integer(n) < min_cells, stringsAsFactors = FALSE)
    d$warning <- ifelse(d$small_cluster, sprintf("n_cells %d < NSO_MIN_CELLS %d: kept, not merged or hidden", d$n_cells, as.integer(min_cells)), "")
    for (s in NSO_STAGES) d[[paste0("score_", s)]] <- unname(S[, s])
    d <- cbind(d, rp[, setdiff(names(rp), "cluster")])
    for (s in NSO_STAGES) d[[paste0("eb_", s)]] <- unname(E[, s])
    d$eb_top_status <- re$score_top_status; d$eb_top_candidate <- re$score_top_candidate; d$eb_top_value <- re$score_top_value
    d$eb_top_minus_second <- re$score_top_minus_second; d$eb_comparison_scope <- re$comparison_scope
    pos <- c("score_top_candidate", "tied")
    d$equal_block_agreement <- ifelse(d$eb_top_status == "unavailable", "equal_block_unavailable",
                                      ifelse(d$score_top_status %in% pos & d$eb_top_status %in% pos, ifelse(d$score_top_candidate == d$eb_top_candidate, "same_top", "different_top"),
                                             ifelse(d$score_top_status == "unassigned_zero" & d$eb_top_status == "unassigned_zero", "both_zero", "not_comparable")))
    d$note <- rep(ATLAS_CL_NOTE, nrow(d)); d }))
}


# ── 2. 실행기 ────────────────────────────────────────────────────────────────
atlas_cl_main <- function(argv, root = NULL, reader = NULL, pins = ATLAS_CL_PINNED) {
  a <- list(attempt = NULL, data_root = NULL, out = NULL); keys <- c("--attempt-dir" = "attempt", "--data-root" = "data_root", "--out-dir" = "out"); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (any(vapply(a[c("data_root", "out")], function(x) is.null(x) || is.na(x) || !nzchar(x), TRUE))) { cat("--data-root, --out-dir 가 필요하다 (--attempt-dir 는 선택).\n"); return(2L) }
  out <- a$out
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  if (is.null(a$attempt) || is.na(a$attempt)) a$attempt <- file.path(root, ATLAS_CL_ATTEMPT)
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", ATLAS_CL_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(ATLAS_CL_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = ATLAS_CL_SOURCES)
    get("atlas_cl_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, reader = reader, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

atlas_cl_build <- function(a, out, snap, log, argv, reader, pins) {
  ck <- .mbvq_checks(); t_start <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  # ── 입력 고정 ──
  mf <- file.path(a$attempt, "manifest.json"); if (!file.exists(mf)) stop("attempt manifest 가 없다: ", mf)
  hm <- .nrv_sha(mf); if (!is.null(pins$attempt) && !identical(hm, pins$attempt)) stop("attempt manifest hash 가 지시서 고정값과 다르다: ", hm)
  fk <- c(lossless = "lossless.rds", cluster_programs = "cluster_programs.tsv", cluster_comparison = "cluster_comparison.tsv",
          posthoc_composition = "posthoc_author_cluster_composition.tsv", input_sha256 = "input_sha256.tsv")
  inp <- vapply(fk, function(f) .mbvq_manifest_file(a$attempt, "manifest.json", f), "")
  if (is.null(reader)) reader <- nrv_h5_reader(a$data_root)
  h5 <- reader$input_files()[["h5ad"]]; if (is.na(h5) || !file.exists(h5)) stop("h5ad 가 없다: ", h5)
  hh <- .nrv_sha(h5); rec <- utils::read.delim(inp[["input_sha256"]], stringsAsFactors = FALSE, colClasses = "character")
  if (!identical(hh, rec$sha256[match("h5ad", rec$input)])) stop("h5ad hash 가 attempt 입력 기록(input_sha256.tsv)과 다르다: ", hh)
  if (!is.null(pins$h5ad) && !identical(hh, pins$h5ad)) stop("h5ad hash 가 계약값과 다르다: ", hh)
  ck$add("inputs_pinned", TRUE, sprintf("attempt manifest %s == instruction%s; lossless, cluster_programs, cluster_comparison, posthoc composition, input_sha256 == that manifest; h5ad == attempt record%s",
                                        substr(hm, 1, 12), if (is.null(pins$attempt)) " (pin disabled)" else "", if (is.null(pins$h5ad)) "" else " and contract"))
  # ── 세포 · 군집 점수 (raw 를 다시 읽지 않는다) ──
  ls <- readRDS(inp[["lossless"]]); md <- ls$metadata; pe <- ls$programs; rm(ls); n <- nrow(md)
  okid <- !anyDuplicated(md$cell_id) && !anyNA(md$cluster92) && !anyNA(md$cluster51) && nrow(pe$cells$mean_log) == n && nrow(pe$cells$n_detected) == n
  ck$add("cell_ids_unique_complete", okid, sprintf("%d unique cell_id; p92 / p51 cluster for every cell; cell-level program matrices have %d rows", n, n))
  if (!okid) stop("세포 ID · 군집 · 세포 행렬이 맞지 않는다.")
  if (!is.null(pins$n_cells)) { ok <- n == pins$n_cells; ck$add("denominator_263159", ok, sprintf("%d cells == pinned %d", n, pins$n_cells)); if (!ok) stop("세포 분모가 고정값과 다르다: ", n) }
  sc <- atlas_cl_scores(pe$clusters)
  oks <- all(vapply(c("p92", "p51"), function(p) { k <- if (p == "p92") md$cluster92 else md$cluster51; z <- sc[sc$partition == p & sc$stage == NSO_STAGES[1], ]
    tb <- table(k); setequal(names(tb), z$cluster) && all(as.integer(tb[z$cluster]) == z$n_cells) }, TRUE))
  ck$add("cluster_sizes_match_cells", oks, "every p92 / p51 cluster size in the program table equals its cell count in metadata (no cluster hidden, merged or added)")
  if (!oks) stop("군집 크기가 세포표와 다르다.")
  tsv <- utils::read.delim(inp[["cluster_programs"]], stringsAsFactors = FALSE, quote = ""); tsv <- tsv[startsWith(tsv$program_id, "marker_"), , drop = FALSE]
  j <- match(paste(sc$partition, sc$cluster, sc$stage), paste(tsv$partition, as.character(tsv$cluster), sub("^marker_", "", tsv$program_id)))
  same <- function(x, y) identical(is.na(x), is.na(y)) && all(abs(x[!is.na(x)] - y[!is.na(y)]) <= 1e-12 * pmax(1, abs(y[!is.na(y)])))
  okt <- !anyNA(j) && nrow(tsv) == nrow(sc) && all(vapply(ATLAS_CL_SCORE_COLS, function(k) same(as.numeric(sc[[k]]), as.numeric(tsv[[k]][j])), TRUE)) &&
    identical(sc$program_status, tsv$status[j])
  ck$add("cluster_programs_tsv_matches_lossless", okt, sprintf("%d cluster x stage marker rows: n_cells, n_blocks, n_blocks_min_cells, pooled, equal-block, det rate and status equal attempt cluster_programs.tsv", nrow(sc)))
  if (!okt) stop("lossless 군집 점수가 cluster_programs.tsv 와 다르다.")
  rc <- atlas_cl_recompute_pooled(md, pe$cells$mean_log); m <- match(paste(sc$partition, sc$cluster, sc$stage), paste(rc$partition, rc$cluster, rc$stage))
  dif <- abs(rc$mean_cells[m] - sc$mean_log_pooled); okr <- !anyNA(m) && identical(is.na(rc$mean_cells[m]), is.na(sc$mean_log_pooled)) && max(c(0, dif), na.rm = TRUE) <= 1e-9
  ck$add("pooled_scores_recomputed_from_cells", okr, sprintf("cluster means of the per-cell panel means equal mean_log_pooled (max |diff| %.3g)", max(c(0, dif), na.rm = TRUE)))
  if (!okr) stop("세포 평균에서 다시 계산한 군집 점수가 mean_log_pooled 와 다르다.")
  det <- atlas_cl_detection(md, pe$cells$n_detected)
  cc <- utils::read.delim(inp[["cluster_comparison"]], stringsAsFactors = FALSE, quote = ""); jd <- match(paste(det$partition, det$cluster, det$stage), paste(cc$partition, as.character(cc$cluster), cc$stage))
  okd <- !anyNA(jd) && same(det$frac_cells_ge1, as.numeric(cc$atlas_frac_ge1[jd])) && same(det$frac_cells_ge2, as.numeric(cc$atlas_frac_ge2[jd]))
  ck$add("detection_fractions_match_attempt", okd, "per cluster x stage fractions of cells with >= 1 / >= 2 panel genes detected equal attempt cluster_comparison.tsv (R64 connectivity sets)")
  if (!okd) stop("검출 세포 비율이 attempt cluster_comparison.tsv 와 다르다.")
  tt <- as.numeric(NRV_THRESHOLDS$value[NRV_THRESHOLDS$name == "tie_tol"]); ck$add("tie_tol_matches_r23", identical(tt, ATLAS_CL_TIE_TOL), sprintf("tie tolerance %g == R23 tie_tol", ATLAS_CL_TIE_TOL))
  # ── 후보 · 단계별 최고 군집 ──
  cand <- atlas_cl_candidates(sc, ATLAS_CL_TIE_TOL, NSO_MIN_CELLS)
  peaks <- do.call(rbind, lapply(c("p92", "p51"), function(p) { S <- atlas_cl_matrix(sc, p, "mean_log_pooled"); E <- atlas_cl_matrix(sc, p, "mean_log_equal_block")
    z <- cand[cand$partition == p, ]; nn <- stats::setNames(z$n_cells, z$cluster)
    pk <- atlas_cl_stage_peaks(S, data.frame(cluster = z$cluster, score_top_candidate = z$score_top_candidate, stringsAsFactors = FALSE), ATLAS_CL_TIE_TOL, nn)
    pe2 <- atlas_cl_stage_peaks(E, data.frame(cluster = z$cluster, score_top_candidate = z$eb_top_candidate, stringsAsFactors = FALSE), ATLAS_CL_TIE_TOL, nn)
    data.frame(partition = p, pk, eb_peak_status = pe2$peak_status, eb_peak_clusters = pe2$peak_clusters, eb_peak_score = pe2$peak_score,
               comparison_note = rep("stage peak = highest cluster for this stage among all clusters; differs from the within-cluster top (cluster_candidates.csv)", nrow(pk)),
               stringsAsFactors = FALSE) }))
  scores <- data.frame(sc[, c("partition", "cluster", "stage", ATLAS_CL_SCORE_COLS, "program_status")], det[match(paste(sc$partition, sc$cluster, sc$stage), paste(det$partition, det$cluster, det$stage)),
                                                                                                          c("frac_cells_ge1", "frac_cells_ge2")], row.names = NULL, stringsAsFactors = FALSE)
  for (p in c("p92", "p51")) { z <- cand[cand$partition == p, ]
    log(p, "clusters:", nrow(z), "|", paste(names(table(z$score_top_status)), table(z$score_top_status), sep = "=", collapse = " "), "| equal-block different_top:", sum(z$equal_block_agreement == "different_top")) }
  # ── 계산 표 → calc_manifest.json ──
  files <- character(0); crt <- logical(0)
  wc <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  wc(scores, "cluster_scores.csv"); wc(cand, "cluster_candidates.csv"); wc(peaks, "stage_peak_clusters.csv")
  p0 <- file.path(out, "lossless.rds")
  saveRDS(list(contract = "neutrophil_atlas_rna_2026-10-19/codex_cluster_instruction.md", attempt_manifest_sha256 = hm, tie_tol = ATLAS_CL_TIE_TOL, min_cells = NSO_MIN_CELLS,
               scores = scores, candidates = cand, stage_peaks = peaks, metadata = md[, c("cell_id", "cluster92", "cluster51", "block")]), p0, version = 3)
  files <- c(files, p0); cmf <- file.path(out, "calc_manifest.json")
  writeLines(jsonlite::toJSON(list(phase = "cluster scores, candidates and stage peaks fixed before the original UMAP and post hoc author context are read", written_at = stamp(),
                                   attempt_manifest_sha256 = hm, files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  calc_files <- files; log("calculation fixed:", length(files), "files")
  # ── 표시 (고정 뒤): 원 UMAP, 세포 ID 순서 확인 ──
  U <- reader$umap(); uid <- as.character(reader$obs("cell_id"))
  if (!identical(uid, md$cell_id)) stop("UMAP(h5 obs) 의 세포 ID 순서가 계산 세포표와 다르다 (재정렬하지 않는다).")
  if (!is.matrix(U) || nrow(U) != n || ncol(U) < 2L || any(!is.finite(U[, 1:2]))) stop("UMAP 모양이 세포 표와 다르다.")
  ck$add("umap_ids_identical_order", TRUE, sprintf("h5 obs cell_id == computed cell order for all %d cells; UMAP used for display only", n)); U <- U[, 1:2]
  reps <- lapply(c(p92 = "p92", p51 = "p51"), function(p) atlas_cl_rep_cells(U, if (p == "p92") md$cluster92 else md$cluster51))
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  fa <- file.path(fig, "figA_cluster_top_candidate_umap.png"); ft <- list(atlas_cl_draw_top(fa, U, md, cand, reps, ko, fam)); files <- c(files, fa)
  for (p in c("p92", "p51")) { fb <- file.path(fig, sprintf("figB_stage_cluster_scores_%s.png", p))
    ft[[length(ft) + 1L]] <- atlas_cl_draw_stage(fb, U, if (p == "p92") md$cluster92 else md$cluster51, atlas_cl_matrix(sc, p, "mean_log_pooled"), peaks[peaks$partition == p, ],
                                                 reps[[p]], cand[cand$partition == p, ], p, ko, fam); files <- c(files, fb) }
  wc(do.call(rbind, ft), "figure_text.csv"); log("figures written"); rm(U)
  # ── 사후 (고정 뒤): attempt_02 의 원저자 구성표로 군집 맥락만 ──
  comp <- utils::read.delim(inp[["posthoc_composition"]], stringsAsFactors = FALSE, quote = "", colClasses = c(cluster = "character"))
  px <- do.call(rbind, lapply(seq_len(nrow(cand)), function(i) { z <- comp[comp$partition == cand$partition[i] & comp$cluster == cand$cluster[i], , drop = FALSE]
    z <- z[order(-z$n_cells, z$author_label, method = "radix"), , drop = FALSE]
    data.frame(partition = cand$partition[i], cluster = cand$cluster[i], n_cells = cand$n_cells[i], score_top_status = cand$score_top_status[i], score_top_candidate = cand$score_top_candidate[i],
               top_author_label = if (nrow(z)) z$author_label[1] else "", top_author_frac = if (nrow(z)) z$frac_of_cluster[1] else NA_real_,
               second_author_label = if (nrow(z) > 1L) z$author_label[2] else "", second_author_frac = if (nrow(z) > 1L) z$frac_of_cluster[2] else NA_real_, stringsAsFactors = FALSE) }))
  px$note <- rep("post hoc context from attempt_02 author composition (read after calc_manifest.json); not used for any score, candidate or selection; not accuracy", nrow(px))
  wc(px, "posthoc_cluster_author_context.csv")
  calc_now <- vapply(calc_files, .nrv_sha, ""); cm0 <- jsonlite::fromJSON(cmf)$files
  ck$add("calc_fixed_before_umap_and_posthoc", identical(unname(calc_now), unname(unlist(cm0[rel(calc_files)]))),
         "UMAP, h5 obs IDs and post hoc author composition were opened only after calc_manifest.json; fixed calculation files unchanged afterwards")
  # ── provenance ──
  allinp <- c(attempt_manifest = mf, stats::setNames(inp, paste0("attempt_", names(inp))), h5ad = h5)
  wc(data.frame(input = names(allinp), path = unname(allinp), sha256 = vapply(seq_along(allinp), function(k) if (allinp[[k]] == h5) hh else .nrv_sha(allinp[[k]]), ""),
                stringsAsFactors = FALSE), "input_sha256.csv")
  wc(snap, "source_snapshot.csv")
  pk <- c("Matrix", "hdf5r", "igraph", "digest", "jsonlite")
  wc(data.frame(item = c("R_executable", "R", "platform", "libPaths", "packages", "title_font", "tie_tol", "min_cells", "random_seed", "command", "working_dir", "started_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pk, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), if (ko) MBVR_TITLE_FAMILY else "default",
                          format(ATLAS_CL_TIE_TOL), as.character(NSO_MIN_CELLS), "none (no random numbers are drawn)", paste(c("66_atlas_cluster_candidates.R", argv), collapse = " "),
                          getwd(), format(t_start, "%Y-%m-%dT%H:%M:%S%z"), stamp()), stringsAsFactors = FALSE), "environment.csv")
  ck$add("tables_roundtrip_ok", all(crt), sprintf("%d CSVs re-read with identical dimensions, names and cell text", length(crt)))
  cks <- ck$table(); wc(cks, "checks.csv")
  smry <- function(p) { z <- cand[cand$partition == p, ]; t1 <- table(factor(unlist(strsplit(z$score_top_candidate[z$score_top_status == "score_top_candidate"], "|", fixed = TRUE)), levels = NSO_STAGES))
    sprintf("%s 군집 %d: 1 위 %s · 동률 %d · 0 %d · 판정불가 %d · 작은 군집 %d · equal-block 1 위 다름 %d", p, nrow(z), paste(names(t1), as.integer(t1), sep = " ", collapse = " · "),
            sum(z$score_top_status == "tied"), sum(z$score_top_status == "unassigned_zero"), sum(z$score_top_status == "unavailable"), sum(z$small_cluster),
            sum(z$equal_block_agreement == "different_top")) }
  pkt <- function(p) { z <- peaks[peaks$partition == p, ]; paste(sprintf("%s %s", z$stage, ifelse(nzchar(z$peak_clusters), z$peak_clusters, z$peak_status)), collapse = " · ") }
  rd <- c(sprintf("# atlas 단계 패널 — 군집 단위 1 위 후보 (%s)", basename(out)), "",
          "- 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_cluster_instruction.md. 입력: attempt_02 (manifest 고정), raw 재계산 없음.",
          "- 점수 = 군집 안 모든 세포의 패널 유전자 logRNA 평균(mean_log_pooled). 보조 = 동일 block 평균(가용 block 없으면 결측).",
          sprintf("- %s.", smry("p92")), sprintf("- %s.", smry("p51")),
          paste0("- 단계별 최고 군집 (p92): ", pkt("p92"), "."), paste0("- 단계별 최고 군집 (p51): ", pkt("p51"), "."),
          "- 군집 내부 1 위(cluster_candidates.csv · 그림 A)와 단계별 최고 군집(stage_peak_clusters.csv · 그림 B)은 서로 다른 비교다.",
          "- score_top_candidate 는 확률 · 확정 cell type · 분화 판정이 아니다. 패널 유전자는 다른 계열(단핵구 · T · 적혈구계 · 전구세포)에서도 발현된다.",
          "- 그림: figures/figA_cluster_top_candidate_umap.png (p92 · p51), figures/figB_stage_cluster_scores_p92.png · _p51.png (6 단계, 같은 척도).",
          "- 원저자 라벨은 posthoc_cluster_author_context.csv 의 사후 맥락뿐이다 (계산 · 선택에 쓰지 않음).")
  p0 <- file.path(out, "README_ko.md"); writeLines(rd, p0); files <- c(files, p0)
  allf <- c(files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), attempt_manifest_sha256 = hm, files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .acl_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.acl_bb, .f))
  quit(status = atlas_cl_main(commandArgs(trailingOnly = TRUE)))
}
