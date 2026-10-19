# =============================================================================
# 46_run_mbv_mono_quant.R — monocyte 정량 발현 재가중치: 실행기
# =============================================================================
#
# 계약: results/claude_bridge/myeloid_branch_validation_2026-10-13/codex_mono_quant_instruction.md (2026-10-15, 사용자 승인 · Codex 설계 · 최종 점수 계약)
# 같은 graph · 군집 · 11 key · 학습 / 보류 block 과, key 마다 저장 monocyte 주 경로의 실제 S, T 를 그대로 둔다. c1 baseline 쌍 표에서 방향별 선 비용만 바꾼다.
#   c_expr(A→B) = c_graph(A,B) · exp(−λ · s_used(A→B))   (c_graph = −log q 그대로, 새 선 · 음수 비용 없음)
#   변형: lambda0 (저장 경로 회귀 확인) · lambda0.5 · lambda1 (주) · lambda2 · control_reverse (λ1, −s) · irf8_holdout (λ1, s = s_M02)
# r ≥ 1 필터 · 오르막 끝점 · 중간 경유점은 섞지 않는다. Yen 차선 10 개가 아니라 그 fold 의 모든 관측 선에서 찾는다.
# 입력은 hash 확인 뒤 읽기만 한다
#   MBV freeze · stream_sets · select · author (각 단계 manifest), NRV pairs · cells · edges (MBV freeze 에 고정한 hash), NRV author 표 · audit_09 세포표 (그 manifest).
# 원 RNA 는 읽지 않는다 (저장 같은 block 경계 집합 합). 원저자 라벨은 경로 계산이 모두 끝난 뒤 사후 대조 · 그림 색에만 읽는다.
# 실행: run-r dcta R/branch_biology/46_run_mbv_mono_quant.R --attempt-dir <MBV attempt> [--out-name compare_03_mono_quant]
#        [--stage-audit audit_09_fine_stage_overlay]
# =============================================================================

MBVQ_SOURCES <- c("05_cell_trajectory_io.R", "10_expression_edge_graph.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                  "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "28_view_nrv.R",
                  "29_audit_nrv.R", "30_mbv_reference.R", "31_mbv_branch.R", "32_run_mbv.R", "33_mbv_graph_view.R", "34_mbv_umap_overlay.R",
                  "35_mbv_fine_stage_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "45_mbv_mono_quant.R",
                  "46_run_mbv_mono_quant.R", "47_mbv_mono_quant_view.R")
MBVQ_ARM <- "monocyte"
MBVQ_EXPECTED_PROGRAMS <- data.frame(program = c("M01", "M02", "M02", "M02"), gene = c("IRF8", "CD14", "CCR2", "FCGR1A"), sign = c(1, 1, 1, 1),
                                     evidence = "conditional", stringsAsFactors = FALSE)       # 계약 3: 참고표 원부호 · conditional 그대로
MBVQ_DIAG_PAIRS <- list(p92 = c("57", "86"))                                     # 계약: 현재 57–86 경계 (p92 군집 번호). 선택되지 않아도 진단 행
MBVQ_GMP_FAMILY <- c("Early GMP", "GMP-Cycle", "GMP-Neut", "GMP-Mono")          # 사후 서술: GMP 계열 원저자 라벨 (성공 판정에 쓰지 않는다)
MBVQ_LEAK_VALUE <- 50                                                            # 보류 block 누설 확인: 보류 칸을 이 값으로 바꿔도 점수 · 비용이 같아야 한다

mbvq_main <- function(argv, root = NULL) {
  a <- list(attempt = NULL, name = "compare_03_mono_quant", stage = "audit_09_fine_stage_overlay"); i <- 1L
  keys <- c("--attempt-dir" = "attempt", "--out-name" = "name", "--stage-audit" = "stage")
  while (i <= length(argv)) {
    k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (k %in% names(keys)) a[[keys[[k]]]] <- v else { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    i <- i + 2L
  }
  if (is.null(a$attempt) || is.na(a$attempt) || !dir.exists(a$attempt) || is.na(a$name) || !grepl("^compare_[A-Za-z0-9_]+$", a$name) ||
      is.na(a$stage) || !grepl("^audit_[A-Za-z0-9_]+$", a$stage)) {
    cat("--attempt-dir(있는 MBV attempt), compare_ 로 시작하는 --out-name, audit_ 로 시작하는 --stage-audit 가 필요하다.\n"); return(2L) }
  out <- file.path(a$attempt, a$name)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), a$name, ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", MBVQ_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    ld <- mdu_load_snapshot(snap, files = MBVQ_SOURCES)
    get("mbvq_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}


# ── 1. 입력 (모두 hash 확인) ─────────────────────────────────────────────────
.mbvq_manifest_file <- function(dir, manifest, f) {
  mf <- file.path(dir, manifest); if (!file.exists(mf)) stop("manifest 가 없다: ", mf)
  h <- jsonlite::fromJSON(mf)$files[[f]]; p <- file.path(dir, f)
  if (is.null(h) || !file.exists(p) || !identical(.nrv_sha(p), h)) stop("파일 hash 가 그 manifest 와 다르다: ", basename(dir), "/", f)
  p
}
.mbvq_stage_file <- function(attempt, stage, f) .mbvq_manifest_file(nrv_stage_dir(attempt, stage), "stage_manifest.json", f)

.mbvq_inputs <- function(a) {
  attempt <- a$attempt; mbv_check_upstream(attempt, "view")
  inp <- c(plan_freeze = .mbvq_stage_file(attempt, "freeze", "plan_freeze.json"))
  fz <- jsonlite::fromJSON(inp[["plan_freeze"]]); m <- list(attempt = attempt, nrv_attempt = fz$args$nrv_attempt)
  inp <- c(inp, freeze = .mbvq_stage_file(attempt, "freeze", "freeze_lossless.rds"), set_sums = .mbvq_stage_file(attempt, "stream_sets", "set_sums.rds"),
           selection = .mbvq_stage_file(attempt, "select", "selection.rds"), pairs = .mbv_nrv_path(m, "structure", "pairs.rds"),
           cells = .mbv_nrv_path(m, "structure", "cells_lossless.rds"), edges = .mbv_nrv_path(m, "structure", "edges.rds"))
  cl <- readRDS(inp[["cells"]])
  list(attempt = attempt, m = m, inp = inp, partition = fz$args$partition, fr = readRDS(inp[["freeze"]]), ss = readRDS(inp[["set_sums"]]),
       sel = readRDS(inp[["selection"]]), pr = readRDS(inp[["pairs"]]), cell_cluster = as.character(cl$cells$cluster), cell_id = as.character(cl$cells$cell_id),
       blocks = as.character(cl$bm$block))
}

.mbvq_checks <- function() {
  rows <- list()
  list(add = function(name, pass, detail) { rows[[length(rows) + 1L]] <<- data.frame(check = name, pass = isTRUE(pass), detail = detail, stringsAsFactors = FALSE); invisible(isTRUE(pass)) },
       table = function() do.call(rbind, rows))
}


# ── 2. 선 보기: 선택 경로 · 진단 선의 유전자별 학습 변화 · 점수 · 비용 분해 (한 행 = 한 방향 선) ─────────
.mbvq_edge_view <- function(ed, k, ge, genes, s_used, cost, score_status, on_path) {
  out <- data.frame(edge = ed$edge[k], from = ed$from[k], to = ed$to[k], on_selected_path = on_path, c_graph = ed$c_graph[k], q = ed$q[k],
                    n_cell_edges = ed$n_cell_edges[k], n_same_block_edges = ed$n_same_block_edges[k], n_valid_train_blocks = ed$n_valid_train_blocks[k],
                    train_boundary_cells = ed$train_boundary_cells[k], s_M01 = ed$s_M01[k], n_M01 = ed$n_M01[k], s_M02 = ed$s_M02[k], n_M02 = ed$n_M02[k],
                    s_main = ed$s_main[k], status_main = ed$status_main[k], s_m02only = ed$s_m02only[k], status_m02only = ed$status_m02only[k],
                    score_status = score_status[k], s_used = s_used[k], c_expr = cost[k], cost_change = cost[k] - ed$c_graph[k], stringsAsFactors = FALSE)
  for (g in genes) { z <- ge[ge$gene == g, , drop = FALSE]; j <- match(out$edge, z$edge)
    out[[paste0("d_", g)]] <- z$d[j]; out[[paste0("detrate_d_", g)]] <- z$detrate_d[j]; out[[paste0("n_valid_blocks_", g)]] <- z$n_valid_blocks[j]
    out[[paste0("det_sum_", g)]] <- z$det_sum[j]; out[[paste0("reason_", g)]] <- z$reason[j] }
  out
}


# ── 3. key x 변형: 점수 · 비용 · 고정 양끝 방향 경로 · 보류 block 사후 측정 ─────────
.mbvq_keys <- function(x, ba, PG, anc, ck, log) {
  mind <- .nrv_thr("gene_min_detected_cells"); minc <- .nrv_thr("set_min_cells"); V <- MBVQ_VARIANTS; genes <- unique(PG$gene)
  dg <- MBVQ_DIAG_PAIRS[[x$partition]]; dg_edges <- if (length(dg) == 2L) c(paste(dg[1], dg[2], sep = ">"), paste(dg[2], dg[1], sep = ">")) else character(0)
  ok <- list(lam0 = TRUE, fix = TRUE, neg = TRUE, newe = TRUE, leak = TRUE, rev = TRUE, e3 = TRUE); n_e3 <- 0L; n_fix <- 0L
  RC <- list(); RE <- list(); DE <- list(); GE <- list(); HG <- list(); TR <- list(); keep <- list(); RE0 <- NULL   # RE0 = 행이 없을 때 쓰는 0 행 틀
  for (en in x$sel$entries) {
    kc <- data.frame(partition = x$partition, key = en$key, holdout_scheme = en$holdout_scheme, fold = as.character(en$fold), stringsAsFactors = FALSE)
    pt <- x$pr$folds[[en$key]]; if (is.null(pt)) stop("fold 쌍 표가 없다: ", en$key)
    two <- function(v) c(v, v)
    ed <- data.frame(from = as.character(c(pt$cluster_a, pt$cluster_b)), to = as.character(c(pt$cluster_b, pt$cluster_a)), c_graph = two(pt$cost), q = two(pt$q),
                     n_cell_edges = two(pt$n_cell_edges), n_same_block_edges = two(pt$n_same_block_edges), stringsAsFactors = FALSE)
    ed$edge <- paste(ed$from, ed$to, sep = ">"); tr <- as.character(en$train); ho <- as.character(en$held)
    bk <- mbvq_ba_subset(ba, ed$edge); vt <- bk$reason[, tr, drop = FALSE] == "ok"
    ed$n_valid_train_blocks <- as.integer(rowSums(vt))
    ed$train_boundary_cells <- as.numeric(rowSums((bk$n_from[, tr, drop = FALSE] + bk$n_to[, tr, drop = FALSE]) * vt))
    ge <- mbvq_gene_effects(bk, tr, mind); ps <- mbvq_program_scores(ge, PG); vs <- mbvq_variant_scores(ps)
    ed <- cbind(ed, ps[match(ed$edge, ps$edge), setdiff(names(ps), "edge")]); rownames(ed) <- NULL
    if (is.null(RE0)) RE0 <- cbind(kc[0, , drop = FALSE], variant = character(0), lambda = numeric(0), row_role = character(0), step = integer(0),
                                   diagnostic_edge = logical(0), .mbvq_edge_view(ed, integer(0), ge, genes, numeric(0), numeric(0), character(0), logical(0)),
                                   stringsAsFactors = FALSE)
    # 보류 block 누설 확인: 보류 칸 (Δ · 검출률 · 검출 수) 을 바꿔도 학습 점수 · 비용이 같아야 한다
    bp <- bk; if (length(ho)) { bp$delta[, ho, ] <- MBVQ_LEAK_VALUE; bp$detrate_delta[, ho, ] <- MBVQ_LEAK_VALUE; bp$det[, ho, ] <- 1e6 }
    ok$leak <- ok$leak && identical(mbvq_variant_scores(mbvq_program_scores(mbvq_gene_effects(bp, tr, mind), PG)), vs)
    rj <- match(paste(ed$to, ed$from, sep = ">"), ed$edge)                       # 역방향 = 같은 자료의 부호 반대
    ok$rev <- ok$rev && !anyNA(rj) && identical(is.na(ed$s_main), is.na(ed$s_main[rj])) && isTRUE(all.equal(ed$s_main, -ed$s_main[rj])) &&
      identical(is.na(ed$s_m02only), is.na(ed$s_m02only[rj])) && isTRUE(all.equal(ed$s_m02only, -ed$s_m02only[rj]))
    # 저장 monocyte 주 경로의 실제 양끝 (후보 풀에서 다시 고르지 않는다; 경로가 없던 key 도 분모에 남긴다)
    st <- en$routes[[MBVQ_ARM]]; st_status <- if (is.null(st$status)) "no_anchor" else st$status; ok0 <- identical(st_status, "route_ok")
    S <- if (ok0) as.character(st$primary[1]) else NA_character_; T <- if (ok0) as.character(st$primary[length(st$primary)]) else NA_character_
    if (ok0 && (!identical(as.character(st$from), S) || !identical(as.character(st$to), T))) stop("저장 경로 양끝이 from · to 와 다르다: ", en$key)
    if (ok0) { n_e3 <- n_e3 + 1L; a0 <- nrv_segment_arrays(x$ss$sums$sb, x$ss$meta$sb, st$primary, x$blocks, genes, minc)    # 기존 E3 정의와 같은 값인가
      j0 <- match(paste(st$primary[-length(st$primary)], st$primary[-1], sep = ">"), ed$edge)
      ok$e3 <- ok$e3 && !anyNA(j0) && isTRUE(all.equal(unname(a0$delta), unname(bk$delta[j0, , , drop = FALSE]))) &&
        isTRUE(all.equal(unname(a0$det), unname(bk$det[j0, , , drop = FALSE]))) }
    paths <- list()
    for (i in seq_len(nrow(V))) { v <- V[i, ]
      s <- vs$s_used[vs$variant == v$variant][match(ed$edge, ps$edge)]; cost <- mbvq_edge_costs(ed$c_graph, s, v$lambda)
      ok$neg <- ok$neg && all(is.finite(cost)) && all(cost >= 0); ed[[paste0("c_expr_", v$variant)]] <- cost
      r <- if (ok0) mbvq_directed_path(data.frame(from = ed$from, to = ed$to, cost = cost, stringsAsFactors = FALSE), S, T) else
        list(status = st_status, path = character(0), cost = NA_real_, n_tied = 0L)
      paths[[v$variant]] <- r$path
      if (v$variant == "lambda0") ok$lam0 <- ok$lam0 && identical(r$status, st_status) && identical(r$path, if (ok0) as.character(st$primary) else character(0)) &&
        (!ok0 || isTRUE(all.equal(r$cost, st$primary_cost, tolerance = 1e-10)))
      if (identical(r$status, "route_ok")) { n_fix <- n_fix + 1L; ok$fix <- ok$fix && identical(r$path[1], S) && identical(r$path[length(r$path)], T) }
      pe <- if (length(r$path) >= 2L) paste(r$path[-length(r$path)], r$path[-1], sep = ">") else character(0); k <- match(pe, ed$edge)
      ok$newe <- ok$newe && !anyNA(k)
      sst <- if (v$score == "irf8_holdout") ed$status_m02only else ed$status_main            # 그 변형이 쓰는 점수의 상태
      RC[[length(RC) + 1L]] <- cbind(kc, variant = v$variant, lambda = v$lambda, score = v$score, role = v$role, status = r$status, S = S, T = T,
        path = paste(r$path, collapse = ">"), stored_status = st_status, stored_path = if (ok0) paste(st$primary, collapse = ">") else "",
        hops = if (length(k)) length(k) else NA_integer_, c_graph_sum = if (length(k)) sum(ed$c_graph[k]) else NA_real_,
        c_expr_sum = if (length(k)) sum(cost[k]) else NA_real_, n_tied = as.integer(r$n_tied),
        n_edges_insufficient = if (length(k)) sum(sst[k] != "ok") else NA_integer_, frac_edges_insufficient = if (length(k)) mean(sst[k] != "ok") else NA_real_,
        min_cell_edges = if (length(k)) min(ed$n_cell_edges[k]) else NA_integer_, n_zero_same_block = if (length(k)) sum(ed$n_same_block_edges[k] == 0L) else NA_integer_,
        min_valid_train_blocks = if (length(k)) min(ed$n_valid_train_blocks[k]) else NA_integer_,
        min_train_boundary_cells = if (length(k)) min(ed$train_boundary_cells[k]) else NA_real_, stringsAsFactors = FALSE)
      kd <- match(dg_edges, ed$edge); kd <- kd[!is.na(kd) & !ed$edge[kd] %in% pe]                    # 진단 선 (선택되지 않았을 때도 남긴다)
      kk <- c(k, kd)
      if (length(kk)) RE[[length(RE) + 1L]] <- cbind(kc, variant = v$variant, lambda = v$lambda,
        row_role = c(rep("selected_path", length(k)), rep(paste0("diagnostic_", paste(dg, collapse = "_")), length(kd))), step = c(seq_along(k), rep(NA_integer_, length(kd))),
        diagnostic_edge = ed$edge[kk] %in% dg_edges, .mbvq_edge_view(ed, kk, ge, genes, s, cost, sst, ed$edge[kk] %in% pe), stringsAsFactors = FALSE)
      if (length(k)) {                                                           # 보류 block 사후 측정 (같은 유전자 · 다른 donor block; 경로 선택 뒤)
        gh <- mbvq_gene_effects(mbvq_ba_subset(bk, pe), ho, mind)
        used <- if (v$lambda == 0) rep(FALSE, nrow(gh)) else if (v$score == "irf8_holdout") gh$gene %in% PG$gene[PG$program == "M02"] else gh$gene %in% PG$gene
        ao <- gh$gene %in% anc; gt <- ge[match(paste(gh$edge, gh$gene), paste(ge$edge, ge$gene)), , drop = FALSE]
        HG[[length(HG) + 1L]] <- cbind(kc, variant = v$variant, lambda = v$lambda, step = match(gh$edge, pe), edge = gh$edge, from = gh$from, to = gh$to, gene = gh$gene,
          program = PG$program[match(gh$gene, PG$gene)], n_heldout_blocks = length(ho), n_valid_blocks_heldout = gh$n_valid_blocks, d_heldout = gh$d,
          detrate_d_heldout = gh$detrate_d, det_sum_heldout = gh$det_sum, evaluable_heldout = gh$evaluable, reason_heldout = gh$reason,
          d_train = gt$d, n_valid_blocks_train = gt$n_valid_blocks, evaluable_train = gt$evaluable, used_in_selection = used, anchor_overlap_gene = ao,
          validation_kind = ifelse(ao, "anchor_dependent_not_independent", ifelse(used, "same_gene_heldout_donor_blocks", "unused_gene_heldout_donor_blocks")),
          note = ifelse(v$score == "irf8_holdout" & gh$gene == "IRF8", "IRF8 excluded from this variant's cost; Harmony/graph dependence not removed", ""),
          stringsAsFactors = FALSE)
      }
      m <- length(r$path); if (m >= 3L) TR[[length(TR) + 1L]] <- cbind(kc, variant = v$variant, position = 2:(m - 1L), A = r$path[1:(m - 2L)], B = r$path[2:(m - 1L)],
                                                                      C = r$path[3:m], stringsAsFactors = FALSE)
    }
    DE[[length(DE) + 1L]] <- cbind(kc, ed, stringsAsFactors = FALSE)
    GE[[length(GE) + 1L]] <- cbind(kc, program = PG$program[match(ge$gene, PG$gene)], sign = PG$sign[match(ge$gene, PG$gene)], ge, stringsAsFactors = FALSE)
    keep[[en$key]] <- list(S = S, T = T, train = tr, held = ho, paths = paths, stored = st)
  }
  rc <- do.call(rbind, RC); l0 <- rc[rc$variant == "lambda0", c("key", "path")]
  rc$changed_vs_lambda0 <- ifelse(rc$status == "route_ok", rc$path != l0$path[match(rc$key, l0$key)], NA)
  nk <- length(x$sel$entries)
  ck$add("boundary_delta_matches_e3", ok$e3, paste0("directed boundary delta / detection arrays on every stored monocyte primary segment equal nrv_segment_arrays (existing E3) for ", n_e3, " keys"))
  ck$add("lambda0_reproduces_stored_route", ok$lam0, paste0("lambda0 directed path (status, path, cost) == stored MBV monocyte primary for ", nk, " keys (lexicographic tie rule)"))
  ck$add("fixed_endpoints", ok$fix, paste0("every variant route starts at the stored S and ends at the stored T of its key (", n_fix, " routes); endpoints not reselected"))
  ck$add("no_negative_costs", ok$neg, "every directed cost c_graph * exp(-lambda * s_used) is finite and >= 0")
  ck$add("no_new_edges", ok$newe, "every selected route edge is an observed pair of that key's c1 fold table")
  ck$add("heldout_not_used_for_costs", ok$leak, paste0("training scores and all variant costs are identical after replacing held-out block delta / detection with ", MBVQ_LEAK_VALUE, " / 1e6"))
  ck$add("reverse_direction_flips_sign", ok$rev, "s_main(B>A) == -s_main(A>B) and s_M02(B>A) == -s_M02(A>B) with identical NA pattern for every observed pair")
  log("keys:", nk, "| route_ok by variant:", paste(names(table(rc$variant[rc$status == "route_ok"])), table(rc$variant[rc$status == "route_ok"]), collapse = ", "),
      "| changed vs lambda0:", paste(tapply(rc$changed_vs_lambda0 %in% TRUE, rc$variant, sum)[V$variant], collapse = "/"))
  list(routes = rc, route_edges = if (length(RE)) do.call(rbind, RE) else RE0, edge_scores = do.call(rbind, DE), gene_effects = do.call(rbind, GE),
       heldout = if (length(HG)) do.call(rbind, HG) else NULL, triples = if (length(TR)) do.call(rbind, TR) else NULL, keep = keep)
}


# ── 4. 새 경로 내부 triple: 위상 구간만 (기존 NRV 세포 edge, 전체 세포). RNA 기울기는 계산하지 않는다 ─────
.mbvq_span <- function(x, tr, ck, log) {
  if (is.null(tr) || !nrow(tr)) { ck$add("internal_span_matches_stored_pooled", TRUE, "no interior triples"); return(NULL) }
  e <- readRDS(x$inp[["edges"]]); cl <- x$cell_cluster; ci <- cl[e$i]; cj <- cl[e$j]
  u <- unique(tr[, c("B", "A", "C")]); res <- vector("list", nrow(u))
  for (r in seq_len(nrow(u))) { abc <- c(u$A[r], u$B[r], u$C[r]); ic <- nrv_internal_coordinate(e[ci %in% abc & cj %in% abc, , drop = FALSE], cl, u$B[r], u$A[r], u$C[r])
    res[[r]] <- data.frame(B = u$B[r], A = u$A[r], C = u$C[r], topological_status = ic$status, n_cells_B = sum(cl == u$B[r]), n_shared = ic$n_shared,
                           n_unreachable = ic$n_unreachable, stringsAsFactors = FALSE) }
  res <- do.call(rbind, res); id <- paste(res$B, res$A, res$C, sep = "|")
  sto <- vapply(id, function(k) { z <- x$sel$icoord[[k]]; if (is.null(z)) NA_character_ else as.character(z$status) }, "", USE.NAMES = FALSE)   # 같은 계산인지 확인만 (복사하지 않는다)
  ck$add("internal_span_matches_stored_pooled", all(is.na(sto) | sto == res$topological_status),
         paste0("recomputed pooled topological status equals the stored select icoord status for ", sum(!is.na(sto)), " of ", nrow(res), " triples (others are new triples)"))
  sp <- merge(tr, res, by = c("B", "A", "C"), sort = FALSE)
  sp <- sp[order(sp$key, match(sp$variant, MBVQ_VARIANTS$variant), sp$position, method = "radix"),
           c("partition", "key", "holdout_scheme", "fold", "variant", "position", "A", "B", "C", "topological_status", "n_cells_B", "n_shared", "n_unreachable")]
  sp$rna_gradient <- "not_evaluated"; sp$interpretation <- "topological_span_only_not_internal_expression_continuity"; rownames(sp) <- NULL
  log("interior triples:", nrow(res), "| status:", paste(names(table(res$topological_status)), table(res$topological_status), collapse = ", "))
  sp
}


# ── 5. 사후 저자 대조 (여기서 처음 저자 라벨 표를 읽는다) ────────────────────
.mbvq_author <- function(x) {
  inp <- c(author_counts = .mbvq_stage_file(x$m$nrv_attempt, "author", "author_cluster_label_counts.tsv.gz"),
           broad_map = .mbvq_stage_file(x$attempt, "author", "broad_label_mapping.csv"))
  cnt <- utils::read.delim(inp[["author_counts"]], stringsAsFactors = FALSE, colClasses = c(cluster = "character", author_label = "character"))
  bm <- utils::read.csv(inp[["broad_map"]], stringsAsFactors = FALSE)
  cnt <- cnt[order(cnt$cluster, -cnt$n, cnt$author_label, method = "radix"), ]; top <- cnt[!duplicated(cnt$cluster), ]; tot <- tapply(cnt$n, cnt$cluster, sum)
  lab <- stats::setNames(top$author_label, top$cluster); frac <- stats::setNames(top$n / as.numeric(tot[top$cluster]), top$cluster)
  fam <- stats::setNames(bm$broad_family[match(top$author_label, bm$author_label)], top$cluster)
  short <- function(k) ifelse(is.na(k) | is.na(lab[k]), NA_character_, sprintf("%s %.0f%%", lab[k], 100 * frac[k]))
  list(lab = lab, frac = frac, fam = fam, short = short, inputs = inp)
}
.mbvq_sp <- function(p) if (is.na(p) || !nzchar(p)) character(0) else strsplit(p, ">", fixed = TRUE)[[1]]
.mbvq_author_cols <- function(rc, au) {
  rc$path_labels <- vapply(rc$path, function(p) paste(vapply(.mbvq_sp(p), function(k) paste(k, au$short(k)), ""), collapse = " > "), "", USE.NAMES = FALSE)
  rc$endpoint_label <- unname(au$short(rc$T))
  rc$n_gmp_family_clusters <- vapply(rc$path, function(p) sum(au$lab[.mbvq_sp(p)] %in% MBVQ_GMP_FAMILY), 0L, USE.NAMES = FALSE)
  rc$n_monocyte_lineage_clusters <- vapply(rc$path, function(p) sum(au$fam[.mbvq_sp(p)] %in% "monocyte_lineage"), 0L, USE.NAMES = FALSE)
  rc
}


# ── 6. 조건 요약 (고정 key 분모, 변화 / 유지 / no_call, 사후 라벨 서술) ─────────
.mbvq_summary <- function(rc, dg) {
  cnt_str <- function(v) { v <- v[!is.na(v) & nzchar(v)]; if (!length(v)) return(""); t <- table(v); o <- order(-as.integer(t), names(t), method = "radix")
    paste(paste0(names(t)[o], " (", as.integer(t)[o], ")"), collapse = "; ") }
  med <- function(v) suppressWarnings(stats::median(v, na.rm = TRUE))
  l1 <- rc[rc$variant == "lambda1", ]
  do.call(rbind, lapply(MBVQ_VARIANTS$variant, function(vn) { z <- rc[rc$variant == vn, ]; y <- l1[match(z$key, l1$key), ]
    data.frame(variant = vn, lambda = z$lambda[1], score = z$score[1], role = z$role[1], n_keys = nrow(z), n_route_ok = sum(z$status == "route_ok"),
               n_no_path = sum(z$status == "no_path"), n_anchor_overlap = sum(z$status == "anchor_overlap"), n_no_anchor = sum(z$status == "no_anchor"),
               n_changed_vs_lambda0 = sum(z$changed_vs_lambda0 %in% TRUE), n_same_as_lambda0 = sum(z$changed_vs_lambda0 %in% FALSE), n_no_call = sum(z$status != "route_ok"),
               n_same_path_as_lambda1 = sum(z$status == "route_ok" & z$path == y$path), paths = cnt_str(z$path), median_hops = med(z$hops),
               median_c_graph_sum = med(z$c_graph_sum), median_c_expr_sum = med(z$c_expr_sum), median_frac_edges_insufficient = med(z$frac_edges_insufficient),
               n_keys_via_diag_edge = if (length(dg) == 2L) sum(grepl(paste0("(^|>)", dg[1], ">", dg[2], "(>|$)"), z$path)) else NA_integer_,
               n_keys_via_gmp_family_label = sum(z$n_gmp_family_clusters > 0L), n_keys_with_monocyte_lineage_label = sum(z$n_monocyte_lineage_clusters > 0L),
               endpoints = cnt_str(ifelse(z$status == "route_ok", paste(z$T, z$endpoint_label), NA)), stringsAsFactors = FALSE) }))
}


# ── 7. 실행 순서 ─────────────────────────────────────────────────────────────
mbvq_build <- function(a, out, snap, log, argv) mbvq_build_from(.mbvq_inputs(a), a, out, snap, log, argv)

#' 읽은 입력 (.mbvq_inputs) 으로 계산 · 표 · 그림 · 확인 · 출처를 쓴다 (경계 조건 테스트가 입력을 직접 넘길 수 있게 나눴다)
mbvq_build_from <- function(x, a, out, snap, log, argv) {
  ck <- .mbvq_checks()
  PG <- mbvq_programs(x$fr$ref$edges, unique(MBVQ_EXPECTED_PROGRAMS$program))
  pm <- identical(PG$program, MBVQ_EXPECTED_PROGRAMS$program) && identical(PG$gene, MBVQ_EXPECTED_PROGRAMS$gene) &&
    isTRUE(all(PG$sign == MBVQ_EXPECTED_PROGRAMS$sign)) && identical(PG$evidence, MBVQ_EXPECTED_PROGRAMS$evidence)
  ck$add("programs_match_reference", pm, paste0("frozen reference edges give ", paste(PG$program, PG$gene, ifelse(PG$sign > 0, "+", "-"), PG$evidence, collapse = "; "),
                                                " (contract: M01 IRF8+, M02 CD14+ CCR2+ FCGR1A+, conditional)"))
  if (!pm) stop("참고표 프로그램이 계약과 다르다: ", paste(PG$program, PG$gene, collapse = ", "))
  anc <- unique(c(x$fr$anchors[["S"]], x$fr$anchors[[MBVQ_ARM]])); PG$anchor_overlap_gene <- PG$gene %in% anc     # S · monocyte T anchor panel 과 겹치는 유전자
  genes <- unique(PG$gene); miss <- setdiff(genes, colnames(x$ss$sums$sb$sum)); if (length(miss)) stop("경계 집합 합에 없는 프로그램 유전자다: ", paste(miss, collapse = ", "))
  pa <- unique(do.call(rbind, lapply(x$pr$folds, function(pt) data.frame(a = as.character(pt$cluster_a), b = as.character(pt$cluster_b), stringsAsFactors = FALSE))))
  de <- unique(data.frame(from = c(pa$a, pa$b), to = c(pa$b, pa$a), stringsAsFactors = FALSE))
  ba <- mbvq_boundary_arrays(x$ss$sums$sb, x$ss$meta$sb, de, x$blocks, genes, .nrv_thr("set_min_cells"))
  log("programs:", paste(PG$program, PG$gene, collapse = ", "), "| anchor-overlap genes:", paste(PG$gene[PG$anchor_overlap_gene], collapse = ", "),
      "| directed edges:", nrow(de), "| blocks:", length(x$blocks), "| boundary cells ok:", sum(ba$reason == "ok"), "of", length(ba$reason))
  kr <- .mbvq_keys(x, ba, PG, anc, ck, log)
  sp <- .mbvq_span(x, kr$triples, ck, log)
  au <- .mbvq_author(x)                                                          # 경로 · 점수 · 구간 계산 뒤에만
  rc <- .mbvq_author_cols(kr$routes, au); dg <- MBVQ_DIAG_PAIRS[[x$partition]]
  re <- kr$route_edges; re$from_label <- unname(au$short(re$from)); re$to_label <- unname(au$short(re$to))
  su <- .mbvq_summary(rc, dg)
  fg <- .mbvq_figures(x, a, out, rc, kr, au, PG, ck, log)
  files <- fg$files
  tabs <- list(route_comparison.tsv = rc, route_edges.tsv = re, directed_edge_scores.tsv.gz = kr$edge_scores, gene_effects.tsv.gz = kr$gene_effects,
               heldout_gene_effects.tsv.gz = kr$heldout, internal_span.tsv = sp, control_summary.tsv = su, programs.tsv = PG, representative.tsv = fg$representative,
               figure_text.tsv = fg$text)
  tabs <- Filter(Negate(is.null), tabs); not_written <- setdiff(c("heldout_gene_effects.tsv.gz", "internal_span.tsv"), names(tabs))   # 경로가 없으면 행도 없다
  rtp <- vapply(names(tabs), function(f) { p <- file.path(out, f); files <<- c(files, p); mbvt_write_tsv(tabs[[f]], p) }, logical(1))
  ck$add("tables_roundtrip_ok", all(rtp), "each TSV (incl. .gz) re-read with the same rows, field counts and text columns")
  ck$add("author_labels_read_after_routes", TRUE, "author label tables are opened only in step 5, after scores, costs, routes and interior spans (code order); figure colours read after that")
  cks <- ck$table(); p <- file.path(out, "checks.tsv"); mbvt_write_tsv(cks, p); files <- c(files, p)
  p <- file.path(out, "mono_quant_lossless.rds")
  saveRDS(list(contract = "codex_mono_quant_instruction.md (2026-10-15)", partition = x$partition, programs = PG, variants = MBVQ_VARIANTS,
               settings = list(set_min_cells = .nrv_thr("set_min_cells"), gene_min_detected_cells = .nrv_thr("gene_min_detected_cells"), min_genes = MBVQ_MIN_GENES,
                               tie_tol = MBVQ_TOL, leak_value = MBVQ_LEAK_VALUE, diag_pairs = dg),
               boundary_arrays = ba, keys = kr$keep, route_comparison = rc, route_edges = re, edge_scores = kr$edge_scores, gene_effects = kr$gene_effects,
               heldout = kr$heldout, internal_span = sp, control_summary = su, representative = fg$representative, checks = cks), p)
  files <- c(files, p)
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  rp <- fg$representative
  s1 <- su[su$variant == "lambda1", ]; s0 <- su[su$variant == "lambda0", ]
  d5 <- re[re$diagnostic_edge %in% TRUE & re$variant == "lambda1", , drop = FALSE]
  rd <- c("# monocyte 정량 발현 재가중치 비교 (탐색적, 비용만 바꿈)", "",
          "- 계약: results/claude_bridge/myeloid_branch_validation_2026-10-13/codex_mono_quant_instruction.md (2026-10-15, 사용자 승인 · Codex 설계). 원 RNA 는 읽지 않았다.",
          "- 같은 graph · 군집 · 11 key · 학습 / 보류 block. key 마다 저장 monocyte 주 경로의 실제 S, T 를 그대로 두고, c1 baseline 쌍 표의 방향별 비용만 바꿨다.",
          "  r ≥ 1 필터 · 오르막 끝점 · 중간 경유점은 쓰지 않았다. 그 fold 의 모든 관측 선 위에서 방향 최단 단순 경로를 찾았다 (동률 = 경로 문자열 사전순).",
          "- 경계 변화 Δ_g^b(A→B) = 평균 x(같은 block 의 B 쪽 경계 세포) − 평균 x(A 쪽 경계 세포), x = log1p(10^4 · count / 총 count). 양쪽 ≥ 5 세포인 학습 block 만.",
          "  d_g = 유효 학습 block Δ 의 동일 가중 평균. 유전자 학습 경계 검출 합 < 5 이거나 유효 block 이 없으면 NA 와 이유 (0 으로 두지 않는다).",
          "- s_M01 = tanh(d_IRF8) (단일 유전자), s_M02 = 평가 가능한 CD14 · CCR2 · FCGR1A 의 tanh(d_g) 평균 (2/3 이상). s = (s_M01 + s_M02) / 2.",
          "  한 프로그램이라도 없으면 s_used = 0 (비용 유지), raw s = NA, 상태 insufficient_boundary_evidence. 결측에 낮은 비용이나 생물학 지지를 준 것이 아니다.",
          "- 비용 c_expr(A→B) = c_graph(A,B) · exp(−λ · s_used(A→B)), c_graph = −log q 그대로. λ = 1 주, 0.5 · 2 사전 민감도, 0 = 저장 경로 회귀 확인.",
          "  반대 방향 대조 = λ1 에서 −s. IRF8 제외 진단 = λ1 에서 s = s_M02 (IRF8 은 그 변형의 비용에서 빠졌지만 원 Harmony · graph 의존은 남는다).",
          "- 이 점수는 경계 사이 정량 차이에 대한 탐색적 단핵구 방향 선호다. 시간 미분 · 분화 법칙 · 단계 배정이 아니고, 확률도 아니다. 경로 비용은 여전히 선 비용의 합이다 (hop 수 편향을 해결했다고 주장하지 않는다).",
          paste0("- 결과 (key 수, 분모 ", s1$n_keys, "): lambda1 경로가 lambda0 과 다른 key ", s1$n_changed_vs_lambda0, ", 같은 key ", s1$n_same_as_lambda0, ", no_call ", s1$n_no_call,
                 ". 변형별 변화 수: ", paste(sprintf("%s %d", su$variant, su$n_changed_vs_lambda0), collapse = ", "), "."),
          paste0("- lambda0 경로 (저장 경로): ", s0$paths), paste0("- lambda1 경로: ", s1$paths),
          if (nrow(d5)) paste0("- 진단 선 ", paste(unique(d5$edge), collapse = " · "), " (lambda1, 관측 key ", length(unique(d5$key)), "/", s1$n_keys, "): 점수 상태 ",
                               paste(names(table(d5$score_status)), table(d5$score_status), collapse = ", "), ", 유효 학습 block 수 최대 ", max(d5$n_valid_train_blocks),
                               ". 경계 근거가 부족하면 비용은 그대로다 (route_edges.tsv 의 diagnostic_edge 행).")
          else if (length(dg)) paste0("- 진단 선 ", paste(dg, collapse = "–"), ": 어느 key 의 fold 쌍 표에도 관측되지 않았다.") else "- 진단 선: 이 분할에는 지정 진단 쌍이 없다 (군집 번호가 p92 기준).",
          paste0("- 대표 key (결과를 보기 전에 정한 규칙: 저장 경로 최빈값을 쓰는 key 중 이름 사전순 첫 번째): ", rp$key[1], ". 그림 두 패널은 같은 key 다."),
          "- 원저자 라벨은 경로 · 점수 · 내부 구간 계산 뒤 사후 대조 (route_comparison 의 path_labels · endpoint_label 등) 와 그림 색에만 썼다. 성공을 자동 판정하지 않는다.",
          "- 새 경로 내부 triple 은 기존 NRV 세포 edge 의 위상 구간 (internal_span.tsv) 만 계산했다. 내부 발현 연속성이 아니다. RNA 기울기는 not_evaluated (새 gene-null 수집 없음).",
          "- 보류 block 측정 (heldout_gene_effects.tsv.gz): 선택에 쓴 유전자의 다른 donor block 값이다 (same_gene_heldout_donor_blocks). M02 유전자는 T anchor 와 겹쳐 독립 검증이 아니다",
          "  (anchor_dependent_not_independent). 보류 block 은 key 마다 2–20 개이고 11 key 는 독립 표본 11 개가 아니다.",
          "- 하지 않은 것: 새 gene-null, 실제 분화 · in-vivo fate · 자동 주석 정확도 개선의 입증, 기존 annotation 덮어쓰기, λ 를 결과 보고 고르기.",
          if (length(not_written)) paste0("- 경로 (또는 내부 triple) 가 없어 행이 없으므로 만들지 않은 표: ", paste(not_written, collapse = ", "), ". no_call 은 route_comparison · control_summary 분모에 남아 있다."),
          "- 그림: figures/fig_mono_quant_two_panel.png (왼쪽 lambda0 = 저장 경로, 오른쪽 lambda1, 같은 key · 같은 좌표 · 같은 창),",
          "  figures/fig_mono_quant_edge_effects.png (대표 key 의 두 경로 선 + 진단 선: 유전자별 학습 변화 d_g, 프로그램 점수, 원 비용 대 새 비용). 그린 글자는 figure_text.tsv.")
  p <- file.path(out, "README_ko.md"); writeLines(rd, p); files <- c(files, p)
  inp <- c(x$inp, au$inputs, fg$inputs)
  pv <- data.frame(input = names(inp), path = unname(inp), sha256 = vapply(inp, .nrv_sha, "", USE.NAMES = FALSE), stringsAsFactors = FALSE)
  p <- file.path(out, "provenance.tsv"); mbvt_write_tsv(pv, p); files <- c(files, p)
  pk <- c("igraph", "Matrix", "jsonlite", "digest"); pv_of <- function(k) tryCatch(as.character(utils::packageVersion(k)), error = function(err) "not_installed")
  env <- data.frame(item = c("R_executable", "R", "platform", "libPaths", paste0("package:", pk), "title_font", "random_seed", "lambdas", "set_min_cells",
                             "gene_min_detected_cells", "program_min_genes", "tie_tol", "leak_check_value", "command", "working_dir", "written_at"),
                    value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"), vapply(pk, pv_of, ""),
                              if (fg$ko) MBVR_TITLE_FAMILY else "default", "none (no random numbers are drawn)", paste(unique(MBVQ_VARIANTS$lambda), collapse = ";"),
                              .nrv_thr("set_min_cells"), .nrv_thr("gene_min_detected_cells"), paste(names(MBVQ_MIN_GENES), MBVQ_MIN_GENES, sep = "=", collapse = ";"),
                              format(MBVQ_TOL), MBVQ_LEAK_VALUE, paste(c("46_run_mbv_mono_quant.R", argv), collapse = " "), getwd(), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
                    stringsAsFactors = FALSE)
  p <- file.path(out, "environment.tsv"); mbvt_write_tsv(env, p); files <- c(files, p)
  p <- file.path(out, "source_snapshot.tsv"); mbvt_write_tsv(snap, p); files <- c(files, p)
  relf <- substring(normalizePath(files), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(output = basename(out), attempt = x$attempt, written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), relf))), auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "), "| lambda1 changed vs lambda0:", s1$n_changed_vs_lambda0, "/", s1$n_keys)
}


# 스크립트로 실행할 때: bootstrap 에 필요한 파일만 먼저 읽는다. 계산은 스냅샷 사본에서 한다.
if (!interactive() && sys.nframe() == 0L) {
  .mbvq_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "27_run_nrv.R"))
    source(file.path(.mbvq_bb, .f))
  quit(status = mbvq_main(commandArgs(trailingOnly = TRUE)))
}
