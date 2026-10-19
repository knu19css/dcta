# =============================================================================
# 26_nrv_evaluate.R — 경계 대비, block 추론, 내부 기울기, 파동 서술, anchor 후보, gene null, 판정
# =============================================================================
#
# 추론의 두 가지 대조 (따로 보고한다)
#   - donor 대조: donor block 이 반복 단위다. 적격 block 의 부호 검정(R17 mdu_sign_test)과 equal-block 평균의 block 재표본 CI.
#     가정: block 은 서로 독립·교환 가능.
#   - gene 대조: 같은 세포·같은 경로에 bin 을 맞춘 무작위 유전자를 넣고, 데이터를 보는 선택(평가 가능성, 구간, 내부 프로그램)을 똑같이 반복한다.
#     가정: bin 안 유전자 교환 가능. 보편적으로 보정된 오류 확률이 아니다.
#   block 수로 미리 거르지 않는다. 계산 가능한 p 를 모두 낸 뒤 family 에 Holm 을 적용한다.
#   경로 내부 위치 순열이나 세 무리 이항(p0 = 1/3) 같은 교환성 없는 교란은 쓰지 않는다. 세 무리 평균은 서술만 한다.
# =============================================================================

#' 경계 대비: 하류(to) 평균 − 상류(from) 평균, block x 유전자. 세포가 min_cells 보다 적은 block 은 NA (행은 지우지 않는다)
nrv_boundary_contrast <- function(sum_from, n_from, sum_to, n_to, min_cells = 5) {
  d <- sum_to / as.numeric(n_to) - sum_from / as.numeric(n_from)
  d[as.numeric(n_from) < min_cells | as.numeric(n_to) < min_cells, ] <- NA; d
}

#' 프로그램 대비: 부호 평균 (평가 가능 유전자만 넘긴다). 같은 유전자가 다른 프로그램에 있어도 프로그램마다 따로 계산한다.
nrv_program_delta <- function(delta, genes, signs) {
  setNames(as.numeric(delta[, genes, drop = FALSE] %*% (signs / length(genes))), rownames(delta))
}

#' block 추론: equal-block 평균, block 재표본 95 % CI, 정확 부호 검정. n_for_strictest_holm_threshold 는 설명 열이다(거르지 않는다).
nrv_block_inference <- function(x, R = 2000L, seed = 20261013L, family_size = 1L, alpha = 0.05) {
  v <- x[!is.na(x)]; n <- length(v)
  ci <- if (n >= 2) .nrv_with_seed(seed, stats::quantile(replicate(R, mean(sample(v, n, replace = TRUE))), c(0.025, 0.975), names = FALSE)) else c(NA_real_, NA_real_)
  list(estimate = if (n) mean(v) else NA_real_, ci_low = ci[1], ci_high = ci[2], n_eligible = n, sign_test = mdu_sign_test(v),
       n_for_strictest_holm_threshold = ceiling(log2(2 * family_size / alpha)))
}

#' 내부 기울기: block 마다 B 세포의 프로그램 점수 s 를 위치 t 에 OLS 회귀한 기울기, 그리고 block 추론
nrv_internal_slope <- function(t, s, block, min_cells = 5) {
  ok <- !is.na(t) & !is.na(s); sp <- split(which(ok), block[ok])
  pb <- do.call(rbind, lapply(names(sp), function(b) { i <- sp[[b]]
    if (length(i) < min_cells || length(unique(t[i])) < 2) return(NULL)
    data.frame(block = b, n_cells = length(i), slope = stats::cov(t[i], s[i]) / stats::var(t[i]), stringsAsFactors = FALSE) }))
  if (is.null(pb)) pb <- data.frame(block = character(0), n_cells = integer(0), slope = numeric(0))
  list(per_block = pb, inference = nrv_block_inference(setNames(pb$slope, pb$block)))
}

#' 세 무리(A 쪽 t < 1/3, 가운데, C 쪽 t ≥ 2/3) 평균. p 값은 없다(교환 가능하지 않으므로 서술만).
nrv_between_descriptive <- function(t, s, block) {
  g <- cut(t, c(-Inf, 1 / 3, 2 / 3, Inf), labels = c("A_side", "middle", "C_side"), right = FALSE)
  m <- tapply(s, list(block, g), mean)
  data.frame(block = rownames(m), mean_A_side = m[, "A_side"], mean_middle = m[, "middle"], mean_C_side = m[, "C_side"], row.names = NULL)
}

#' 파동 서술: site 별 equal-block 평균, block 을 함께 재표본한 정점 빈도, CI 기반 진입·퇴장 site.
#' 끝점 정규화 50 % 시작점은 쓰지 않는다. 탐색적 서술량이며 동시구간 보정 시간 추정치가 아니다.
#' 결측 (Codex 재현 수정)
#'   - block x site 의 NA 는 그 칸만 뺀다. n_blocks_site 에 site 별 값이 있는 block 수를 남긴다.
#'   - 재표본에서 값이 있는 site 가 하나도 없으면 그 반복은 정점 없음으로 센다(n_draws_no_peak). 정점 빈도의 분모는 R 그대로다.
#'   - 진입·퇴장 CI 에서 두 site 를 비교할 수 없는 반복은 -Inf 로 두어 근거로 세지 않는다(보수적).
nrv_wave_profile <- function(Y, direction, R = 2000L, seed = 20261013L) {
  site <- colnames(Y); nb <- setNames(as.integer(colSums(!is.na(Y))), site); none <- setNames(rep(NA_real_, length(site)), site)
  if (!nrow(Y) || all(is.na(Y))) return(list(status = "all_missing", site_mean = none, n_blocks_site = nb, peak_freq = none,
                                              n_draws_no_peak = as.integer(R), R = as.integer(R), entry_site = NA_character_, exit_site = NA_character_,
                                              site_ci_low = none, site_ci_high = none, n_draws_site_missing = setNames(rep(as.integer(R), length(site)), site)))
  cm <- function(M) { m <- colMeans(M, na.rm = TRUE); m[is.nan(m)] <- NA_real_; m }
  B <- .nrv_with_seed(seed, replicate(R, cm(Y[sample.int(nrow(Y), nrow(Y), replace = TRUE), , drop = FALSE])))
  dirB <- direction * matrix(B, nrow = R, ncol = length(site), byrow = TRUE); y <- cm(Y)
  pk_draw <- apply(dirB, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))
  peak <- tabulate(pk_draw[!is.na(pk_draw)], length(site)) / R
  lo <- function(D) apply(D, 2, function(x) { x[is.na(x)] <- -Inf; stats::quantile(x, 0.025, names = FALSE) })
  ent <- which(lo(dirB - dirB[, 1]) > MDU_NUM_TOL); ent <- ent[ent > 1]
  pk <- which.max(direction * y); ex <- which(lo(dirB[, pk] - dirB) > MDU_NUM_TOL); ex <- ex[ex > pk]
  Bs <- dirB / direction                                                   # 원래 척도의 site 평균 재표본 (그림용 서술 CI)
  ci <- apply(Bs, 2, function(x) if (all(is.na(x))) c(NA_real_, NA_real_) else stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  list(status = "ok", site_mean = y, n_blocks_site = nb, peak_freq = setNames(peak, site), n_draws_no_peak = sum(is.na(pk_draw)), R = as.integer(R),
       entry_site = if (length(ent)) site[min(ent)] else NA_character_, exit_site = if (length(ex)) site[min(ex)] else NA_character_,
       site_ci_low = setNames(ci[1, ], site), site_ci_high = setNames(ci[2, ], site), n_draws_site_missing = setNames(as.integer(colSums(is.na(Bs))), site))
}

#' anchor 점수: (군집 x block) 점수의 equal-block 평균과, 공통 block 추출 재표본에서의 1 위 빈도.
#' 최고점은 그 상태의 존재 증거가 아니다 → anchor_status 는 candidate_only (적격 block 부족이면 anchor_ineligible).
#' 동률 (Codex 재현 수정): 1 위 점수와 tol(승자 동률 1e-12) 안의 후보는 그 반복의 질량 1 을 똑같이 나눈다. 행 순서로 승자를 고르지 않는다.
#'   적격 후보 점수가 하나도 없는 반복은 n_draws_no_candidate 로 세고, top_freq 의 분모는 R 그대로다.
nrv_anchor_scores <- function(score_cb, n_cb, min_cells = 20, min_blocks = 3, R = 2000L, seed = 20261013L, tol = MDU_TIE_TOL) {
  el <- n_cb >= min_cells & !is.na(score_cb); ne <- rowSums(el)
  sc <- vapply(seq_len(nrow(score_cb)), function(k) if (ne[k]) mean(score_cb[k, el[k, ]]) else NA_real_, numeric(1))
  ok <- ne >= min_blocks; S0 <- score_cb; S0[!el] <- 0
  mass <- .nrv_with_seed(seed, replicate(R, {
    w <- tabulate(sample.int(ncol(score_cb), ncol(score_cb), replace = TRUE), ncol(score_cb))   # 공통 block 추출: 공분산 보존
    den <- as.numeric(el[ok, , drop = FALSE] %*% w); d <- as.numeric(S0[ok, , drop = FALSE] %*% w) / den; d[den == 0] <- NA
    m <- numeric(nrow(score_cb))
    if (any(!is.na(d))) { top <- !is.na(d) & d >= max(d, na.rm = TRUE) - tol; m[which(ok)[top]] <- 1 / sum(top) }
    m }))
  mass <- matrix(mass, nrow = nrow(score_cb))
  data.frame(cluster = rownames(score_cb), score = sc, n_eligible_blocks = as.integer(ne), top_freq = rowSums(mass) / R,
             n_draws_no_candidate = sum(colSums(mass) == 0), anchor_status = ifelse(ok, "candidate_only", "anchor_ineligible"), stringsAsFactors = FALSE)
}

#' 학습 값으로 구간을 고른 뒤(문헌 방향으로 가장 큰 구간, nrv_pick_max) 보류 block 에서 검정한다. 선택은 보류 값을 보지 않는다.
#' ... 는 nrv_block_inference 로 넘긴다(R, seed, family_size).
nrv_select_then_test <- function(train_seg, heldout_seg_block, direction, ...) {
  sel <- nrv_pick_max(train_seg, direction)
  list(selected = sel, heldout = if (is.na(sel)) NULL else nrv_block_inference(direction * heldout_seg_block[sel, ], ...))
}

# ── 판정 조립 ───────────────────────────────────────────────────────────────
#' 필수 프로그램: 층마다 eval 유전자가 하나라도 있는 프로그램 전부 (E10 제외). 유전자 1 개면 단일 유전자 근거로 따로 평가한다.
nrv_required_programs <- function(programs, min_genes = 2) {
  k <- stats::aggregate(list(n_eval_genes = rep(1L, nrow(programs))), programs[, c("program_id", "edge_id")], sum)
  k <- k[k$edge_id != "E10", ]; k$layer <- unname(NRV_EDGE_LAYER[k$edge_id]); k$program_evaluable_by_rule <- k$n_eval_genes >= min_genes
  k[order(k$program_id, method = "radix"), ]
}

#' 순서 포함 범위: 필수 프로그램 전부가 배정 구간에서 지지되고 위치가 문헌 순서와 충돌하지 않을 때만 complete.
#' 비었거나, 평가 불가가 있거나, 위치·순위가 NA 이면 incomplete (빈집합 참 금지).
nrv_coverage <- function(required, phase_status, phase_position, lit_rank) {
  if (!length(required)) return(list(status = "incomplete", reasons = "no_required_programs"))
  st <- phase_status[required]; okst <- c("phase_supported", "phase_supported_single_gene")
  if (anyNA(st) || any(!st %in% okst)) return(list(status = "incomplete", reasons = paste0("not_supported:", required[is.na(st) | !st %in% okst])))
  pos <- phase_position[required]; rk <- lit_rank[required]
  if (anyNA(pos) || anyNA(rk)) return(list(status = "incomplete", reasons = "phase_position_or_rank_missing"))
  bad <- outer(rk, rk, "<") & outer(pos, pos, ">")                         # 문헌상 앞선 프로그램이 뒤에 오면 충돌 (동률 허용)
  if (any(bad, na.rm = TRUE)) return(list(status = "order_conflict", reasons = "phase_order_conflict"))
  list(status = "complete_order_consistent", reasons = character(0), single_gene = required[st == "phase_supported_single_gene"])
}

#' full route 라벨: 다섯 근거를 모두 요구한다. 최고 라벨은 supported_anchor_conditional (anchor 는 후보일 뿐).
nrv_full_route_label <- function(structural, terminal_ok, internal_status, segment_observed, required_phase_status, coverage_by_layer,
                                 block_structure_all_span_ok) {
  rc <- character(0)
  if (structural != "route_ok") rc <- c(rc, structural)
  if (!isTRUE(terminal_ok)) rc <- c(rc, "terminal_anchor_unavailable")
  if (!length(internal_status)) rc <- c(rc, "internal_not_applicable") else if (any(internal_status != "internally_continuous"))
    rc <- c(rc, paste0("internal:", unique(internal_status[internal_status != "internally_continuous"])))
  if (!length(segment_observed) || !all(segment_observed)) rc <- c(rc, "unobserved_segment")
  cov <- vapply(coverage_by_layer, function(x) x$status, character(1))
  if (any(cov != "complete_order_consistent")) rc <- c(rc, paste0("coverage_", names(cov)[cov != "complete_order_consistent"], ":", cov[cov != "complete_order_consistent"]))
  conflict <- any(required_phase_status == "phase_reversed") || any(cov == "order_conflict") || any(internal_status == "internally_reversed")
  label <- if (conflict) "conflicting" else if (length(rc)) "insufficient_evidence" else if (isTRUE(block_structure_all_span_ok)) "supported_anchor_conditional" else "supported_anchor_conditional_pooled_geometry"
  sg <- unlist(lapply(coverage_by_layer, function(x) x$single_gene))
  if (startsWith(label, "supported") && length(sg)) label <- paste0(label, "_single_gene_component")
  list(label = label, reason_codes = rc, single_gene_programs = unique(sg))
}

#' 문헌 edge 10 행. 대응(mapping)은 nrv_combined_mapping 의 계산 결과를 따른다. 형태 단계 정확 대응은 not_confirmed.
#' candidate_stage_set = 대응한 관측 부분경로(subroute; 여러 개면 "|"), n_observed_edges = 그 부분경로의 관측 edge 수 (Codex 재현 수정: 열 계약)
nrv_edge_assessment <- function(ref, roles, program_results, mapping = NULL) {
  e <- ref$edges
  if (!is.null(mapping) && !all(c("edge_id", "mapping_status", "subroute", "n_observed_edges", "combined_status") %in% names(mapping)))
    stop("대응 표에 nrv_combined_mapping 의 열(edge_id, mapping_status, subroute, n_observed_edges, combined_status)이 없다.")
  do.call(rbind, lapply(seq_len(nrow(e)), function(k) {
    id <- e$edge_id[k]; r <- roles[roles$edge_id == id, ]; pr <- program_results[startsWith(program_results$program_id, paste0(id, "_")), ]
    mp <- if (!is.null(mapping)) mapping[mapping$edge_id == id, ] else NULL
    ms <- if (id == "E10") "out_of_scope" else if (!is.null(mp) && nrow(mp)) mp$mapping_status[1] else "unmapped"
    circ <- nrow(r) > 0 && all(r$role == "anchor_circular"); cs <- if (!is.null(mp) && nrow(mp)) mp$combined_status[1] else NA_character_
    verdict <- if (id == "E10") "out_of_scope" else if (circ || !nrow(pr)) "insufficient_evidence" else if (ms != "single") "unresolved_mapping" else
      if (identical(cs, "conflicting")) "conflicting" else if (!is.na(cs) && startsWith(cs, "supported")) cs else "insufficient_evidence"
    data.frame(edge_id = id, reference_support = e$reference_support[k], reference_scope = unname(NRV_EDGE_LAYER[id]),
               mapping_status = ms, verdict = verdict,
               reason_code = if (id == "E10") "bm_to_pb_out_of_scope" else if (circ) "anchor_circular" else if (!nrow(pr)) "holdout_unavailable" else if (ms != "single") paste0("mapping_", ms) else "evaluated",
               program_match = if (nrow(pr)) paste(pr$program_id, pr$status, sep = ":", collapse = ";") else "",
               candidate_stage_set = if (!is.null(mp) && nrow(mp)) mp$subroute[1] else "",
               n_observed_edges = if (!is.null(mp) && nrow(mp)) as.integer(mp$n_observed_edges[1]) else NA_integer_,
               exact_stage_mapping = "not_confirmed",
               observed_reverse_trend = identical(cs, "conflicting") || any(pr$status == "conflicting"),   # 대응이 없어도 반대 방향 프로그램은 표시 (Codex 검토)
               stringsAsFactors = FALSE)
  }))
}

#' 경로 주장은 한 키(분할, holdout, fold, 배정) 안에서만 한다
nrv_assert_single_key <- function(results) {
  k <- unique(results[, c("partition", "holdout_scheme", "fold", "allocation")])
  if (nrow(k) != 1L) stop("경로 주장은 한 키(분할, holdout, fold, 배정) 안에서만 한다. 여러 키가 섞였다.")
  invisible(TRUE)
}

#' 세 조건을 같은 분모(10 문헌 edge)로 쌓는다. marker_only 의 연결 근거는 not_evaluated.
nrv_condition_table <- function(key, marker_only, route, combined) {
  ids <- sprintf("E%02d", 1:10)
  for (d in list(marker_only, route, combined)) if (!setequal(d$edge_id, ids) || nrow(d) != 10L) stop("세 조건의 분모(10 문헌 edge)가 같지 않다.")
  mk <- function(d, cond, conn) data.frame(key, condition = cond, edge_id = d$edge_id, status = d$status, connection_evidence = conn, stringsAsFactors = FALSE)
  rbind(mk(marker_only, "marker_only", "not_evaluated"), mk(route, "anchor_route", "evaluated"), mk(combined, "combined", "evaluated"))
}

#' 전환 집계: 같은 가설에서 비교 가능한 상태만 짝짓는다. marker 의 not_evaluated 를 이득으로 빼지 않는다.
nrv_transitions <- function(cond_marker, cond_combined) {
  m <- merge(cond_marker[, c("edge_id", "status")], cond_combined[, c("edge_id", "status")], by = "edge_id", suffixes = c("_marker", "_combined"))
  sup <- startsWith(m$status_combined, "supported")
  kind <- ifelse(m$status_marker == "marker_pair_present" & sup, "confirmed_by_route",
          ifelse(m$status_marker == "marker_pair_present" & m$status_combined == "conflicting", "rejected_by_route",
          ifelse(m$status_marker != "marker_pair_present" & sup, "route_support_without_marker_pair",
          ifelse(m$status_marker == "marker_pair_present", "kept_unresolved", "not_comparable"))))
  as.data.frame(table(kind = factor(kind, levels = c("confirmed_by_route", "rejected_by_route", "route_support_without_marker_pair",
                                                     "kept_unresolved", "not_comparable"))))
}

# ── 세 조건과 대응 (6.4) ────────────────────────────────────────────────────
#' 선택에 쓰인 유전자 (배정 x 조건). 같은 조건의 판정에서 전역 제외한다.
nrv_selection_genes <- function(ref, allocation, condition = c("anchor_route", "combined", "marker_only_full")) {
  condition <- match.arg(condition); an <- unique(unlist(nrv_anchor_genes(allocation)))
  ro <- nrv_gene_roles(ref, allocation); ev <- unique(ro$gene[ro$role == "evaluation"])
  switch(condition, anchor_route = an, combined = unique(c(an, setdiff(ref$markers$gene, ev))), marker_only_full = unique(c(an, ref$markers$gene)))
}

#' masked node panel = node marker 중 그 배정의 eval 유전자가 아닌 것. 비면 held_masked_panel_empty. 같은 패널인 node 는 구별 불가로 표시.
nrv_masked_panels <- function(ref, allocation) {
  ro <- nrv_gene_roles(ref, allocation); ev <- unique(ro$gene[ro$role == "evaluation"]); sp <- split(ref$markers$gene, ref$markers$node_id)
  mk <- vapply(sp, function(g) paste(sort(setdiff(g, ev)), collapse = "|"), character(1))
  same <- vapply(names(mk), function(n) paste(setdiff(names(mk)[mk == mk[[n]] & nzchar(mk)], n), collapse = "|"), character(1))
  data.frame(node_id = names(sp), full_panel = vapply(sp, paste, character(1), collapse = "|"), masked_panel = unname(mk),
             status = ifelse(nzchar(mk), "ok", "held_masked_panel_empty"), indistinguishable_with = unname(same), row.names = NULL, stringsAsFactors = FALSE)
}

#' marker-only 라벨: node 점수의 equal-block 평균으로 0 점 unassigned 포함 승자, 동률은 후보 집합 그대로, 공통 block 추출 빈도
#' 결측 (Codex 재현 수정): 점수 NA 를 0 점으로 바꾸지 않는다. node 마다 "세포 ≥ min_cells 이고 점수가 유한한 block" 만 분모에 넣고,
#'   그런 block 이 없으면 그 node 는 그 군집에서 평가 불가(NA)다. unavailable_nodes 에 남긴다(예: 배정 A 의 빈 MPP masked panel).
nrv_marker_only <- function(score_node_cb, n_cb, R = 2000L, seed = 20261013L, min_cells = 20) {
  nodes <- names(score_node_cb); cl <- rownames(n_cb); big <- n_cb >= min_cells
  eb <- function(w) { B <- vapply(nodes, function(nd) { M <- score_node_cb[[nd]]; el <- big & is.finite(M); M[!el] <- 0
                                                        den <- as.numeric(el %*% w); ifelse(den > 0, as.numeric(M %*% w) / den, NA_real_) }, numeric(length(cl)))
                      matrix(B, nrow = length(cl), dimnames = list(cl, nodes)) }
  winners <- function(row) { if (all(is.na(row))) return("unavailable"); m <- mdu_winner_mass(setNames(row, nodes)); names(m)[!is.na(m) & m > 0] }
  base <- eb(rep(1, ncol(n_cb))); cs <- apply(base, 1, function(r) paste(winners(r), collapse = "|"))
  freq <- .nrv_with_seed(seed, rowMeans(replicate(R, {
    w <- tabulate(sample.int(ncol(n_cb), ncol(n_cb), replace = TRUE), ncol(n_cb))
    apply(eb(w), 1, function(r) paste(winners(r), collapse = "|")) == cs })))
  srt <- t(apply(cbind(base, unassigned = 0), 1, sort, decreasing = TRUE, na.last = TRUE))
  lab <- ifelse(grepl("|", cs, fixed = TRUE), "tie", cs)
  data.frame(cluster = cl, label = lab, candidate_set = unname(cs), margin = srt[, 1] - srt[, 2], label_freq = as.numeric(freq),
             status = ifelse(lab == "tie", "tie", ifelse(lab %in% c("unassigned", "unavailable"), lab, "labelled")),
             unavailable_nodes = unname(apply(base, 1, function(r) paste(nodes[is.na(r)], collapse = "|"))), stringsAsFactors = FALSE)
}

#' marker-only 문헌 edge 가설: 관측 쌍 중 u ∈ 후보(a), v ∈ 후보(b). 관측 그래프의 쌍 목록에 조건부라는 사실을 열로 공개한다.
nrv_marker_edge_hypotheses <- function(pt, cand, ref_edges) {
  has <- function(cl, x) vapply(cl, function(k) x %in% cand[[k]], logical(1))
  allc <- unique(unlist(cand))
  do.call(rbind, lapply(seq_len(nrow(ref_edges)), function(k) {
    u <- ref_edges$from[k]; v <- ref_edges$to[k]
    hit <- which((has(pt$cluster_a, u) & has(pt$cluster_b, v)) | (has(pt$cluster_a, v) & has(pt$cluster_b, u)))
    st <- if (!(u %in% allc) || !(v %in% allc)) "labels_unassigned" else if (length(hit)) "marker_pair_present" else "marker_pair_absent"
    data.frame(edge_id = ref_edges$edge_id[k], status = st, pair_ids = paste(pt$cluster_a[hit], pt$cluster_b[hit], sep = "|", collapse = ";"),
               candidate_universe = "observed_cluster_pairs", connection_evidence = "not_evaluated", stringsAsFactors = FALSE)
  }))
}

#' 결합 대응: 이미 고정한 경로 안의 최소 연속 부분경로에 문헌 edge 를 대응한다. 한 문헌 edge 가 여러 관측 edge 에 걸칠 수 있다.
#' bundled 는 한 군집의 동률 후보 집합에 u, v 가 모두 있을 때만이다 (없는 상태를 만들지 않는다).
nrv_combined_mapping <- function(route, cand, ref_edges, phase, internal) {
  if (is.null(internal)) internal <- character(0)
  inC <- function(x) vapply(route, function(k) x %in% cand[[k]], logical(1))
  do.call(rbind, lapply(seq_len(nrow(ref_edges)), function(k) {
    id <- ref_edges$edge_id[k]; u <- ref_edges$from[k]; v <- ref_edges$to[k]; iu <- which(inC(u)); iv <- which(inC(v))
    out <- function(ms, sub = "", ne = NA_integer_, cs = "insufficient_evidence", rs = "") data.frame(edge_id = id, mapping_status = ms, subroute = sub,
      n_observed_edges = ne, combined_status = cs, reason = rs, stringsAsFactors = FALSE)
    both <- intersect(iu, iv)
    if (length(both)) return(out("bundled", paste(route[both], collapse = "|"), 0L, rs = "tied_candidate_set"))
    spans <- do.call(rbind, lapply(iu, function(i) { j <- iv[iv > i]; if (!length(j)) return(NULL); j <- min(j)
      if (any(c(iu, iv) > i & c(iu, iv) < j)) return(NULL); data.frame(i = i, j = j) }))
    if (is.null(spans)) return(out("unmapped"))
    if (nrow(spans) > 1L) return(out("multiple", paste(vapply(seq_len(nrow(spans)), function(r) paste(route[spans$i[r]:spans$j[r]], collapse = ">"), ""), collapse = "|")))
    i <- spans$i; j <- spans$j; sub <- paste(route[i:j], collapse = ">"); mid <- route[setdiff(i:j, c(i, j))]
    ph <- phase[phase$edge_id == id, ]; inside <- !is.na(ph$segment_index) & ph$segment_index >= i & ph$segment_index <= j - 1L
    ist <- if (length(mid)) unname(internal[mid]) else character(0)
    cs <- if (!nrow(ph)) "insufficient_evidence" else if (any(ph$status[inside] == "phase_reversed") || any(ist == "internally_reversed", na.rm = TRUE)) "conflicting" else
      if (all(inside) && all(ph$status %in% c("phase_supported", "phase_supported_single_gene")) && (!length(mid) || isTRUE(all(ist == "internally_continuous"))))
        paste0(if (all(ph$conditional_only)) "supported_conditional" else "supported", if (any(ph$status == "phase_supported_single_gene")) "_single_gene_component" else "") else "insufficient_evidence"
    out("single", sub, as.integer(j - i), cs, if (cs == "insufficient_evidence") "programs_or_internal_not_all_supported" else "")
  }))
}

#' gene null 꼬리 비율. null 이 모두 NA 면 p = 1 이 아니라 NA + 사유.
nrv_gene_null_p <- function(obs, null) {
  v <- null[!is.na(null)]
  if (is.na(obs) || !length(v)) return(structure(NA_real_, reason = if (is.na(obs)) "observed_missing" else "null_all_missing"))
  (1 + sum(v >= obs)) / (1 + length(v))
}

# ── 실행기에서 관측과 gene null 이 같이 쓰는 통계 ─────────────────────────────
# gene null 은 아래 함수를 관측과 똑같이 다시 부른다. 그래서 평가 가능성 필터, 구간 σ 선택, 내부 프로그램 선택이 같은 규칙으로 반복된다.

#' 선택 규칙 하나: 문헌 방향 값(direction x 값)이 가장 큰 이름. 동률이면 순서상 앞의 것(고정 경로 순서, 프로그램 목록 순서). 값이 없으면 NA.
nrv_pick_max <- function(x, direction = 1) {
  v <- direction * x
  if (!length(v) || all(is.na(v))) return(NA_character_)
  names(x)[which.max(v)]
}

#' 프로그램 하나의 경로 통계
#'   arr = list(delta = 구간 x block x 유전자 경계 대비 Δμ, det = 같은 모양의 경계 합집합 검출 세포 수)
#'   평가 가능 유전자: 학습 block 의 검출 세포 합 ≥ min_det. 평가 가능 유전자가 min_genes 이상인 구간만 쓴다.
#'   구간 값 = 평가 가능 유전자의 부호 평균 Δ (block 마다) → 학습 block 의 equal-block 평균으로 σ 를 고른다 → 보류 block 의 Δ_P(σ)
nrv_route_program_stats <- function(arr, genes, signs, train, held, min_genes = 2, min_det = 5) {
  seg <- dimnames(arr$delta)[[1]]; blocks <- dimnames(arr$delta)[[2]]; k <- length(genes)
  det_tr <- matrix(apply(arr$det[, train, genes, drop = FALSE], c(1, 3), sum), length(seg), k, dimnames = list(seg, genes))
  ok <- det_tr >= min_det; n_ok <- setNames(as.integer(rowSums(ok)), seg)
  dP <- matrix(NA_real_, length(seg), length(blocks), dimnames = list(seg, blocks))
  for (s in seg[n_ok >= min_genes & n_ok > 0]) {
    D <- matrix(arr$delta[s, , genes, drop = FALSE], length(blocks), k)
    dP[s, ] <- as.numeric(D[, ok[s, ], drop = FALSE] %*% (signs[ok[s, ]] / n_ok[[s]]))
  }
  tr <- dP[, train, drop = FALSE]; has <- rowSums(!is.na(tr)) > 0
  train_seg <- setNames(rep(NA_real_, length(seg)), seg); train_seg[has] <- rowMeans(tr[has, , drop = FALSE], na.rm = TRUE)
  sigma <- nrv_pick_max(train_seg)
  ho <- if (is.na(sigma)) setNames(rep(NA_real_, length(held)), held) else dP[sigma, held]
  list(n_genes = n_ok, train_seg = train_seg, dP = dP, sigma = sigma, heldout = ho, estimate = if (all(is.na(ho))) NA_real_ else mean(ho, na.rm = TRUE))
}

#' 집합 합으로 구한 block 별 OLS 기울기 = cov(t, s) / var(t)  (nrv_internal_slope 의 세포 수준 계산과 같다)
#'   n = 세포 수, t1 = Σt, t2 = Σt², s1 = Σs, st = Σ t·s.  n < min_cells 이거나 t 의 분산이 수치적 0 이면 NA.
nrv_slope_from_sums <- function(n, t1, t2, s1, st, min_cells = 5) {
  vt <- t2 - t1^2 / n; out <- (st - t1 * s1 / n) / vt
  out[n < min_cells | is.na(vt) | vt <= MDU_NUM_TOL] <- NA_real_
  out
}

#' 내부 군집 B 의 프로그램 통계: 프로그램마다 block 별 기울기 → 학습 block 의 equal-block 평균이 가장 큰 프로그램 → 보류 block 기울기
#'   isum = list(n, t1, t2 (block 이름 벡터), S1 = Σx, St = Σ t·x (block x 유전자)), progs = list(이름 = list(genes, signs, min_genes))
#'   유전자가 min_genes 보다 적은 프로그램은 후보에서 뺀다.
nrv_internal_program_stats <- function(isum, progs, train, held, min_cells = 5) {
  blocks <- names(isum$n)
  progs <- progs[vapply(progs, function(p) length(p$genes) >= p$min_genes && length(p$genes) > 0, logical(1))]
  sl <- matrix(NA_real_, length(progs), length(blocks), dimnames = list(names(progs), blocks))
  for (p in names(progs)) {
    w <- progs[[p]]$signs / length(progs[[p]]$genes); g <- progs[[p]]$genes
    sl[p, ] <- nrv_slope_from_sums(isum$n, isum$t1, isum$t2, as.numeric(isum$S1[blocks, g, drop = FALSE] %*% w),
                                   as.numeric(isum$St[blocks, g, drop = FALSE] %*% w), min_cells)
  }
  tr <- sl[, train, drop = FALSE]; has <- rowSums(!is.na(tr)) > 0
  train_stat <- setNames(rep(NA_real_, nrow(sl)), rownames(sl)); train_stat[has] <- rowMeans(tr[has, , drop = FALSE], na.rm = TRUE)
  sel <- nrv_pick_max(train_stat)
  ho <- if (is.na(sel)) setNames(rep(NA_real_, length(held)), held) else sl[sel, held]
  list(slopes = sl, train_stat = train_stat, selected = sel, heldout = ho, estimate = if (all(is.na(ho))) NA_real_ else mean(ho, na.rm = TRUE))
}
