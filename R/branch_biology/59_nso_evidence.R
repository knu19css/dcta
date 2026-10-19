# =============================================================================
# 59_nso_evidence.R — 호중구 성숙 단계 관측(NSO): 실제 RNA 근거 · 관측 연결 요약 (계산만, 원저자 라벨 · UMAP 없음)
# =============================================================================
#
# 계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md (2026-10-19), 계획 Task 2. Ruling: claude_progress.md (7, 9, 10)
#   - 원 RNA 를 chunk 로 한 번 읽어(R25 nrv_log_chunk: x = log1p(10^4 · count / 전체 유전자 library)) 참고 유전자의 세포별 값과 군집 x block 합을 모은다.
#   - 프로그램 = 단계 marker 6 개 + 단계 사이 변화(E05–E09) core / conditional x up / down. 값은 관측 유전자 logRNA 의 동일 가중 평균(서술값)이다.
#     입력에 유전자가 없으면 unavailable (0 점으로 바꾸지 않는다). 입력에 있으나 count 가 모두 0 이면 observed_zero (생물학적 부재로 읽지 않는다).
#   - 관측 연결: 원 그래프의 양수 edge 만 쓴다(R24). 군집 쌍 표 · 경계 세포 · 같은 block 경계 대비(R24 집합, R26 nrv_boundary_contrast).
#     대비 = 저장 순서의 cluster_b 경계 평균 − cluster_a 경계 평균. 반대 방향은 부호 반전일 뿐 시간 · 분화 방향 검증이 아니다.
#   - B 내부 연속성: 군집 내부 edge 만의 연결 성분으로 R24 nrv_internal_coordinate 와 같은 상태를 모든 관측 삼중쌍에 판정한다 (hop 거리 없음).
#   - 통계적 유의성 · 부호 cutoff · 단계 라벨 · 새 경로 · 새 edge 를 만들지 않는다. 평균이 높은 프로그램을 단계 정답으로 바꾸지 않는다.
#
# 주요 함수
#   nso_gene_evidence(each_chunk, gene_map, groups)              세포별 logRNA · 검출, 군집 x block n · 평균 · 검출률
#   nso_program_evidence(gene_evidence, reference)                프로그램 정의 · 군집 x block · 군집 · 세포 단계 marker 서술표
#   nso_observed_connections(edges, cells, gene_evidence)         관측 군집 쌍 · 경계 RNA 대비 (유전자 · 프로그램)
#   nso_internal_continuity(edges, cluster)                       B 내부 연결 성분과 삼중쌍 상태
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph) })

NSO_MIN_CELLS <- as.numeric(NRV_THRESHOLDS$value[NRV_THRESHOLDS$name == "set_min_cells"])   # R23 운영 기준 (block 한쪽 세포 수 하한)
NSO_ORIENTATION <- paste("delta = mean over boundary cells of cluster_b (cells with a real edge to cluster_a) minus the same for cluster_a;",
                         "stored order cluster_a < cluster_b by radix name. The reverse orientation is the sign flip only:",
                         "it is not a time, velocity or differentiation-direction test.")


# ── 1. 세포 · 군집 x block 유전자 근거 ────────────────────────────────────────
#' groups = 이름 붙은 list (예: list(p92 = data.frame(cluster, block), p51 = ...)), 세포 순서 = raw 행 순서.
#' keep_cells = 원 count 를 그대로 남길 세포 번호 (실제 RNA spot-check 용; 같은 symbol 의 여러 열은 R25 와 같이 합친다)
#' 반환: genes(입력 여부), cell(library, log, det: 세포 x 입력 유전자 희소행렬), groups(긴 표), group_mats(행렬 사본), kept(cell, counts, library), n_cells, scope
nso_gene_evidence <- function(each_chunk, gene_map, groups, keep_cells = integer(0)) {
  if (!length(groups) || is.null(names(groups)) || any(!nzchar(names(groups)))) stop("groups 는 이름 붙은 list 여야 한다 (예: list(p92 = ..., p51 = ...)).")
  n <- nrow(groups[[1]])
  for (p in names(groups)) { g <- groups[[p]]
    if (!all(c("cluster", "block") %in% names(g))) stop("groups$", p, " 에 cluster · block 열이 없다.")
    if (nrow(g) != n) stop("groups 사이 세포 수가 다르다: ", p)
    if (anyNA(g$cluster) || anyNA(g$block)) stop("groups$", p, " 에 결측이 있다.") }
  genes <- names(gene_map); have <- genes[lengths(gene_map) > 0]
  H <- lapply(groups, function(g) { k <- paste(g$cluster, g$block, sep = "\r"); uk <- sort(unique(k), method = "radix")
    Matrix::sparseMatrix(i = match(k, uk), j = seq_along(k), x = 1, dims = c(length(uk), length(k)), dimnames = list(uk, NULL)) })
  S <- lapply(H, function(h) matrix(0, nrow(h), length(have), dimnames = list(rownames(h), have))); D <- S
  lib <- rep(NA_real_, n); parts <- list(); seen <- 0L
  keep <- sort(unique(as.integer(keep_cells))); if (anyNA(keep) || any(keep < 1L | keep > n)) stop("keep_cells 가 세포 범위(1..", n, ")를 벗어났다.")
  kept <- list()
  each_chunk(function(X, rows) {
    if (max(rows) > n) stop("raw 세포 수가 groups 세포 수보다 많다.")
    lc <- nrv_log_chunk(X, gene_map); lib[rows] <<- lc$library; seen <<- seen + length(rows)
    parts[[length(parts) + 1L]] <<- list(rows = rows, log = lc$log, det = lc$det)
    for (p in names(H)) { Hc <- H[[p]][, rows, drop = FALSE]; S[[p]] <<- S[[p]] + as.matrix(Hc %*% lc$log); D[[p]] <<- D[[p]] + as.matrix(Hc %*% lc$det) }
    hit <- which(rows %in% keep)
    if (length(hit)) kept[[length(kept) + 1L]] <<- list(cell = rows[hit], library = lc$library[hit],
                                                         counts = as.matrix(methods::as(X[hit, , drop = FALSE], "CsparseMatrix") %*% .nrv_P(gene_map, ncol(X))))
  })
  if (seen != n || anyNA(lib)) stop("raw 세포 수가 groups 세포 수와 다르다: ", seen, " vs ", n)
  kp <- if (length(kept)) { cnt <- do.call(rbind, lapply(kept, `[[`, "counts")); colnames(cnt) <- have
    list(cell = as.integer(unlist(lapply(kept, `[[`, "cell"))), counts = cnt, library = unlist(lapply(kept, `[[`, "library"))) } else list(cell = integer(0), counts = NULL, library = numeric(0))
  if (!identical(as.integer(unlist(lapply(parts, `[[`, "rows"))), seq_len(n))) stop("chunk 세포 순서가 1..n 이 아니다.")
  Lm <- do.call(rbind, lapply(parts, `[[`, "log")); Dm <- do.call(rbind, lapply(parts, `[[`, "det")); colnames(Lm) <- colnames(Dm) <- have
  mats <- lapply(names(H), function(p) { sp <- strsplit(rownames(H[[p]]), "\r", fixed = TRUE)
    list(cluster = vapply(sp, `[`, "", 1L), block = vapply(sp, `[`, "", 2L), n = as.numeric(Matrix::rowSums(H[[p]])), sum = S[[p]], det = D[[p]]) }); names(mats) <- names(H)
  grp <- do.call(rbind, lapply(names(mats), function(p) { m <- mats[[p]]; k <- length(have); if (!k) return(NULL)
    data.frame(partition = p, cluster = rep(m$cluster, k), block = rep(m$block, k), gene = rep(have, each = length(m$n)), n_cells = rep(as.integer(m$n), k),
               sum_log = as.numeric(m$sum), mean_log = as.numeric(m$sum) / rep(m$n, k), n_detected = as.numeric(m$det), det_rate = as.numeric(m$det) / rep(m$n, k),
               stringsAsFactors = FALSE) }))
  if (is.null(grp)) grp <- data.frame(partition = character(0), cluster = character(0), block = character(0), gene = character(0), n_cells = integer(0), sum_log = numeric(0),
                                      mean_log = numeric(0), n_detected = numeric(0), det_rate = numeric(0), stringsAsFactors = FALSE)
  list(genes = data.frame(gene = genes, in_input = lengths(gene_map) > 0, raw_cols = vapply(gene_map, function(x) paste(x, collapse = "|"), ""),
                          status = ifelse(lengths(gene_map) > 0, "in_input", "missing_gene"), row.names = NULL, stringsAsFactors = FALSE),
       cell = list(library = lib, log = Lm, det = Dm), groups = grp, group_mats = mats, kept = kp, n_cells = as.integer(n),
       scope = sprintf("reference genes (%d listed, %d in input) x all %d cells; x = log1p(1e4 * count / all-gene library) via R25 nrv_log_chunk; group = cluster x block",
                       length(genes), length(have), n))
}


# ── 2. 프로그램 근거 ─────────────────────────────────────────────────────────
#' 프로그램 정의: 단계 marker 6 개(marker_<단계>) + 단계 사이 변화(<edge>_<core|conditional>_<up|down>, R23 program_id 와 같은 형식)
.nso_programs <- function(reference, have) {
  mk <- reference$stage_markers; ch <- reference$changes
  a <- do.call(rbind, lapply(NSO_STAGES, function(s) { g <- mk$gene[mk$stage_id == s]
    data.frame(program_id = paste0("marker_", s), program_kind = "stage_marker", stage_id = s, edge_id = "", from = "", to = "", direction = "", evidence = "",
               genes = paste(g, collapse = "|"), shared_genes = paste(mk$gene[mk$stage_id == s & nzchar(mk$shared_with)], collapse = "|"), stringsAsFactors = FALSE) }))
  k <- unique(ch[, c("edge_id", "from", "to", "evidence", "direction")])
  b <- if (nrow(k)) do.call(rbind, lapply(seq_len(nrow(k)), function(r) { w <- ch$edge_id == k$edge_id[r] & ch$evidence == k$evidence[r] & ch$direction == k$direction[r]
    data.frame(program_id = paste(k$edge_id[r], k$evidence[r], k$direction[r], sep = "_"), program_kind = paste("transition", k$evidence[r], k$direction[r], sep = "_"),
               stage_id = "", edge_id = k$edge_id[r], from = k$from[r], to = k$to[r], direction = k$direction[r], evidence = k$evidence[r],
               genes = paste(ch$gene[w], collapse = "|"), shared_genes = "", stringsAsFactors = FALSE) })) else NULL
  d <- rbind(a, b); gl <- strsplit(d$genes, "|", fixed = TRUE)
  d$n_expected <- as.integer(lengths(gl)); d$n_observed <- as.integer(vapply(gl, function(g) sum(g %in% have), 0L))
  d$missing_genes <- vapply(gl, function(g) paste(setdiff(g, have), collapse = "|"), "")
  d$single_gene_program <- d$n_expected == 1L
  d$status <- ifelse(d$n_expected == 0L, "no_reference_genes", ifelse(d$n_observed == 0L, "no_genes_in_input", "ok"))
  rownames(d) <- NULL; d
}
.nso_prog_genes <- function(def, k, have) intersect(strsplit(def$genes[k], "|", fixed = TRUE)[[1]], have)
.nso_unavailable_status <- function(def_status) ifelse(def_status == "no_reference_genes", "unavailable_no_reference_genes", "unavailable_missing_genes")

#' 프로그램 서술표. 반환: programs(정의), groups(군집 x block), clusters(군집: pooled · equal-block), cells(단계 marker 6 개: mean_log · n_detected · unavailable)
#' min_cells: equal-block 평균에 넣을 block 의 최소 세포 수 (R23 set_min_cells)
nso_program_evidence <- function(gene_evidence, reference, min_cells = NSO_MIN_CELLS) {
  have <- colnames(gene_evidence$cell$log); def <- .nso_programs(reference, have)
  grp <- list(); clu <- list()
  for (p in names(gene_evidence$group_mats)) { m <- gene_evidence$group_mats[[p]]; ng <- length(m$n)
    ucl <- fdt_cluster_order_safe(m$cluster); ci <- match(m$cluster, ucl); ncl <- as.numeric(tapply(m$n, factor(m$cluster, levels = ucl), sum))
    nbk <- as.integer(tapply(m$n, factor(m$cluster, levels = ucl), length)); nbk_ok <- as.integer(tapply(m$n >= min_cells, factor(m$cluster, levels = ucl), sum))
    for (k in seq_len(nrow(def))) { g <- .nso_prog_genes(def, k, have)
      if (!length(g)) { st <- .nso_unavailable_status(def$status[k])
        grp[[length(grp) + 1L]] <- data.frame(partition = p, cluster = m$cluster, block = m$block, program_id = def$program_id[k], n_cells = as.integer(m$n), n_expected = def$n_expected[k],
                                              n_observed = 0L, n_detected_genes = NA_integer_, detected_genes = "", mean_log = NA_real_, mean_det_rate = NA_real_, status = st, stringsAsFactors = FALSE)
        clu[[length(clu) + 1L]] <- data.frame(partition = p, cluster = ucl, program_id = def$program_id[k], n_cells = as.integer(ncl), n_blocks = nbk, n_blocks_min_cells = nbk_ok,
                                              mean_log_pooled = NA_real_, mean_log_equal_block = NA_real_, mean_det_rate_pooled = NA_real_, status = st, stringsAsFactors = FALSE)
        next }
      mean_g <- m$sum[, g, drop = FALSE] / m$n; det_g <- m$det[, g, drop = FALSE]
      ml <- rowMeans(mean_g); dr <- rowMeans(det_g / m$n); nd <- as.integer(rowSums(det_g > 0))
      dg <- vapply(seq_len(ng), function(r) paste(g[det_g[r, ] > 0], collapse = "|"), "")
      grp[[length(grp) + 1L]] <- data.frame(partition = p, cluster = m$cluster, block = m$block, program_id = def$program_id[k], n_cells = as.integer(m$n), n_expected = def$n_expected[k],
                                            n_observed = length(g), n_detected_genes = nd, detected_genes = dg, mean_log = unname(ml), mean_det_rate = unname(dr),
                                            status = ifelse(nd == 0L, "observed_zero", "observed"), stringsAsFactors = FALSE)
      cs <- rowsum(m$sum[, g, drop = FALSE], factor(m$cluster, levels = ucl), reorder = TRUE); cd <- rowsum(det_g, factor(m$cluster, levels = ucl), reorder = TRUE)
      ok <- m$n >= min_cells; eb <- vapply(seq_along(ucl), function(c0) { w <- ci == c0 & ok; if (any(w)) mean(ml[w]) else NA_real_ }, 0)
      clu[[length(clu) + 1L]] <- data.frame(partition = p, cluster = ucl, program_id = def$program_id[k], n_cells = as.integer(ncl), n_blocks = nbk, n_blocks_min_cells = nbk_ok,
                                            mean_log_pooled = unname(rowMeans(cs / ncl)), mean_log_equal_block = eb, mean_det_rate_pooled = unname(rowMeans(cd / ncl)),
                                            status = ifelse(rowSums(cd > 0) == 0, "observed_zero", "observed"), stringsAsFactors = FALSE) } }
  # 세포: 단계 marker 6 개만 (단계 사이 변화는 세포 하나의 값이 아니라 연결 대비로 본다)
  Lc <- gene_evidence$cell$log; Dc <- gene_evidence$cell$det; n <- nrow(Lc); ids <- paste0("marker_", NSO_STAGES)
  cm <- matrix(NA_real_, n, length(ids), dimnames = list(NULL, ids)); cn <- cm
  for (k in match(ids, def$program_id)) { g <- .nso_prog_genes(def, k, have); if (!length(g)) next
    cm[, def$program_id[k]] <- as.numeric(Lc[, g, drop = FALSE] %*% rep(1 / length(g), length(g))); cn[, def$program_id[k]] <- as.numeric(Matrix::rowSums(Dc[, g, drop = FALSE])) }
  un <- def[def$program_id %in% ids & def$status != "ok", , drop = FALSE]
  unv <- if (nrow(un)) paste(un$program_id, un$status, sep = ":", collapse = ";") else "none"
  bind <- function(x) { if (!length(x)) return(NULL); out <- do.call(rbind, x); rownames(out) <- NULL; out }
  list(programs = def, groups = bind(grp), clusters = bind(clu), cells = list(mean_log = cm, n_detected = cn, unavailable = rep(unv, n)), min_cells = min_cells)
}

#' 군집 ID 순서 (숫자면 숫자 순서, 아니면 radix) — R54 fdt_cluster_order 와 같은 규칙을 여기서 쓴다 (R54 를 적재하지 않아도 되게)
fdt_cluster_order_safe <- function(ids) {
  ids <- unique(as.character(ids)); n <- suppressWarnings(as.numeric(ids))
  if (!anyNA(n)) ids[order(n)] else sort(ids, method = "radix")
}


# ── 3. 관측 연결 ─────────────────────────────────────────────────────────────
.nso_feature_template <- function() data.frame(cluster_a = character(0), cluster_b = character(0), n_cell_edges = integer(0), weight_sum = numeric(0), q = numeric(0),
                                               n_boundary_a = integer(0), n_boundary_b = integer(0), n_same_block_edges = integer(0), feature_kind = character(0),
                                               feature = character(0), mean_boundary_a = numeric(0), mean_boundary_b = numeric(0), det_rate_boundary_a = numeric(0),
                                               det_rate_boundary_b = numeric(0), delta_pooled = numeric(0), delta_equal_block = numeric(0), n_eligible_blocks = integer(0),
                                               n_blocks_with_boundary = integer(0), status = character(0), stringsAsFactors = FALSE)

#' 관측 군집 쌍의 경계 RNA 대비 (유전자 + 프로그램). cells = data.frame(cluster, block, study) (세포 순서 = gene_evidence 세포 순서).
#' 경계 집합: pool|a>b = a 의 세포 중 b 와 실제 edge 가 있는 세포 (R24 nrv_set_definitions), sb|a>b|block = 같은 block 이웃만.
#' status 는 같은 block 대비(equal-block 평균)의 가용성이다: ok(적격 block ≥ 2) · single_block · no_eligible_block · no_same_block_boundary.
#' 입력에 없는 유전자 · 유전자 없는 프로그램은 값 NA 와 missing_gene / unavailable_* 상태 (0 으로 바꾸지 않는다).
nso_observed_connections <- function(edges, cells, gene_evidence, programs = NULL, min_cells = NSO_MIN_CELLS) {
  cl <- as.character(cells$cluster); bl <- as.character(cells$block); st <- as.character(cells$study); n <- length(cl)
  Lc <- gene_evidence$cell$log; Dc <- gene_evidence$cell$det; have <- colnames(Lc); genes_all <- gene_evidence$genes$gene
  if (nrow(Lc) != n) stop("세포 수가 gene evidence 와 다르다: ", n, " vs ", nrow(Lc))
  meta <- list(orientation = NSO_ORIENTATION, min_cells = min_cells,
               boundary_definition = "boundary cells of a toward b = cells of a with at least one positive graph edge to a cell of b (R24 sets); pooled = all blocks; equal-block = mean over blocks with >= min_cells boundary cells on both sides",
               scope = "observed cluster pairs only (n_cell_edges >= 1); no new edges, no UMAP distance")
  pt <- nrv_pair_table(edges, cl, bl, st)
  if (!nrow(pt)) return(list(pairs = pt, features = .nso_feature_template(), block_deltas = list(), meta = meta))
  sd <- nrv_set_definitions(edges, data.frame(cluster = cl, block = bl, stringsAsFactors = FALSE))
  agg <- function(sets) { if (!length(sets)) return(NULL); H <- nrv_set_matrix(sets, n)
    list(sum = as.matrix(H %*% Lc), det = as.matrix(H %*% Dc), n = setNames(as.numeric(Matrix::rowSums(H)), rownames(H))) }
  P <- agg(sd$pool$sets); B <- agg(sd$sb$sets); sbm <- sd$sb$meta
  def <- programs; if (is.null(def)) def <- data.frame(program_id = character(0), genes = character(0), stringsAsFactors = FALSE)
  out <- vector("list", nrow(pt)); bdel <- vector("list", nrow(pt))
  for (r in seq_len(nrow(pt))) { a <- pt$cluster_a[r]; b <- pt$cluster_b[r]
    fa <- paste0("pool|", a, ">", b); fb <- paste0("pool|", b, ">", a)
    row <- function(M, k) stats::setNames(as.numeric(M[k, , drop = FALSE]), colnames(M))   # 유전자가 1 개여도 이름을 잃지 않게
    ma <- row(P$sum, fa) / P$n[[fa]]; mb <- row(P$sum, fb) / P$n[[fb]]; da <- row(P$det, fa) / P$n[[fa]]; db <- row(P$det, fb) / P$n[[fb]]
    bk <- if (is.null(sbm)) character(0) else sort(sbm$block[sbm$from == a & sbm$to == b], method = "radix")
    dl <- NULL
    if (length(bk)) { ia <- paste("sb", paste0(a, ">", b), bk, sep = "|"); ib <- paste("sb", paste0(b, ">", a), bk, sep = "|")
      if (!all(ib %in% names(B$n))) stop("같은 block 경계 집합이 한쪽에만 있다: ", a, "–", b)
      dl <- nrv_boundary_contrast(B$sum[ia, , drop = FALSE], B$n[ia], B$sum[ib, , drop = FALSE], B$n[ib], min_cells); rownames(dl) <- bk }
    bdel[r] <- list(dl)                                                       # NULL 이어도 자리를 지운다 (bdel[[r]] <- NULL 은 원소를 없앤다)
    blk_stat <- function(v) { ne <- sum(!is.na(v)); list(eq = if (ne) mean(v, na.rm = TRUE) else NA_real_, ne = as.integer(ne),
      status = if (!length(bk)) "no_same_block_boundary" else if (!ne) "no_eligible_block" else if (ne == 1L) "single_block" else "ok") }
    rows <- list()
    for (g in genes_all) {
      if (!g %in% have) { rows[[length(rows) + 1L]] <- list("gene", g, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_integer_, length(bk), "missing_gene"); next }
      s <- blk_stat(if (is.null(dl)) numeric(0) else dl[, g])
      rows[[length(rows) + 1L]] <- list("gene", g, ma[[g]], mb[[g]], da[[g]], db[[g]], mb[[g]] - ma[[g]], s$eq, s$ne, length(bk), s$status) }
    for (k in seq_len(nrow(def))) { g <- intersect(strsplit(def$genes[k], "|", fixed = TRUE)[[1]], have)
      if (!length(g)) { st0 <- if (!is.null(def$status) && identical(def$status[k], "no_reference_genes")) "unavailable_no_reference_genes" else "unavailable_missing_genes"
        rows[[length(rows) + 1L]] <- list("program", def$program_id[k], NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_integer_, length(bk), st0); next }
      s <- blk_stat(if (is.null(dl)) numeric(0) else rowMeans(dl[, g, drop = FALSE]))
      rows[[length(rows) + 1L]] <- list("program", def$program_id[k], mean(ma[g]), mean(mb[g]), mean(da[g]), mean(db[g]), mean(mb[g] - ma[g]), s$eq, s$ne, length(bk), s$status) }
    v <- function(j) vapply(rows, function(x) x[[j]], if (j %in% c(1, 2, 11)) "" else if (j %in% c(9, 10)) 0L else 0)
    out[[r]] <- data.frame(cluster_a = a, cluster_b = b, n_cell_edges = as.integer(pt$n_cell_edges[r]), weight_sum = pt$weight_sum[r], q = pt$q[r],
                           n_boundary_a = as.integer(pt$n_boundary_a[r]), n_boundary_b = as.integer(pt$n_boundary_b[r]), n_same_block_edges = as.integer(pt$n_same_block_edges[r]),
                           feature_kind = v(1), feature = v(2), mean_boundary_a = v(3), mean_boundary_b = v(4), det_rate_boundary_a = v(5), det_rate_boundary_b = v(6),
                           delta_pooled = v(7), delta_equal_block = v(8), n_eligible_blocks = v(9), n_blocks_with_boundary = v(10), status = v(11), stringsAsFactors = FALSE) }
  names(bdel) <- paste(pt$cluster_a, pt$cluster_b, sep = "|")
  f <- do.call(rbind, out); rownames(f) <- NULL
  list(pairs = pt, features = f, block_deltas = bdel, meta = meta)
}


# ── 4. B 내부 연속성 ─────────────────────────────────────────────────────────
#' 군집 내부 edge 만으로 만든 연결 성분(군집을 넘지 않는다)과, 관측 삼중쌍 A–B–C (A, C = B 의 서로 다른 이웃 군집, A < C) 의 상태.
#' 상태 정의는 R24 nrv_internal_coordinate 와 같다: no_boundary (E_A 또는 E_C 가 빔) · no_internal_span (한 경계가 다른 경계에 포함)
#'   · internal_disconnected (E_A 와 E_C 가 같은 성분에 없다) · span_ok. hop 거리 · 위치 t 는 계산하지 않는다.
nso_internal_continuity <- function(edges, cluster) {
  cluster <- as.character(cluster); n <- length(cluster); ci <- cluster[edges$i]; cj <- cluster[edges$j]; inside <- ci == cj
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (any(inside)) g <- igraph::add_edges(g, as.vector(rbind(edges$i[inside], edges$j[inside])))
  memb <- igraph::components(g)$membership
  x <- c(edges$i[!inside], edges$j[!inside]); y <- c(cj[!inside], ci[!inside])          # 세포 x 가 군집 y 와 실제 edge 를 가진다
  bnd <- lapply(split(x, paste(cluster[x], y, sep = "\r")), function(v) sort(unique(v)))
  ucl <- fdt_cluster_order_safe(cluster); tri <- list()
  summ <- do.call(rbind, lapply(ucl, function(B) {
    keys <- names(bnd)[startsWith(names(bnd), paste0(B, "\r"))]; nb <- sort(substring(keys, nchar(B) + 2L), method = "radix")
    mB <- memb[cluster == B]; tb <- table(mB)
    if (length(nb) >= 2L) { cmb <- utils::combn(length(nb), 2L)
      st <- vapply(seq_len(ncol(cmb)), function(k) { EA <- bnd[[paste0(B, "\r", nb[cmb[1, k]])]]; EC <- bnd[[paste0(B, "\r", nb[cmb[2, k]])]]
        if (!length(EA) || !length(EC)) "no_boundary" else if (all(EA %in% EC) || all(EC %in% EA)) "no_internal_span" else
          if (!length(intersect(memb[EA], memb[EC]))) "internal_disconnected" else "span_ok" }, "")
      cnt <- function(k, f) vapply(seq_len(ncol(cmb)), function(j) f(bnd[[paste0(B, "\r", nb[cmb[k, j]])]], j), 0)
      tri[[length(tri) + 1L]] <<- data.frame(cluster_B = B, cluster_A = nb[cmb[1, ]], cluster_C = nb[cmb[2, ]],
        n_EA = as.integer(cnt(1, function(E, j) length(E))), n_EC = as.integer(cnt(2, function(E, j) length(E))),
        n_shared_cells = as.integer(vapply(seq_len(ncol(cmb)), function(j) length(intersect(bnd[[paste0(B, "\r", nb[cmb[1, j]])]], bnd[[paste0(B, "\r", nb[cmb[2, j]])]])), 0L)),
        n_components_EA = as.integer(cnt(1, function(E, j) length(unique(memb[E])))), n_components_EC = as.integer(cnt(2, function(E, j) length(unique(memb[E])))),
        status = st, stringsAsFactors = FALSE)
    } else st <- character(0)
    data.frame(cluster = B, n_cells = sum(cluster == B), n_internal_edges = sum(inside & ci == B), n_components = length(tb), largest_component_cells = as.integer(max(tb)),
               n_singleton_cells = as.integer(sum(tb == 1L)), n_neighbors = length(nb), n_triples = length(st), n_span_ok = sum(st == "span_ok"),
               n_internal_disconnected = sum(st == "internal_disconnected"), n_no_internal_span = sum(st == "no_internal_span"), stringsAsFactors = FALSE) }))
  tr <- if (length(tri)) do.call(rbind, tri) else data.frame(cluster_B = character(0), cluster_A = character(0), cluster_C = character(0), n_EA = integer(0), n_EC = integer(0),
                                                             n_shared_cells = integer(0), n_components_EA = integer(0), n_components_EC = integer(0), status = character(0),
                                                             stringsAsFactors = FALSE)
  rownames(summ) <- NULL; rownames(tr) <- NULL
  list(clusters = summ, triples = tr,
       definition = "components of the within-cluster subgraph (cross-cluster edges removed); status as R24 nrv_internal_coordinate: no_boundary / no_internal_span / internal_disconnected / span_ok")
}
