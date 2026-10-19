# =============================================================================
# 24_nrv_structure.R — 고정 그래프의 구조: donor block, 양수 edge, 군집 쌍 표, 경계 집합, fold, 경로, 내부 좌표
# =============================================================================
#
# 원칙
#   - 원그래프의 양의 edge 만 쓴다. 새 edge, 가짜 세포를 만들지 않는다.
#   - 군집 ID 는 문자열 이름으로 다룬다("0", "2", "10" 을 숫자 index 로 쓰지 않는다).
#   - 연결 강도 q 는 상대 연결 비중이다. 분화 확률이나 신뢰확률이 아니다.
#   - 경로와 내부 좌표는 그래프만으로 정한다(유전자를 쓰지 않는다).
#
# 식
#   w_ab = 군집 a, b 사이 세포 edge 가중치 합, s_a = a 의 군집 밖 가중치 합(그 분석 조건의 전체 BM 군집 그래프)
#   p_ab = w_ab / s_a,  q_ab = w_ab / sqrt(s_a s_b),  비용 c_ab = 0.5 (log s_a + log s_b) - log w_ab  (q = 1 이면 0)
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph) })

# 반복 표본: 같은 사람의 반복 채취를 같은 block 으로 묶는다 (정확한 ID)
NRV_REPEAT_BLOCKS <- c(Granja2019_BMMC_D1T1 = "Granja2019_BMMC_D1", Granja2019_BMMC_D1T2 = "Granja2019_BMMC_D1",
                       Oetjen2018_C = "Oetjen2018_C_block", Oetjen2018_Ck = "Oetjen2018_C_block",
                       Oetjen2018_S = "Oetjen2018_S_block", Oetjen2018_Sk = "Oetjen2018_S_block")

#' donor_id → donor block. 반복 표본 ID 가 자료에 없거나 block 수가 기대와 다르면 멈춘다.
nrv_blocks <- function(donor, expected_n = NULL, require_all = TRUE) {
  donor <- as.character(donor)
  if (anyNA(donor) || any(!nzchar(donor))) stop("donor_id 에 결측이 있다.")
  miss <- setdiff(names(NRV_REPEAT_BLOCKS), unique(donor))
  if (require_all && length(miss)) stop("반복 표본 ID 가 자료에 없다: ", paste(miss, collapse = ", "))
  b <- donor; hit <- donor %in% names(NRV_REPEAT_BLOCKS); b[hit] <- NRV_REPEAT_BLOCKS[donor[hit]]
  if (!is.null(expected_n) && length(unique(b)) != expected_n) stop("block 수가 ", expected_n, " 이 아니다: ", length(unique(b)))
  unname(b)
}

#' 상삼각(i < j)의 양수 edge 만. 저장된 0 과 대각은 제외하고 개수를 attr("dropped") 에 남긴다.
#' 음수·비유한·중복·비대칭은 오류다 (R10 ee_upper_edges 에는 양수 필터가 없어 새로 감싼다).
nrv_positive_upper_edges <- function(G) {
  if (methods::is(G, "TsparseMatrix") && anyDuplicated(paste(G@i, G@j))) stop("그래프에 중복 항목이 있다.")
  G <- methods::as(G, "CsparseMatrix")
  if (nrow(G) != ncol(G)) stop("그래프가 정사각 행렬이 아니다.")
  if (any(!is.finite(G@x))) stop("그래프에 비유한 가중치가 있다.")
  if (any(G@x < 0)) stop("그래프에 음수 가중치가 있다.")
  if (!Matrix::isSymmetric(G)) stop("그래프가 대칭이 아니다.")
  n_diag <- sum(Matrix::diag(G) != 0)
  up <- Matrix::summary(Matrix::triu(G, k = 1))
  n_zero <- sum(up$x == 0)
  e <- data.frame(i = as.integer(up$i), j = as.integer(up$j), w = as.numeric(up$x))[up$x > 0, , drop = FALSE]
  rownames(e) <- NULL; attr(e, "dropped") <- c(diagonal = n_diag, stored_zero = n_zero); e
}

.nrv_lv <- function(x) sort(unique(as.character(x)), method = "radix")       # 지역 설정과 무관한 이름 순서

#' 관측 군집 쌍 표. cell_mask 가 있으면 양끝이 모두 마스크 안인 edge 만 남기고, 분자와 분모를 그 edge 로 다시 계산한다.
nrv_pair_table <- function(e, cluster, block, study, cell_mask = NULL) {
  cluster <- as.character(cluster); block <- as.character(block); study <- as.character(study)
  if (!is.null(cell_mask)) e <- e[cell_mask[e$i] & cell_mask[e$j], , drop = FALSE]
  a <- cluster[e$i]; b <- cluster[e$j]; cr <- a != b; e <- e[cr, , drop = FALSE]; a <- a[cr]; b <- b[cr]
  lv <- .nrv_lv(cluster); sw <- match(a, lv) > match(b, lv)
  ca <- ifelse(sw, b, a); cb <- ifelse(sw, a, b); ia <- ifelse(sw, e$j, e$i); ib <- ifelse(sw, e$i, e$j)
  s <- tapply(c(e$w, e$w), c(a, b), sum)                                     # 군집 밖 가중치 합 (내부 edge 는 이미 제외)
  key <- paste(ca, cb, sep = "\r"); sp <- split(seq_len(nrow(e)), factor(key, levels = unique(key)))
  sizes <- if (is.null(cell_mask)) table(cluster) else table(cluster[cell_mask])
  out <- do.call(rbind, lapply(sp, function(ix) {
    sb <- block[ia[ix]] == block[ib[ix]]
    data.frame(cluster_a = ca[ix[1]], cluster_b = cb[ix[1]], n_cell_edges = length(ix), weight_sum = sum(e$w[ix]),
               n_same_block_edges = sum(sb), n_same_study_edges = sum(study[ia[ix]] == study[ib[ix]]),
               n_blocks_same_block = length(unique(block[ia[ix]][sb])), n_boundary_a = length(unique(ia[ix])),
               n_boundary_b = length(unique(ib[ix])), stringsAsFactors = FALSE)
  }))
  if (is.null(out)) return(data.frame(cluster_a = character(0), cluster_b = character(0), n_cell_edges = integer(0), weight_sum = numeric(0),
                                      n_same_block_edges = integer(0), n_same_study_edges = integer(0), n_blocks_same_block = integer(0),
                                      n_boundary_a = integer(0), n_boundary_b = integer(0), size_a = integer(0), size_b = integer(0),
                                      s_a = numeric(0), s_b = numeric(0), p_ab = numeric(0), p_ba = numeric(0), q = numeric(0), cost = numeric(0)))
  rownames(out) <- NULL
  out$size_a <- as.integer(sizes[out$cluster_a]); out$size_b <- as.integer(sizes[out$cluster_b])
  out$s_a <- as.numeric(s[out$cluster_a]); out$s_b <- as.numeric(s[out$cluster_b])
  out$p_ab <- out$weight_sum / out$s_a; out$p_ba <- out$weight_sum / out$s_b
  out$q <- out$weight_sum / sqrt(out$s_a * out$s_b)
  out$cost <- 0.5 * (log(out$s_a) + log(out$s_b)) - log(out$weight_sum)
  out$cost[abs(out$cost) < 1e-15] <- 0                                       # 반올림 잡음만 0 으로 정리 (q = 1)
  out
}

#' 경계 집합: from 군집 세포 중 to 군집과 실제 edge 가 있는 세포. same_block 이면 같은 block 이웃만 센다.
nrv_boundary_sets <- function(e, cluster, block, from, to, same_block = FALSE) {
  cluster <- as.character(cluster); block <- as.character(block)
  ii <- c(e$i, e$j); jj <- c(e$j, e$i)                                       # 양방향
  sets <- lapply(seq_along(from), function(k) {
    keep <- cluster[ii] == from[k] & cluster[jj] == to[k]
    if (same_block) keep <- keep & block[ii] == block[jj]
    sort(unique(ii[keep]))
  })
  setNames(sets, paste(from, to, sep = ">"))
}

#' 실행기가 한 번에 쓰는 집합 정의 (그래프만). 결과는 nrv_boundary_sets 와 같은 정의다.
#'   cb   : 군집 x block 의 세포                                   id "cb|<군집>|<block>"
#'   sb   : 같은 block 경계 B^b_{i→j} (방향 있는 관측 쌍 x block)   id "sb|<i>><j>|<block>"
#'   pool : 전체 경계 B_{i→j} (block 을 가리지 않음, 서술용)          id "pool|<i>><j>"
#' 반환: list(cb, sb, pool) = 각각 list(sets = 이름 붙은 세포 번호 목록, meta = data.frame(set_id, cluster 또는 from/to, block, n_cells))
nrv_set_definitions <- function(e, cells) {
  cl <- as.character(cells$cluster); bl <- as.character(cells$block)
  if (any(grepl("[|>]", c(cl, bl)))) stop("군집·block 이름에 집합 id 구분자(| 또는 >)가 있다.")
  grp <- function(cell, key) {                                               # (세포, 집합 키) 중복 제거 → 키마다 정렬된 세포 번호
    uk <- sort(unique(key), method = "radix"); k <- match(key, uk)
    d <- !duplicated(cell + (k - 1) * (length(cl) + 1))
    s <- split(cell[d], factor(k[d], levels = seq_along(uk))); names(s) <- uk; lapply(s, sort)
  }
  part <- function(ids, k) vapply(strsplit(ids, "|", fixed = TRUE), `[`, character(1), k)
  ft <- function(ids, k) vapply(strsplit(part(ids, 2L), ">", fixed = TRUE), `[`, character(1), k)
  cb <- grp(seq_along(cl), paste("cb", cl, bl, sep = "|"))
  x <- e$i; y <- e$j; cross <- cl[x] != cl[y]; c1 <- c(x[cross], y[cross]); c2 <- c(y[cross], x[cross])
  dir <- paste(cl[c1], cl[c2], sep = ">"); same <- bl[c1] == bl[c2]
  pool <- grp(c1, paste("pool", dir, sep = "|")); sb <- grp(c1[same], paste("sb", dir[same], bl[c1[same]], sep = "|"))
  list(cb = list(sets = cb, meta = data.frame(set_id = names(cb), cluster = part(names(cb), 2L), block = part(names(cb), 3L),
                                              n_cells = lengths(cb, use.names = FALSE), stringsAsFactors = FALSE)),
       sb = list(sets = sb, meta = data.frame(set_id = names(sb), from = ft(names(sb), 1L), to = ft(names(sb), 2L), block = part(names(sb), 3L),
                                              n_cells = lengths(sb, use.names = FALSE), stringsAsFactors = FALSE)),
       pool = list(sets = pool, meta = data.frame(set_id = names(pool), from = ft(names(pool), 1L), to = ft(names(pool), 2L),
                                                  n_cells = lengths(pool, use.names = FALSE), stringsAsFactors = FALSE)))
}

#' 집합 x 세포 지시 희소행렬 (집합 수준 합을 행렬 곱으로 구하는 데 쓴다)
nrv_set_matrix <- function(sets, n_cells) {
  Matrix::sparseMatrix(i = rep(seq_along(sets), lengths(sets)), j = unlist(sets), x = 1, dims = c(length(sets), n_cells),
                       dimnames = list(names(sets), NULL))
}

#' fold: 층(study x sorting 구분) 안에서 seed 로 섞은 뒤 돌아가며 배정한다. 각 block 은 정확히 한 fold 에서 보류된다.
nrv_folds <- function(block_meta, k = 5, seed = 20261013) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  bm <- block_meta[order(block_meta$block, method = "radix"), ]
  bm$key <- stats::runif(nrow(bm)); bm$stratum <- paste(bm$study, bm$sorting_class)
  bm <- bm[order(bm$stratum, bm$key, method = "radix"), ]
  data.frame(block = bm$block, fold = ((seq_len(nrow(bm)) - 1L) %% k) + 1L, stringsAsFactors = FALSE)
}

# ── 경로 (유전자 없음) ──────────────────────────────────────────────────────
# 비용 0 인 edge(q = 1)도 명시적 edge 로 둔다. vertex 이름 = 군집 ID 문자열.
.nrv_graph <- function(pt, extra = character(0)) {
  v <- data.frame(name = .nrv_lv(c(pt$cluster_a, pt$cluster_b, extra)))
  igraph::graph_from_data_frame(data.frame(from = pt$cluster_a, to = pt$cluster_b, cost = pt$cost), directed = FALSE, vertices = v)
}

# 동률 최단 경로 전부 (개수 상한 없음, Codex 재현 수정: 예전에는 Yen 50 개에서만 골랐다).
# s 에서 출발해 "지금까지 비용 + 이 edge + 이웃에서 t 까지 최소 비용 ≤ 최소 비용 + tol" 인 edge 만 따라간다.
# 노드는 반복하지 않으므로 비용 0 edge 의 고리에서도 끝난다. Yen 개수 상한은 후보 경로 보고에만 쓴다.
.nrv_tied_paths <- function(g, s, t, tol = MDU_TIE_TOL) {
  w <- igraph::E(g)$cost; el <- igraph::as_edgelist(g, names = TRUE)
  to_t <- igraph::distances(g, to = t, weights = w)[, 1]; D <- to_t[[s]]
  nb <- split(c(el[, 2], el[, 1]), c(el[, 1], el[, 2])); nw <- split(c(w, w), c(el[, 1], el[, 2]))
  out <- list()
  walk <- function(path, cost) {
    u <- path[length(path)]
    if (u == t) { out[[length(out) + 1L]] <<- path; return(invisible(NULL)) }
    for (k in seq_along(nb[[u]])) {
      v <- nb[[u]][k]; c2 <- cost + nw[[u]][k]
      if (!v %in% path && c2 + to_t[[v]] <= D + tol) walk(c(path, v), c2)
    }
  }
  walk(s, 0); out
}

# 주 경로: (S, T) 쌍 가운데 최소 비용과 승자 동률(1e-12) 안의 쌍들, 그 동률 최단 경로 전부에서 경로 문자열 사전순 최소
.nrv_best_path <- function(g, S, T, tol = MDU_TIE_TOL) {
  cand <- expand.grid(s = S, t = T, stringsAsFactors = FALSE)
  cand$cost <- mapply(function(s, t) igraph::distances(g, v = s, to = t, weights = igraph::E(g)$cost)[1, 1], cand$s, cand$t)
  cand <- cand[is.finite(cand$cost), , drop = FALSE]
  if (!nrow(cand)) return(list(status = "no_path", primary = NULL, pair_costs = cand))
  tc <- cand[cand$cost <= min(cand$cost) + tol, , drop = FALSE]
  tied <- unlist(lapply(seq_len(nrow(tc)), function(k) .nrv_tied_paths(g, tc$s[k], tc$t[k], tol)), recursive = FALSE)
  best <- tied[[order(vapply(tied, paste, character(1), collapse = ">"), method = "radix")[1]]]
  from <- best[1]; to <- best[length(best)]
  list(status = "route_ok", primary = best, primary_cost = cand$cost[cand$s == from & cand$t == to], from = from, to = to,
       pair_costs = cand, n_tied_primary = length(tied))
}

#' S 후보 → T 후보 최소 비용 경로. 상태: route_ok / no_path / anchor_overlap / no_anchor.
#' 개수 제한(k_alt)과 동률 규칙은 계산 예산이지 생물학적 기준이 아니다. 발현 결과를 본 뒤 다시 고르지 않는다.
nrv_route <- function(pt, S, T, k_alt = 10L) {
  if (!length(S) || !length(T)) return(list(status = "no_anchor", primary = NULL, overlap = character(0)))
  S <- as.character(S); T <- as.character(T); g <- .nrv_graph(pt, c(S, T)); ov <- intersect(S, T)
  if (length(ov)) {                                                           # 길이 0 경로는 회수 성공이 아니다
    S2 <- setdiff(S, ov); T2 <- setdiff(T, ov)
    nt <- if (length(S2) && length(T2)) .nrv_best_path(g, S2, T2) else list(status = "not_computable")
    return(list(status = "anchor_overlap", mapping_status = "unresolved_mapping", overlap = ov, primary = NULL, nontrivial = nt))
  }
  b <- .nrv_best_path(g, S, T)
  if (b$status != "route_ok") return(c(b, list(overlap = character(0))))
  nrv_check_path(pt, b$primary)
  ks <- igraph::k_shortest_paths(g, from = b$from, to = b$to, k = k_alt, weights = igraph::E(g)$cost)$vpaths
  alt <- data.frame(rank = seq_along(ks), path = vapply(ks, function(p) paste(igraph::as_ids(p), collapse = ">"), character(1)),
                    cost = vapply(ks, function(p) sum(igraph::E(g, path = p)$cost), numeric(1)), stringsAsFactors = FALSE)
  c(b, list(alternatives = alt, overlap = character(0)))
}

#' 경로의 모든 인접 쌍이 실제 관측 쌍(세포 edge ≥ 1)이고 노드가 반복되지 않는지
nrv_check_path <- function(pt, path) {
  if (anyDuplicated(path)) stop("경로에 반복 노드가 있다: ", paste(path, collapse = ">"))
  k <- paste(pt$cluster_a, pt$cluster_b, sep = "\r"); k2 <- paste(pt$cluster_b, pt$cluster_a, sep = "\r")
  for (r in seq_len(length(path) - 1L)) {
    p <- paste(path[r], path[r + 1L], sep = "\r"); hit <- c(which(k == p), which(k2 == p))
    if (!length(hit) || pt$n_cell_edges[hit[1]] < 1L) stop("경로 인접 쌍에 실제 edge 가 없다: ", path[r], "–", path[r + 1L])
  }
  invisible(TRUE)
}

#' 저장 순서(cluster_a, cluster_b) 기준 부호: 경로가 a → b 이면 +1, b → a 이면 -1
nrv_orient <- function(pt, from, to) {
  if (any(pt$cluster_a == from & pt$cluster_b == to)) return(1)
  if (any(pt$cluster_a == to & pt$cluster_b == from)) return(-1)
  stop("관측 쌍이 아니다: ", from, "–", to)
}

# ── B 내부 좌표 (그래프만) ──────────────────────────────────────────────────
#' B 의 유도 부분그래프에서 A 쪽 경계 E_A, C 쪽 경계 E_C 까지의 hop 거리와 위치 t = h_A / (h_A + h_C).
#' 공유 세포(E_A ∩ E_C)는 t = NA. 상태: span_ok / no_internal_span / internal_disconnected / no_boundary.
nrv_internal_coordinate <- function(e, cluster, B, A, C) {
  cluster <- as.character(cluster); cb <- which(cluster == B)
  ii <- c(e$i, e$j); jj <- c(e$j, e$i)
  EA <- sort(unique(ii[cluster[ii] == B & cluster[jj] == A])); EC <- sort(unique(ii[cluster[ii] == B & cluster[jj] == C]))
  if (!length(EA) || !length(EC)) return(list(status = "no_boundary", cells = NULL, n_shared = 0L, n_unreachable = NA_integer_))
  inside <- e[cluster[e$i] == B & cluster[e$j] == B, , drop = FALSE]
  loc <- function(x) match(x, cb); nb <- length(cb)
  hop_from <- function(src) {                                                # 출발 집합 하나만 super-source 로 붙인다
    g <- igraph::make_empty_graph(n = nb + 1L, directed = FALSE)
    g <- igraph::add_edges(g, c(rbind(loc(inside$i), loc(inside$j)), rbind(nb + 1L, loc(src))))
    as.numeric(igraph::distances(g, v = nb + 1L, to = seq_len(nb))[1, ]) - 1
  }
  hA <- hop_from(EA); hC <- hop_from(EC); shared <- cb %in% intersect(EA, EC)
  t <- ifelse(shared | !is.finite(hA) | !is.finite(hC), NA_real_, hA / (hA + hC))
  status <- if (all(EA %in% EC) || all(EC %in% EA)) "no_internal_span" else if (!any(is.finite(hC[loc(EA)]))) "internal_disconnected" else "span_ok"
  list(status = status, cells = data.frame(cell = cb, h_A = hA, h_C = hC, t = t, shared = shared),
       n_shared = sum(shared), n_unreachable = sum(!is.finite(hA) | !is.finite(hC)))
}

#' 보류 block 자체의 내부 구조: 같은 block 의 edge 만 남기고, 다른 block 세포는 경계에서 뺀다.
#' block 안에서 끊긴 것은 다른 block 이웃이 빠진 조건부 자료 부족일 수 있다(생물학적 부재가 아니다).
nrv_block_structure <- function(e, cluster, block, B, A, C) {
  cluster <- as.character(cluster); block <- as.character(block)
  ub <- sort(unique(block[cluster == B]), method = "radix")
  do.call(rbind, lapply(ub, function(b) {
    keep <- block[e$i] == b & block[e$j] == b
    cl_b <- ifelse(block == b, cluster, "__other_block__")
    ic <- nrv_internal_coordinate(e[keep, , drop = FALSE], cl_b, B, A, C)
    ii <- c(e$i[keep], e$j[keep]); jj <- c(e$j[keep], e$i[keep])
    data.frame(block = b, block_structure = ic$status, n_cells_B = sum(cl_b == B),
               n_EA = length(unique(ii[cl_b[ii] == B & cl_b[jj] == A])), n_EC = length(unique(ii[cl_b[ii] == B & cl_b[jj] == C])),
               stringsAsFactors = FALSE)
  }))
}
