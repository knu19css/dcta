# =============================================================================
# 70_gse253355_graph.R — GSE253355: (원 그래프 · 군집이 없을 때) 새 표현 그래프 구성과, 단계별 최고 군집 사이 실제 연결 감사
# =============================================================================
#
# 계약: codex_gse253355_instruction.md 4 · 6 항, codex_gse253355_local_addendum.md, codex_gse253355_review.md (군집 요약 경로 ≠ 실제 세포 경로).
#   새 그래프 (원 그래프 · UMAP · seurat_clusters 가 객체에 없을 때만, 단계 점수 전에, 패널 · 라벨 · DAG 를 받지 않는다):
#     raw count 의 vst 형 HVG(Seurat FindVariableFeatures vst 식을 Matrix 로 재구현, Seurat 함수 아님) → log1p(1e4 · count / library) 를 유전자마다 척도화
#     (상한 10) → 공분산 고유분해 PCA → 자기 포함 정확 kNN(유클리드) → SNN = 이웃 집합 Jaccard (1/15 미만 제거) → igraph Louvain(resolution, seed) →
#     표시 좌표 = igraph DrL 그래프 레이아웃 (dcta 에 uwot 이 없어 UMAP 이 아니다). 모두 새 산출물이며 원 그래프 · 원 UMAP 이라고 부르지 않는다.
#   연결 감사: 그래프에 실제로 있는 edge 만 센다(만들지 않는다). 군집 요약(quotient) 경로와 실제 세포 연결(연결 성분 · 최단 세포 경로)을 따로 적고,
#     군집 요약에만 경로가 있으면 cluster_path_only 로 표시한다. 같은 군집이 두 단계의 최고면 same_cluster 로 적고 자기 전이를 만들지 않는다.
#     연결은 분화를 입증하지 않고, 연결이 없다는 것만으로 생물학적 부재를 말하지 않는다.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph) })

GSX_GRAPH_PARAMS <- list(n_hvg = 2000L, n_pcs = 30L, k = 20L, prune = 1 / 15, resolution = 0.5, seed = 42L, scale_max = 10, loess_span = 0.3)

# ── 1. 새 표현 그래프 ─────────────────────────────────────────────────────────
#' vst 형 HVG (Seurat vst 식): 유전자 평균 · 분산 → log10 분산 ~ log10 평균 loess → 기대 sd 로 표준화(상한 sqrt(N)) 한 값의 분산 순위
gsx_hvg <- function(counts, n = 2000L, span = 0.3) {
  N <- ncol(counts); G <- nrow(counts); mu <- Matrix::rowMeans(counts); ri <- counts@i + 1L
  sq <- counts; sq@x <- sq@x^2; v <- (Matrix::rowSums(sq) - N * mu^2) / (N - 1); rm(sq)
  ok <- is.finite(v) & v > 0; fit <- stats::loess(log10(v[ok]) ~ log10(mu[ok]), span = span)
  ve <- rep(NA_real_, G); ve[ok] <- 10^stats::fitted(fit); sdv <- sqrt(ve); clip <- sqrt(N)
  Z <- counts; zx <- pmin((counts@x - mu[ri]) / sdv[ri], clip); zx[!is.finite(zx)] <- 0; Z@x <- zx^2; zz <- Matrix::rowSums(Z); rm(Z, zx)
  nnz <- tabulate(ri, G); vs <- (zz + (N - nnz) * (mu / sdv)^2) / (N - 1); vs[!ok] <- 0
  if (any(!is.finite(vs[ok]))) stop("vst 기대 분산(loess)을 계산하지 못한 유전자가 있다 (유전자 수 ", sum(ok), "): HVG 를 고르지 않고 멈춘다.")
  cand <- which(ok); top <- cand[order(-vs[cand], rownames(counts)[cand], method = "radix")][seq_len(min(as.integer(n), length(cand)))]   # 분산 0 유전자는 후보가 아니다
  list(genes = rownames(counts)[top], table = data.frame(gene = rownames(counts), mean = mu, variance = v, variance_expected = ve, variance_standardized = vs,
                                                         selected = seq_len(G) %in% top, row.names = NULL, stringsAsFactors = FALSE))
}

#' 척도화 PCA (공분산 고유분해, 결정적 부호: 각 축에서 절댓값이 가장 큰 적재값을 양수로)
gsx_pca <- function(counts, genes, n_pcs = 30L, scale_max = 10) {
  lib <- Matrix::colSums(counts); idx <- match(genes, rownames(counts)); if (anyNA(idx)) stop("HVG 가 counts 에 없다.")
  L <- methods::as(counts[idx, , drop = FALSE] %*% Matrix::Diagonal(x = 1e4 / lib), "CsparseMatrix"); L@x <- log1p(L@x)
  X <- as.matrix(Matrix::t(L)); rm(L); n <- nrow(X); mu <- colMeans(X); X <- sweep(X, 2, mu)
  sdv <- sqrt(colSums(X^2) / (n - 1)); keep <- sdv > 0; X <- X[, keep, drop = FALSE]; X <- sweep(X, 2, sdv[keep], "/"); X[X > scale_max] <- scale_max
  np <- min(as.integer(n_pcs), ncol(X), n - 1L); ev <- eigen(crossprod(X) / (n - 1), symmetric = TRUE)
  V <- ev$vectors[, seq_len(np), drop = FALSE]; sg <- apply(V, 2, function(z) sign(z[which.max(abs(z))])); V <- sweep(V, 2, sg, "*")
  S <- X %*% V; rownames(S) <- colnames(counts); colnames(S) <- paste0("PC", seq_len(np)); rownames(V) <- genes[keep]; colnames(V) <- colnames(S)
  list(scores = S, loadings = V, sdev = sqrt(pmax(ev$values[seq_len(np)], 0)), genes = genes[keep], n_pcs = np)
}

#' 정확 kNN (유클리드, 자기 자신이 첫 이웃, 나머지는 거리 · 번호 순). Y = 세포 x 차원
gsx_knn <- function(Y, k = 20L, chunk = 500L) {
  n <- nrow(Y); k <- min(as.integer(k), n); sq <- rowSums(Y^2); idx <- matrix(0L, n, k)
  for (s in seq(1L, n, by = chunk)) { r <- s:min(n, s + chunk - 1L)
    D <- outer(sq[r], sq, "+") - 2 * tcrossprod(Y[r, , drop = FALSE], Y)
    for (t in seq_along(r)) { d <- D[t, ]; d[r[t]] <- -Inf
      thr <- if (k < n) sort.int(d, partial = k)[k] else Inf; cand <- which(d <= thr); idx[r[t], ] <- cand[order(d[cand], cand)][seq_len(k)] } }
  idx
}

#' SNN: 자기 포함 이웃 집합의 Jaccard = 공유 수 / (2k - 공유 수), prune 미만은 0 (대칭 dgCMatrix, 대각 = 1)
gsx_snn <- function(nn, prune = 1 / 15) {
  n <- nrow(nn); k <- ncol(nn); K <- Matrix::sparseMatrix(i = rep(seq_len(n), k), j = as.vector(nn), x = 1, dims = c(n, n))
  S <- methods::as(methods::as(Matrix::tcrossprod(K), "generalMatrix"), "CsparseMatrix"); S@x <- S@x / (2 * k - S@x); S@x[S@x < prune] <- 0
  Matrix::drop0(S)
}

#' 양수 상삼각 edge (i < j, 대각 제외). 비대칭이면 두 방향 중 큰 값 (실제로 있는 방향의 edge 만, 새로 만들지 않는다)
gsx_upper_edges <- function(G) {
  T <- Matrix::summary(methods::as(G, "TsparseMatrix")); T <- T[T$x > 0 & T$i != T$j, , drop = FALSE]
  if (!nrow(T)) return(data.frame(i = integer(0), j = integer(0), w = numeric(0)))
  a <- pmin(T$i, T$j); b <- pmax(T$i, T$j); n <- nrow(G); key <- (a - 1) * n + (b - 1); o <- order(key, -T$x); keep <- !duplicated(key[o])
  e <- data.frame(i = as.integer(a[o][keep]), j = as.integer(b[o][keep]), w = as.numeric(T$x[o][keep])); e[order(e$i, e$j), , drop = FALSE]
}

#' Louvain (igraph, 가중 SNN, 대각 제외, seed 고정). 군집 번호 = 크기 순 0, 1, … (동률은 첫 세포 번호)
gsx_louvain <- function(S, resolution = 0.5, seed = 42L) {
  n <- nrow(S); e <- gsx_upper_edges(S); g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (nrow(e)) g <- igraph::add_edges(g, as.vector(rbind(e$i, e$j)), attr = list(weight = e$w))
  set.seed(seed); m <- as.integer(igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight, resolution = resolution)))
  sz <- tabulate(m); first <- vapply(seq_along(sz), function(z) min(which(m == z)), 0L); rk <- order(-sz, first)
  as.character(match(m, rk) - 1L)
}

#' 표시 좌표: igraph DrL 레이아웃 (가중 SNN, seed 고정). UMAP 이 아니고 계산에 쓰지 않는다.
gsx_layout <- function(S, seed = 42L) {
  n <- nrow(S); e <- gsx_upper_edges(S); g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (nrow(e)) g <- igraph::add_edges(g, as.vector(rbind(e$i, e$j)), attr = list(weight = e$w))
  set.seed(seed); L <- igraph::layout_with_drl(g, weights = igraph::E(g)$weight)
  colnames(L) <- c("layout_1", "layout_2"); rownames(L) <- rownames(S); L
}


# ── 2. 연결 감사 ─────────────────────────────────────────────────────────────
#' 군집 쌍 edge 표 (R24 nrv_pair_table 재사용; block = donor). 실제 edge 가 없는 쌍은 표에 없다.
gsx_cluster_edges <- function(e, cluster, donor) nrv_pair_table(e, as.character(cluster), as.character(donor), rep("all", length(cluster)))

#' 연결 성분 (R62 nca_components 재사용). keep 이 있으면 그 세포끼리의 유도 그래프
gsx_components <- function(e, n, keep = NULL) nca_components(n, e, keep = keep)

.gsx_cell_graph <- function(e, n) { g <- igraph::make_empty_graph(n = n, directed = FALSE); if (nrow(e)) g <- igraph::add_edges(g, as.vector(rbind(e$i, e$j))); g }

#' 세포 집합 A → C 의 최단 실제 세포 경로 (계산용 임시 출발점을 A 전체에 잇고 BFS; 보고하는 경로에는 그 출발점이 없다)
.gsx_cell_path <- function(g, A, C) {
  n <- igraph::vcount(g); g2 <- igraph::add_vertices(g, 1L); s <- n + 1L; g2 <- igraph::add_edges(g2, as.vector(rbind(s, A)))
  d <- as.numeric(igraph::distances(g2, v = s, to = C)); if (all(!is.finite(d))) return(list(hops = NA_integer_, cells = integer(0)))
  tgt <- C[which.min(d)]; p <- as.integer(igraph::shortest_paths(g2, from = s, to = tgt, output = "vpath")$vpath[[1]]); p <- p[p != s]
  list(hops = as.integer(min(d) - 1), cells = p)
}

#' 단계별 최고 군집 쌍마다: 같은 군집 · 직접 edge · 군집 요약 경로 · 실제 세포 연결(성분 · 최단 세포 경로). keep = 유도 그래프 세포(없으면 전체)
gsx_peak_paths <- function(e, cluster, donor, peaks, n, keep = NULL) {
  cluster <- as.character(cluster); donor <- as.character(donor); if (is.null(keep)) keep <- rep(TRUE, n)
  ek <- e[keep[e$i] & keep[e$j], , drop = FALSE]; g <- .gsx_cell_graph(ek, n); comp <- igraph::components(g)$membership
  a <- cluster[ek$i]; b <- cluster[ek$j]; x <- a != b; kc <- fdt_cluster_order_safe(unique(cluster[keep]))
  qg <- igraph::make_empty_graph(n = length(kc), directed = FALSE)
  if (any(x)) { ia <- match(a[x], kc); ib <- match(b[x], kc); key <- unique(pmin(ia, ib) * (length(kc) + 1L) + pmax(ia, ib))
    qg <- igraph::add_edges(qg, as.vector(rbind(key %/% (length(kc) + 1L), key %% (length(kc) + 1L)))) }
  st <- peaks$stage; pcs <- lapply(as.character(peaks$peak_clusters), function(z) if (nzchar(z)) strsplit(z, "|", fixed = TRUE)[[1]] else character(0)); rows <- list(); paths <- list()
  for (x1 in seq_len(length(st) - 1L)) for (x2 in (x1 + 1L):length(st)) { A <- pcs[[x1]]; C <- pcs[[x2]]
    r <- data.frame(from_stage = st[x1], to_stage = st[x2], consecutive = x2 == x1 + 1L, from_clusters = paste(A, collapse = "|"), to_clusters = paste(C, collapse = "|"),
                    from_cluster_cells = if (length(A)) sum(cluster %in% A) else NA_integer_, from_cells_in_graph = if (length(A)) sum(cluster %in% A & keep) else NA_integer_,
                    to_cluster_cells = if (length(C)) sum(cluster %in% C) else NA_integer_, to_cells_in_graph = if (length(C)) sum(cluster %in% C & keep) else NA_integer_,
                    relation = "", direct_cell_edges = NA_integer_, direct_weight = NA_real_, direct_endpoint_cells_from = NA_integer_, direct_endpoint_cells_to = NA_integer_,
                    direct_donors = "", quotient_status = "", quotient_hops = NA_integer_, quotient_path = "", cells_connected = NA, cell_hops = NA_integer_,
                    cell_path_clusters = "", path_agreement = "", stringsAsFactors = FALSE)
    if (!length(A) || !length(C)) { r$relation <- "unavailable"; rows[[length(rows) + 1L]] <- r; next }
    if (length(intersect(A, C))) { r$relation <- "same_cluster"; rows[[length(rows) + 1L]] <- r; next }
    Ac <- which(cluster %in% A & keep); Cc <- which(cluster %in% C & keep)
    if (!length(Ac) || !length(Cc)) { r$relation <- "peak_not_in_subset"; rows[[length(rows) + 1L]] <- r; next }
    r$relation <- "different_clusters"
    h <- (cluster[ek$i] %in% A & cluster[ek$j] %in% C) | (cluster[ek$i] %in% C & cluster[ek$j] %in% A)
    r$direct_cell_edges <- sum(h); r$direct_weight <- sum(ek$w[h])
    ea <- unique(c(ek$i[h][cluster[ek$i[h]] %in% A], ek$j[h][cluster[ek$j[h]] %in% A])); ec <- unique(c(ek$i[h][cluster[ek$i[h]] %in% C], ek$j[h][cluster[ek$j[h]] %in% C]))
    r$direct_endpoint_cells_from <- length(ea); r$direct_endpoint_cells_to <- length(ec); r$direct_donors <- paste(sort(unique(donor[c(ea, ec)])), collapse = "|")
    best <- NULL
    for (ca in intersect(A, kc)) for (cc in intersect(C, kc)) { p <- suppressWarnings(igraph::shortest_paths(qg, from = match(ca, kc), to = match(cc, kc), output = "vpath")$vpath[[1]])
      if (length(p) && (is.null(best) || length(p) < length(best))) best <- as.integer(p) }
    if (is.null(best)) r$quotient_status <- "no_quotient_path" else { r$quotient_status <- "quotient_path"; r$quotient_hops <- length(best) - 1L; r$quotient_path <- paste(kc[best], collapse = ">") }
    r$cells_connected <- any(comp[Ac] %in% comp[Cc])
    if (r$cells_connected) { cp <- .gsx_cell_path(g, Ac, Cc); r$cell_hops <- cp$hops; r$cell_path_clusters <- paste(rle(cluster[cp$cells])$values, collapse = ">")
      paths[[paste(st[x1], st[x2], sep = ">")]] <- cp$cells }
    r$path_agreement <- if (r$cells_connected) "cell_path" else if (r$quotient_status == "quotient_path") "cluster_path_only" else "no_path"
    rows[[length(rows) + 1L]] <- r }
  d <- do.call(rbind, rows); rownames(d) <- NULL
  d$note <- rep("quotient path = cluster summary graph only, not a continuous cell trajectory; cell connectivity uses the actual cell graph; neither proves differentiation", nrow(d))
  attr(d, "cell_paths") <- paths; d                                                          # 실제 세포 경로 번호 (무시되는 결과에만 저장)
}
