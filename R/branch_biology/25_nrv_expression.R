# =============================================================================
# 25_nrv_expression.R — 원 RNA 를 chunk 로 읽어 집합 수준 합을 모으고, 유전자 대조를 고정한다
# =============================================================================
#
# 식 (R17 과 같다)
#   x_cg = log1p(10^4 · count_cg / L_c),  L_c = raw/X 의 전체 유전자 합,  같은 symbol 의 여러 열은 log 전에 합친다
#   검출 d_cg = 1[count_cg > 0]
#   집합 S 의 합 = Σ_{c∈S} x_cg  → 집합 지시행렬 H 와 chunk 행렬 곱으로 누적한다(전체 dense RNA 를 만들지 않는다)
#
# 유전자 대조 두 가지
#   - anchor 대조 유전자: legacy module_v1 select_controls 규칙을 R 로 옮긴 것(평균 순위 25 bin, 차지한 bin 마다 50 개).
#     Python 난수열까지 재현하지는 않는다.
#   - gene null 집합: 대상 유전자마다 (전역 평균 log 10 분위 x 검출률 10 분위) 같은 bin 에서 하나씩 뽑는다.
#     이것은 bin 을 맞춘 무작위 집합과의 비교량이다. 공발현·조절관계까지 교환 가능하게 만들지 않는다.
#
# chunk 반복자: nrv_chunks_matrix(X) 와 nrv_chunks_h5(path) 는 같은 인터페이스 function(callback(X_chunk, rows)) 를 돌려준다.
# =============================================================================

suppressPackageStartupMessages(library(Matrix))

#' symbol → raw 열 번호 목록 (없는 유전자는 attr "missing")
nrv_gene_map <- function(symbols, genes) {
  symbols <- as.character(symbols); m <- setNames(lapply(genes, function(g) which(symbols == g)), genes)
  attr(m, "missing") <- genes[lengths(m) == 0]; m
}

.nrv_P <- function(gene_map, n_cols) {
  ok <- lengths(gene_map) > 0
  Matrix::sparseMatrix(i = unlist(gene_map[ok]), j = rep(seq_len(sum(ok)), lengths(gene_map[ok])), x = 1,
                       dims = c(n_cols, sum(ok)), dimnames = list(NULL, names(gene_map)[ok]))
}

#' chunk 하나의 log 발현과 검출 (세포 x 유전자, 희소)
nrv_log_chunk <- function(X, gene_map) {
  X <- methods::as(X, "CsparseMatrix")
  if (length(X@x) && (any(!is.finite(X@x)) || any(X@x < 0) || any(X@x != round(X@x)))) stop("raw count 가 아니다(비유한, 음수, 정수 아님).")
  lib <- Matrix::rowSums(X); if (any(lib <= 0)) stop("library 합이 0 인 세포가 있다.")
  cnt <- methods::as(X %*% .nrv_P(gene_map, ncol(X)), "CsparseMatrix")       # 중복 symbol 은 log 전에 합친다
  det <- cnt; det@x <- as.numeric(det@x > 0)
  L <- methods::as(Matrix::Diagonal(x = 1e4 / lib) %*% cnt, "CsparseMatrix"); L@x <- log1p(L@x)
  list(library = unname(lib), log = L, det = det)
}

#' 메모리 행렬의 chunk 반복자 (테스트용, 같은 인터페이스)
nrv_chunks_matrix <- function(X, chunk = 5000L) function(callback) {
  n <- nrow(X); for (s in seq(1L, n, by = chunk)) { r <- s:min(n, s + chunk - 1L); callback(X[r, , drop = FALSE], r) }
}

#' h5ad raw/X(CSR) 의 chunk 반복자
nrv_chunks_h5 <- function(h5_path, chunk = 5000L) function(callback) {
  f <- hdf5r::H5File$new(h5_path, "r"); on.exit(f$close_all()); rx <- f[["raw/X"]]
  if (!identical(as.character(hdf5r::h5attr(rx, "encoding-type")), "csr_matrix")) stop("raw/X 가 csr_matrix 가 아니다.")
  shp <- as.integer(hdf5r::h5attr(rx, "shape")); ptr <- as.numeric(rx[["indptr"]]$read()); ds <- rx[["data"]]; di <- rx[["indices"]]
  for (s0 in seq(0, shp[1] - 1, by = chunk)) {
    s1 <- min(s0 + chunk, shp[1]); a <- ptr[s0 + 1]; b <- ptr[s1 + 1]
    x <- if (b > a) as.numeric(ds[(a + 1):b]) else numeric(0); j <- if (b > a) as.integer(di[(a + 1):b]) + 1L else integer(0)
    X <- Matrix::sparseMatrix(i = rep.int(seq_len(s1 - s0), diff(ptr[(s0 + 1):(s1 + 1)])), j = j, x = x, dims = c(s1 - s0, shp[2]))
    callback(X, (s0 + 1):s1)
  }
}

#' 집합 수준 합: sum(집합 x 유전자), det(검출 세포 수), n(집합 크기)
nrv_stream_set_sums <- function(each_chunk, gene_map, H) {
  g <- names(gene_map)[lengths(gene_map) > 0]; S <- matrix(0, nrow(H), length(g), dimnames = list(rownames(H), g)); D <- S
  each_chunk(function(X, rows) {
    lc <- nrv_log_chunk(X, gene_map); Hc <- H[, rows, drop = FALSE]
    S <<- S + as.matrix(Hc %*% lc$log); D <<- D + as.matrix(Hc %*% lc$det)
  })
  list(sum = S, det = D, n = Matrix::rowSums(H))
}

#' 정해진 세포만 세포 수준 log 발현으로 모은다(경로 내부 군집 세포용, 희소)
nrv_read_cells <- function(each_chunk, gene_map, cells) {
  keep <- sort(unique(cells)); parts <- list()
  each_chunk(function(X, rows) {
    hit <- which(rows %in% keep); if (!length(hit)) return(invisible(NULL))
    lc <- nrv_log_chunk(X[hit, , drop = FALSE], gene_map); rownames(lc$log) <- rows[hit]
    parts[[length(parts) + 1L]] <<- lc$log
  })
  M <- do.call(rbind, parts); M[as.character(keep), , drop = FALSE]
}

#' 전체 유전자 전역 통계 (raw 열 단위): 평균 log 발현, 검출률
nrv_gene_global_stats <- function(each_chunk, symbols) {
  ng <- length(symbols); s <- numeric(ng); d <- numeric(ng); n <- 0
  each_chunk(function(X, rows) {
    X <- methods::as(X, "CsparseMatrix"); lib <- Matrix::rowSums(X)
    L <- methods::as(Matrix::Diagonal(x = 1e4 / lib) %*% X, "CsparseMatrix"); L@x <- log1p(L@x)
    s <<- s + Matrix::colSums(L); d <<- d + Matrix::colSums(X > 0); n <<- n + nrow(X)
  })
  data.frame(col = seq_len(ng), symbol = symbols, mean_log = s / n, detect_frac = d / n, stringsAsFactors = FALSE)
}

.nrv_with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed); force(expr)
}

#' anchor 대조 유전자 (legacy select_controls 규칙): Q 가 차지한 bin 마다 Q·제외 목록 밖에서 최대 ctrl_size 개
nrv_select_controls <- function(gs, Q_cols, exclude_cols, n_bins = 25L, ctrl_size = 50L, seed = 20261013L) {
  pool <- gs[is.finite(gs$mean_log) & gs$mean_log > 0, ]; n_items <- round(nrow(pool) / (n_bins - 1L))
  pool$bin <- rank(pool$mean_log, ties.method = "min") %/% n_items
  qb <- sort(unique(pool$bin[match(Q_cols, pool$col)])); if (anyNA(qb)) stop("신호 유전자가 발현 pool 에 없다.")
  .nrv_with_seed(seed, do.call(rbind, lapply(qb, function(b) {
    cand <- pool$col[pool$bin == b & !(pool$col %in% c(Q_cols, exclude_cols))]
    if (!length(cand)) stop("bin ", b, " 에 대조 후보가 없다.")
    data.frame(bin = b, col = if (length(cand) > ctrl_size) sort(cand[sample.int(length(cand), ctrl_size)]) else cand)
  })))
}

#' gene null 집합 (R x 대상 수 정수 행렬). attr "bin" 에 열별 bin 을 남긴다.
nrv_null_sets <- function(gs, target_cols, exclude_cols, R = 1000L, n_bins = 10L, seed = 20261013L) {
  br <- function(v) stats::quantile(v, seq(0, 1, length.out = n_bins + 1L), names = FALSE)
  bin <- paste(findInterval(gs$mean_log, br(gs$mean_log), rightmost.closed = TRUE, all.inside = TRUE),
               findInterval(gs$detect_frac, br(gs$detect_frac), rightmost.closed = TRUE, all.inside = TRUE))
  ok <- !(gs$col %in% exclude_cols); tb <- bin[match(target_cols, gs$col)]
  N <- .nrv_with_seed(seed, t(vapply(seq_len(R), function(r) {
    used <- integer(0)
    for (k in seq_along(target_cols)) {
      cand <- setdiff(gs$col[ok & bin == tb[k]], used); if (!length(cand)) stop("null bin 에 후보가 없다: ", tb[k])
      used <- c(used, cand[sample.int(length(cand), 1L)])
    }
    used
  }, integer(length(target_cols)))))
  if (length(target_cols) == 1L) N <- matrix(N, ncol = 1L)
  attr(N, "bin") <- setNames(bin, gs$col); N
}

#' 여러 (집합 행렬, 유전자 묶음) 을 raw 한 번 읽기로 모은다.
#'   maps    = 이름 붙은 gene_map 목록 (예: ref = 참고 유전자 symbol, pool = raw 열 "col<번호>")
#'   targets = 이름 붙은 list(H = 집합 x 세포 행렬(가중치 가능), map = maps 의 이름, det = 검출 세포 수도 셀지)
#'   결과는 target 마다 nrv_stream_set_sums 와 같은 list(sum, det, n). det = FALSE 이면 det 는 NULL, n = rowSums(H).
#'   유전자가 많으면 gene_block 개씩 나눠 곱해 임시 dense 행렬을 작게 둔다(전체 dense RNA 는 만들지 않는다).
nrv_stream_multi <- function(each_chunk, maps, targets, gene_block = 500L) {
  used <- unique(vapply(targets, function(tg) tg$map, character(1)))
  if (!all(used %in% names(maps))) stop("target 이 모르는 유전자 묶음을 쓴다: ", paste(setdiff(used, names(maps)), collapse = ", "))
  zero <- function(tg) { m <- maps[[tg$map]]; g <- names(m)[lengths(m) > 0]; matrix(0, nrow(tg$H), length(g), dimnames = list(rownames(tg$H), g)) }
  acc <- lapply(targets, function(tg) { e <- new.env(); e$sum <- zero(tg); e$det <- if (isTRUE(tg$det)) zero(tg) else NULL; e })   # 제자리 누적
  each_chunk(function(X, rows) {
    lc <- lapply(maps[used], function(m) nrv_log_chunk(X, m))
    for (k in names(targets)) {
      tg <- targets[[k]]; e <- acc[[k]]; Hc <- tg$H[, rows, drop = FALSE]; l <- lc[[tg$map]]; ng <- ncol(l$log)
      if (!ng) next
      for (s in seq(1L, ng, by = gene_block)) {
        j <- s:min(ng, s + gene_block - 1L)
        e$sum[, j] <- e$sum[, j, drop = FALSE] + as.matrix(Hc %*% l$log[, j, drop = FALSE])
        if (isTRUE(tg$det)) e$det[, j] <- e$det[, j, drop = FALSE] + as.matrix(Hc %*% l$det[, j, drop = FALSE])
      }
    }
  })
  setNames(lapply(names(targets), function(k) list(sum = acc[[k]]$sum, det = acc[[k]]$det, n = Matrix::rowSums(targets[[k]]$H))), names(targets))
}
