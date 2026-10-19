# =============================================================================
# 05_cell_trajectory_io.R — 세포 단위 경로 비교용 입력 어댑터 (R 전용)
# =============================================================================
#
# 목적
#   저장된 세포 이웃 그래프(scipy save_npz 형식 .npz)와 BoneMarrowMap h5ad 를
#   R 에서 직접 읽어, 세포 ID 순서를 엄격히 맞춘 적합 입력을 만든다.
#   Python/reticulate 를 쓰지 않는다. 형식이 예상과 다르면 추측하지 않고 멈춘다.
#
# 적합 입력(ct_load_inputs)이 돌려주는 것
#   cell_id     세포 ID (graph_ids.txt 순서 = h5ad obs 순서여야 한다)
#   latent      Harmony 50차원에서 [1, S.Score, G2M.Score] 최소제곱 잔차를 뺀 좌표
#               (기존 dcta/projection_cluster_context.py::residualize_latent 와 같은 정의,
#                z-scaling 없음). 유전자 수준 회귀가 아니다.
#   graph       세포 x 세포 대칭 sparse 그래프 (graph_base.npz 그대로)
#   cluster_id  저장된 탐색 파티션의 unit (문자열). 원래 membership 은 바꾸지 않는다.
#   donor, cycle(S, G2M, phase), hashes
#   저자 라벨과 UMAP 은 돌려주지 않는다. 그림 단계만 ct_load_display() 로 따로 읽는다.
#
# 예시
#   inp <- ct_load_inputs("data/external/hay_hca_source.h5ad",
#                         ".../graph_base.npz", ".../graph_ids.txt",
#                         ".../res4.0_seed0/membership.tsv.gz")
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(hdf5r)
})


# ------------------------------------------------------------------ NPY / NPZ

#' .npy 파일 하나를 읽는다 (numpy 형식 1.0/2.0, C order 만).
#'
#' 지원 dtype: <i4, <i8, <f8, <f4, |S<n> (바이트 문자열). little-endian 만 허용한다.
#' 모르는 dtype, big-endian, fortran order, 잘린 파일은 오류로 멈춘다.
ct_read_npy <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", 6L)
  if (length(magic) < 6L || !identical(magic, as.raw(c(0x93, charToRaw("NUMPY"))))) {
    stop("NPY magic 이 아니다: ", path)
  }
  ver <- as.integer(readBin(con, "raw", 2L))
  hlen <- if (ver[1] == 1L) readBin(con, "integer", 1L, size = 2L, signed = FALSE, endian = "little") else
    readBin(con, "integer", 1L, size = 4L, endian = "little")
  header <- rawToChar(readBin(con, "raw", hlen))
  descr <- sub(".*'descr':\\s*'([^']+)'.*", "\\1", header)
  fortran <- grepl("'fortran_order':\\s*True", header)
  shape_txt <- sub(".*'shape':\\s*\\(([^)]*)\\).*", "\\1", header)
  dims <- as.numeric(strsplit(gsub("\\s", "", shape_txt), ",")[[1]])
  dims <- dims[!is.na(dims)]
  n <- if (length(dims)) prod(dims) else 1
  if (fortran) stop("fortran_order 배열은 지원하지 않는다: ", path)
  if (startsWith(descr, ">")) stop("big-endian dtype 은 지원하지 않는다 (", descr, "): ", path)
  read_n <- function(what, size, signed = TRUE) {
    x <- readBin(con, what, n, size = size, signed = signed, endian = "little")
    if (length(x) != n) stop("NPY 자료가 잘렸다 (", length(x), " / ", n, "): ", path)
    x
  }
  if (descr %in% c("<i4", "=i4")) {
    out <- as.numeric(read_n("integer", 4L))
  } else if (descr %in% c("<i8", "=i8")) {
    lo_hi <- readBin(con, "integer", 2 * n, size = 4L, endian = "little")
    if (length(lo_hi) != 2 * n) stop("NPY 자료가 잘렸다: ", path)
    lo <- lo_hi[c(TRUE, FALSE)]; hi <- lo_hi[c(FALSE, TRUE)]
    lo_u <- ifelse(lo < 0, lo + 2^32, lo)
    out <- hi * 2^32 + lo_u
    if (any(abs(out) > 2^53)) stop("int64 값이 double 정밀도를 넘는다: ", path)
  } else if (descr %in% c("<f8", "=f8")) {
    out <- read_n("double", 8L)
  } else if (descr %in% c("<f4", "=f4")) {
    out <- read_n("double", 4L)
  } else if (grepl("^\\|S[0-9]+$", descr)) {
    w <- as.integer(sub("^\\|S", "", descr))
    b <- readBin(con, "raw", w * n)
    if (length(b) != w * n) stop("NPY 자료가 잘렸다: ", path)
    out <- vapply(seq_len(n), function(i) {
      r <- b[((i - 1L) * w + 1L):(i * w)]; rawToChar(r[r != as.raw(0)])
    }, character(1))
  } else {
    stop("지원하지 않는 NPY dtype: ", descr, " (", path, ")")
  }
  out
}


#' scipy.sparse.save_npz 형식(.npz = zip 안의 data/indices/indptr/format/shape)을
#' Matrix 의 CsparseMatrix 로 읽는다. CSR 과 CSC 만 허용한다.
#' 포인터/인덱스/길이를 모두 검사하고, 맞지 않으면 멈춘다.
ct_read_npz_sparse <- function(path) {
  members <- utils::unzip(path, list = TRUE)$Name
  need <- c("data.npy", "indices.npy", "indptr.npy", "format.npy", "shape.npy")
  if (!all(need %in% members)) {
    stop("npz 에 필요한 배열이 없다: ", paste(setdiff(need, members), collapse = ", "))
  }
  d <- tempfile("npz_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  utils::unzip(path, files = need, exdir = d)
  rd <- function(nm) ct_read_npy(file.path(d, nm))
  fmt <- rd("format.npy"); shape <- rd("shape.npy")
  x <- rd("data.npy"); idx <- rd("indices.npy"); ptr <- rd("indptr.npy")
  if (!fmt %in% c("csr", "csc")) stop("지원하지 않는 sparse format: ", fmt)
  if (length(shape) != 2L) stop("shape 는 길이 2 여야 한다.")
  n_major <- if (fmt == "csr") shape[1] else shape[2]
  n_minor <- if (fmt == "csr") shape[2] else shape[1]
  if (length(ptr) != n_major + 1) stop("indptr 길이가 맞지 않는다: ", length(ptr), " != ", n_major + 1)
  if (ptr[1] != 0 || any(diff(ptr) < 0)) stop("indptr 가 0 에서 시작하는 단조 증가가 아니다.")
  if (length(idx) != length(x) || ptr[length(ptr)] != length(x)) {
    stop("data/indices 길이가 indptr 끝값과 맞지 않는다 (", length(x), ", ", length(idx), ", ",
         ptr[length(ptr)], ").")
  }
  if (length(idx) && (any(idx < 0) || any(idx >= n_minor) || any(idx != round(idx)))) {
    stop("indices 가 범위를 벗어났거나 정수가 아니다.")
  }
  # CSC 이면 p = 열 포인터, i = 행 인덱스. CSR 이면 전치의 CSC 로 만든 뒤 다시 전치한다.
  if (fmt == "csc") {
    G <- Matrix::sparseMatrix(i = idx + 1, p = ptr, x = x, dims = shape)
  } else {
    G <- Matrix::t(Matrix::sparseMatrix(i = idx + 1, p = ptr, x = x, dims = rev(shape)))
  }
  methods::as(G, "CsparseMatrix")
}


#' 세포 그래프 검사: n x n 정사각, 유한·음이 아닌 값, 대칭.
#' 대각 성분은 바꾸지 않고 개수만 attr 로 남긴다.
ct_validate_graph <- function(G, n_cells) {
  if (nrow(G) != ncol(G)) stop("그래프가 정사각 행렬이 아니다: ", nrow(G), " x ", ncol(G))
  if (nrow(G) != n_cells) stop("그래프 크기 ", nrow(G), " 가 세포 수 ", n_cells, " 와 다르다.")
  if (any(!is.finite(G@x)) || any(G@x < 0)) stop("그래프 값에 비유한값 또는 음수가 있다.")
  asym <- if (length(G@x)) max(abs(G - Matrix::t(G))) else 0
  if (asym > 1e-12) stop("그래프가 대칭이 아니다 (최대 차이 ", format(asym), ").")
  attr(G, "n_diagonal_nonzero") <- sum(Matrix::diag(G) != 0)
  G
}


# ------------------------------------------------------------------ ID / H5

#' 두 세포 ID 벡터가 중복 없이 같은 순서인지 확인한다. 암묵적으로 재정렬하지 않는다.
ct_check_ids <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  if (anyDuplicated(a) || anyDuplicated(b)) stop("세포 ID 에 중복이 있다.")
  if (length(a) != length(b) || !identical(a, b)) {
    stop("세포 ID 순서/길이가 다르다 (", length(a), " vs ", length(b), "). 암묵적으로 재정렬하지 않는다.")
  }
  TRUE
}


#' h5ad 의 (cell, dim) 행렬을 cell x dim R 행렬로 읽는다.
#' hdf5r 은 차원을 뒤집어 보여 주므로 R 에서는 dim x cell 이 보이고, 이를 전치한다.
#' 세포 수와 맞지 않으면 방향을 추측하지 않고 멈춘다.
ct_read_h5_matrix <- function(f, path, n_cells) {
  ds <- f[[path]]
  dm <- ds$dims
  if (length(dm) != 2L || dm[2] != n_cells) {
    stop(path, " 의 방향/세포 수가 맞지 않는다 (R 에서 보이는 차원 ", paste(dm, collapse = " x "),
         ", 기대 세포 수 ", n_cells, ").")
  }
  t(ds$read())
}


#' h5ad obs 열을 읽는다. categorical(categories + codes)이면 문자열로 바꾼다.
ct_h5_obs <- function(f, name) {
  o <- f[[paste0("obs/", name)]]
  if (inherits(o, "H5Group")) {
    cats <- as.character(o[["categories"]]$read()); codes <- as.integer(o[["codes"]]$read())
    if (any(codes >= length(cats) | codes < -1L, na.rm = TRUE)) {
      stop("obs/", name, " 의 categorical code 가 범위를 벗어났다 (범주 ", length(cats), "개).")
    }
    # 원래 길이의 NA 벡터를 먼저 만들고, 유효 코드 위치만 채운다 (code -1 = 결측; 위치가 밀리지 않는다)
    out <- rep(NA_character_, length(codes))
    ok <- !is.na(codes) & codes >= 0L
    out[ok] <- cats[codes[ok] + 1L]
    return(out)
  }
  o$read()
}


#' 세포주기 점수로 잠재 좌표를 잔차화한다: H - X (X'X)^{-1} X'H, X = [1, S, G2M].
ct_residualize <- function(H, S, G2M) {
  if (!is.matrix(H)) stop("H 는 행렬이어야 한다.")
  if (length(S) != nrow(H) || length(G2M) != nrow(H)) stop("세포주기 점수 길이가 세포 수와 다르다.")
  if (any(!is.finite(H)) || any(!is.finite(S)) || any(!is.finite(G2M))) {
    stop("잠재 좌표 또는 세포주기 점수에 비유한값이 있다.")
  }
  X <- cbind(1, as.numeric(S), as.numeric(G2M))
  # rank-aware 최소제곱 잔차: 공선/상수 열이 있어도 lm.fit / np.linalg.lstsq 잔차와 같다
  q <- qr(X)
  R <- qr.resid(q, H)
  if (any(!is.finite(R))) stop("잔차 좌표가 유한하지 않다.")
  unname(R)
}


# ------------------------------------------------------------------ 전체 입력

#' 적합 입력을 만든다. 저자 라벨과 UMAP 은 읽지 않는다.
ct_load_inputs <- function(source_path, graph_path, ids_path, membership_path) {
  hashes <- c(source = digest::digest(source_path, algo = "sha256", file = TRUE),
              graph_base.npz = digest::digest(graph_path, algo = "sha256", file = TRUE),
              graph_ids.txt = digest::digest(ids_path, algo = "sha256", file = TRUE),
              membership = digest::digest(membership_path, algo = "sha256", file = TRUE))
  gids <- readLines(ids_path); gids <- gids[nzchar(gids)]
  f <- hdf5r::H5File$new(source_path, "r"); on.exit(f$close_all())
  idx_key <- hdf5r::h5attr(f[["obs"]], "_index")
  obs_ids <- as.character(f[[paste0("obs/", idx_key)]]$read())
  ct_check_ids(gids, obs_ids)
  n <- length(gids)
  H <- ct_read_h5_matrix(f, "obsm/X_harmony", n)
  S <- as.numeric(ct_h5_obs(f, "S.Score")); G2M <- as.numeric(ct_h5_obs(f, "G2M.Score"))
  latent <- ct_residualize(H, S, G2M)
  donor <- ct_h5_obs(f, "donor_id"); phase <- ct_h5_obs(f, "CyclePhase")
  graph <- ct_validate_graph(ct_read_npz_sparse(graph_path), n)
  memb <- utils::read.delim(gzfile(membership_path), stringsAsFactors = FALSE, colClasses = c(cell = "character"))
  ct_check_ids(gids, memb$cell)
  list(cell_id = gids, latent = latent, graph = graph,
       cluster_id = as.character(memb$unit), donor = as.character(donor),
       cycle = data.frame(S = S, G2M = G2M, phase = as.character(phase), stringsAsFactors = FALSE),
       hashes = hashes,
       latent_definition = "obsm/X_harmony minus least-squares fit on [1, S.Score, G2M.Score]; no scaling")
}


#' 그림 전용 표시 자료: 원래 저자 UMAP 과 전체 저자 주석. 적합 함수에 넘기지 않는다.
#' author_label = obs/AuthorCellType (전체 저자 주석, Codex R1). ontology 형식의 obs/cell_type 은
#' 저자 라벨 이름으로 쓰지 않는다. 있으면 AuthorCellType_Broad 를 author_label_broad 로 따로 둔다.
ct_load_display <- function(source_path, cell_id) {
  f <- hdf5r::H5File$new(source_path, "r"); on.exit(f$close_all())
  idx_key <- hdf5r::h5attr(f[["obs"]], "_index")
  ct_check_ids(cell_id, as.character(f[[paste0("obs/", idx_key)]]$read()))
  U <- ct_read_h5_matrix(f, "obsm/X_umap", length(cell_id))
  if (!f$exists("obs/AuthorCellType")) stop("obs/AuthorCellType 이 없다 (cell_type 으로 대신하지 않는다).")
  out <- data.frame(cell_id = cell_id, umap_1 = U[, 1], umap_2 = U[, 2],
                    author_label = ct_h5_obs(f, "AuthorCellType"), stringsAsFactors = FALSE)
  out$author_label_broad <- if (f$exists("obs/AuthorCellType_Broad")) ct_h5_obs(f, "AuthorCellType_Broad") else NA_character_
  attr(out, "label_source") <- "obs/AuthorCellType"
  out
}
