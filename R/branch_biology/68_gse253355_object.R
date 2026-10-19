# =============================================================================
# 68_gse253355_object.R — GSE253355 Seurat 객체를 SeuratObject 없이 슬롯 그대로 읽기 (raw count · meta.data · reduction · graph · 점검 기록)
# =============================================================================
#
# 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_gse253355_instruction.md 1 항. Ruling: 같은 폴더 claude_gse253355_progress.md
#   - dcta 에는 Seurat / SeuratObject 가 없다(설치하지 않는다). S4 객체의 슬롯은 속성이므로 attr() 로 읽는다. Seurat 함수 · 재정규화 · 통합 assay 를 쓰지 않는다.
#   - raw count 만 쓴다: v5 Assay5 의 counts layer(표본별 분할 layer 는 LogMap 의 유전자 · 세포 이름으로 붙인다) 또는 v3/v4 Assay 의 counts 슬롯. data · scale.data 는 쓰지 않는다.
#   - 분할 layer 의 유전자가 서로 다르면 채우지 않고 멈춘다. 세포 순서는 meta.data 순서다. 정수가 아닌 값 · 음수가 있으면 raw count 가 아니므로 멈춘다.
#   - .part (다운로드 중) 는 읽지 않는다. gzip 한 번 / 두 번(gzcon 으로 한 번 더 푼다) / 무압축을 구분해 기록한다.
# =============================================================================

suppressPackageStartupMessages(library(Matrix))

.gsx_attr <- function(x, name) attr(x, name, exact = TRUE)
.gsx_dim <- function(m) { d <- dim(m); if (is.null(d)) d <- .gsx_attr(m, "Dim"); as.integer(d) }
.gsx_cls <- function(x) { k <- class(x); p <- attr(k, "package"); paste0(k[1], if (!is.null(p)) paste0(" (", p, ")") else "") }

#' RDS 열기: .part 거절, gzip 한 번 · 두 번 · 무압축을 머리 바이트로 판단한다. 반환: list(object, format)
gsx_open_rds <- function(path) {
  if (grepl("\\.part$", path)) stop("다운로드 중인 .part 파일은 읽지 않는다: ", path)
  if (!file.exists(path)) stop("입력 파일이 없다: ", path)
  gz <- as.raw(c(0x1f, 0x8b)); con <- file(path, "rb"); m1 <- readBin(con, "raw", 2L); close(con)
  if (!identical(m1, gz)) return(list(object = readRDS(path), format = "none"))
  g <- gzfile(path, "rb"); m2 <- readBin(g, "raw", 2L); close(g)
  if (!identical(m2, gz)) return(list(object = readRDS(path), format = "gzip"))
  con <- gzcon(gzfile(path, "rb")); on.exit(close(con)); list(object = readRDS(con), format = "double_gzip")
}

#' 희소 행렬 슬롯(i · p · x · Dim · Dimnames)을 dgCMatrix 로 (Graph 등 dgCMatrix 를 확장한 미정의 클래스 포함)
gsx_as_dgc <- function(m) {
  if (identical(class(m)[1], "dgCMatrix")) return(m)
  i <- .gsx_attr(m, "i"); p <- .gsx_attr(m, "p"); x <- .gsx_attr(m, "x"); d <- .gsx_attr(m, "Dim"); dn <- .gsx_attr(m, "Dimnames")
  if (is.null(i) || is.null(p) || is.null(d)) stop("dgCMatrix 형 희소 행렬 슬롯이 아니다: ", .gsx_cls(m))
  methods::new("dgCMatrix", i = as.integer(i), p = as.integer(p), x = as.numeric(if (is.null(x)) rep(1, length(i)) else x), Dim = as.integer(d),
               Dimnames = if (is.null(dn)) list(NULL, NULL) else dn)
}

#' LogMap (논리 행렬 + 이름) → 보통 논리 행렬
.gsx_logmap <- function(x) {
  d <- dim(x); dn <- dimnames(x); if (is.null(d) || is.null(dn) || is.null(dn[[1]]) || is.null(dn[[2]])) stop("LogMap 의 이름(유전자 · 세포 x layer)이 없다.")
  matrix(as.logical(as.vector(x)), d[1], d[2], dimnames = dn)
}

#' raw count (유전자 x 세포 dgCMatrix, 세포 순서 = meta.data). 반환: list(counts, layers, source)
gsx_counts <- function(obj, assay = "RNA") {
  A <- .gsx_attr(obj, "assays")[[assay]]; if (is.null(A)) stop("assay 가 없다: ", assay)
  cells_obj <- rownames(.gsx_attr(obj, "meta.data"))
  if (!is.null(.gsx_attr(A, "layers"))) {
    L <- .gsx_attr(A, "layers"); ln <- grep("^counts(\\.|$)", names(L), value = TRUE)
    if (!length(ln)) stop("counts layer 가 없다 (raw count 만 쓴다): ", paste(names(L), collapse = ", "))
    fm <- .gsx_logmap(.gsx_attr(A, "features")); cm <- .gsx_logmap(.gsx_attr(A, "cells"))
    miss <- setdiff(ln, intersect(colnames(fm), colnames(cm))); if (length(miss)) stop("LogMap 에 없는 counts layer: ", paste(miss, collapse = ", "))
    feats <- lapply(ln, function(l) rownames(fm)[fm[, l]])
    if (length(unique(feats)) != 1L) stop("분할 counts layer 들의 유전자가 서로 다르다 (0 으로 채우지 않고 멈춘다): ", paste(ln, lengths(feats), sep = "=", collapse = ", "))
    parts <- lapply(ln, function(l) { m <- gsx_as_dgc(L[[l]]); d <- .gsx_dim(m)
      if (d[1] != length(feats[[1]]) || d[2] != sum(cm[, l])) stop("layer 크기가 LogMap 과 다르다: ", l, " ", paste(d, collapse = "x"))
      dimnames(m) <- list(feats[[1]], rownames(cm)[cm[, l]]); m })
    X <- if (length(parts) == 1L) parts[[1]] else do.call(cbind, parts); src <- paste0("assays$", assay, "@layers$", ln)
  } else if (!is.null(.gsx_attr(A, "counts"))) {
    X <- gsx_as_dgc(.gsx_attr(A, "counts")); ln <- "counts"; src <- paste0("assays$", assay, "@counts")
  } else stop("raw count 슬롯(layers$counts* 또는 counts)을 찾지 못했다: ", .gsx_cls(A))
  if (is.null(rownames(X)) || is.null(colnames(X))) stop("counts 행렬에 유전자 · 세포 이름이 없다.")
  if (anyDuplicated(rownames(X))) stop("유전자 이름이 중복이다: ", paste(utils::head(unique(rownames(X)[duplicated(rownames(X))]), 5), collapse = ", "))
  if (anyDuplicated(colnames(X))) stop("세포 ID 가 중복이다.")
  if (!is.null(cells_obj)) { if (!setequal(colnames(X), cells_obj)) stop("counts 세포와 meta.data 세포가 다르다.")
    if (!identical(colnames(X), cells_obj)) X <- X[, cells_obj, drop = FALSE] }
  if (length(X@x) && (any(!is.finite(X@x)) || any(X@x < 0) || any(X@x != round(X@x)))) stop("counts 가 raw count 가 아니다 (음수 · 비정수 · 비유한 값).")
  list(counts = X, layers = ln, source = src)
}

gsx_meta <- function(obj) { md <- .gsx_attr(obj, "meta.data"); if (!is.data.frame(md) || is.null(rownames(md))) stop("meta.data 가 없다."); md }

#' reduction 좌표 (세포 x 차원, 행 이름 = 세포 ID)
gsx_umap <- function(obj, name) {
  R <- .gsx_attr(obj, "reductions")[[name]]; if (is.null(R)) stop("reduction 이 없다: ", name, " (있는 것: ", paste(names(.gsx_attr(obj, "reductions")), collapse = ", "), ")")
  E <- .gsx_attr(R, "cell.embeddings"); if (is.null(E)) stop("cell.embeddings 가 없다: ", name)
  E <- as.matrix(E); if (is.null(rownames(E))) stop("reduction 좌표에 세포 이름이 없다: ", name); E
}

#' 원 그래프 (세포 x 세포 dgCMatrix)
gsx_graph <- function(obj, name) {
  G <- .gsx_attr(obj, "graphs")[[name]]; if (is.null(G)) stop("그래프가 없다: ", name, " (있는 것: ", paste(names(.gsx_attr(obj, "graphs")), collapse = ", "), ")")
  gsx_as_dgc(G)
}

#' 그래프 규칙 기록: 정사각 · 대칭 · 음수 없음 · 대각 · 양수 상삼각 edge 수
gsx_graph_checks <- function(G) {
  d <- .gsx_dim(G); dg <- Matrix::diag(G); up <- Matrix::triu(G, k = 1)
  list(n_rows = d[1], n_cols = d[2], square = d[1] == d[2], symmetric = isTRUE(Matrix::isSymmetric(G)), nonnegative = !length(G@x) || all(G@x >= 0),
       n_diagonal = as.integer(sum(dg != 0)), diagonal_range = if (any(dg != 0)) paste(signif(range(dg[dg != 0]), 4), collapse = "–") else "none",
       n_upper_edges = as.integer(sum(up@x > 0)), weight_range = if (length(G@x)) paste(signif(range(G@x), 4), collapse = "–") else "none", n_stored = length(G@x))
}

#' 세포 ID 정렬: counts · meta.data · UMAP · 그래프가 같은 순서여야 한다 (재정렬하지 않는다)
gsx_check_ids <- function(cells, meta = NULL, umap = NULL, graph = NULL) {
  chk <- list(meta = if (!is.null(meta)) rownames(meta), umap = if (!is.null(umap)) rownames(umap), graph_rows = if (!is.null(graph)) rownames(graph),
              graph_cols = if (!is.null(graph)) colnames(graph))
  for (k in names(chk)) if (!is.null(chk[[k]]) && !identical(as.character(chk[[k]]), as.character(cells))) stop("세포 ID 순서가 다르다 (", k, "): 재정렬하지 않고 멈춘다.")
  TRUE
}

#' 객체 점검 기록: class · 슬롯 · assay / layer 크기 · meta.data 열 · reduction · graph (값을 바꾸지 않는다)
gsx_inspect <- function(obj) {
  sl <- function(x) setdiff(names(attributes(x)), "class")
  as_ <- .gsx_attr(obj, "assays"); md <- .gsx_attr(obj, "meta.data"); rd <- .gsx_attr(obj, "reductions"); gr <- .gsx_attr(obj, "graphs")
  assays <- lapply(names(as_), function(a) { A <- as_[[a]]; L <- .gsx_attr(A, "layers")
    list(name = a, class = .gsx_cls(A), slots = sl(A),
         layers = if (!is.null(L)) lapply(names(L), function(l) list(name = l, class = .gsx_cls(L[[l]]), dim = .gsx_dim(L[[l]]))) else
           lapply(intersect(c("counts", "data", "scale.data"), sl(A)), function(l) list(name = l, class = .gsx_cls(.gsx_attr(A, l)), dim = .gsx_dim(.gsx_attr(A, l))))) })
  list(object_class = .gsx_cls(obj), object_slots = sl(obj), version = tryCatch(as.character(.gsx_attr(obj, "version")), error = function(e) ""),
       project_name = as.character(.gsx_attr(obj, "project.name")), active_assay = as.character(.gsx_attr(obj, "active.assay")),
       n_cells = nrow(md), assays = assays,
       reductions = lapply(names(rd), function(r) { E <- .gsx_attr(rd[[r]], "cell.embeddings"); list(name = r, class = .gsx_cls(rd[[r]]), dim = .gsx_dim(E),
                                                                                                       key = as.character(.gsx_attr(rd[[r]], "key")), assay = as.character(.gsx_attr(rd[[r]], "assay.used"))) }),
       graphs = lapply(names(gr), function(g) list(name = g, class = .gsx_cls(gr[[g]]), dim = .gsx_dim(gr[[g]]), n_stored = length(.gsx_attr(gr[[g]], "x")),
                                                  assay = as.character(.gsx_attr(gr[[g]], "assay.used")))))
}

#' meta.data 열 표: 이름 · 형 · 결측 · 고유값 수 · 예시 (계산에 쓰지 않는 서술 기록)
gsx_meta_columns <- function(md) {
  do.call(rbind, lapply(names(md), function(k) { v <- md[[k]]; u <- unique(as.character(v)); u <- u[!is.na(u)]
    data.frame(column = k, class = class(v)[1], n_missing = sum(is.na(v)), n_unique = length(u),
               examples = substr(paste(utils::head(sort(u), 8), collapse = " | "), 1, 300), stringsAsFactors = FALSE) }))
}
