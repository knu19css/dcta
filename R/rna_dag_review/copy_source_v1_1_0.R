# =============================================================================
# copy_source_v1_1_0.R — 문헌 기반 조혈 DAG v1.1.0 의 작은 표를 검토용으로 명시 복사하고 SHA-256 목록을 남긴다
# =============================================================================
# 사용: Rscript R/rna_dag_review/copy_source_v1_1_0.R --source-dir <v1.1.0 표 폴더> [--out-dir data/priors/rna_dag_review/source_v1_1_0]
#   byte_identical : 바이트 그대로 복사 (nodes, edges, edge_evidence, taxonomy, gate_membership)
#   column_subset  : evidence.tsv 에서 검토자 기록 열(reviewed_by)을 빼고, review_level 값 끝의 검토자 접미사를 떼어 문헌 검토 수준만 남긴다
#   withheld       : 공개 복사본에 넣지 않는다 (자유 서술에 내부 작업 메모가 있고 빌더가 쓰지 않음). 원본 SHA 만 기록
# 실측 세포 자료 · 원 RNA 는 복사하지 않는다. 같은 이름의 다른 파일이 이미 있으면 덮어쓰지 않고 멈춘다.
# =============================================================================

RDR_COPY_BYTE_IDENTICAL <- c("nodes.tsv", "edges.tsv", "edge_evidence.tsv", "taxonomy.tsv", "gate_membership.tsv")
RDR_COPY_WITHHELD <- c("adult_evidence_audit.tsv")
RDR_REVIEW_LEVELS <- c("abstract_only", "abstract_plus_search", "abstract_plus_captions", "full_text_partial", "full_text_with_captions")

.rdr_sha <- function(path) digest::digest(file = path, algo = "sha256")

.rdr_read_raw_tsv <- function(path) utils::read.delim(path, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
                                                      colClasses = "character", na.strings = character(0), fileEncoding = "UTF-8")

.rdr_write_raw_tsv <- function(x, path) {
  cells <- unlist(x, use.names = FALSE)
  if (anyNA(cells) || any(grepl("[\t\n\r]", c(cells, names(x))))) stop("cannot write ", basename(path), ": NA, tab or newline inside a cell")
  lines <- c(paste(names(x), collapse = "\t"), do.call(paste, c(unname(as.list(x)), sep = "\t")))
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
}

#' review_level 값에서 문헌 검토 수준만 남긴다 (예: "<level>_<검토자 접미사>" -> "<level>"). 알 수 없는 값은 멈춘다.
rdr_normalize_review_level <- function(x) {
  out <- x
  for (i in which(!x %in% RDR_REVIEW_LEVELS)) {
    hit <- RDR_REVIEW_LEVELS[startsWith(x[i], paste0(RDR_REVIEW_LEVELS, "_"))]
    if (!length(hit)) stop("unknown review_level value: ", x[i])
    out[i] <- hit[which.max(nchar(hit))]
  }
  out
}

#' evidence.tsv -> 공개 부분집합 (reviewed_by 열 제거, review_level 정규화). 나머지 열 · 행 · 값은 그대로.
rdr_evidence_public_subset <- function(ev) {
  if (!all(c("evidence_id", "review_level", "reviewed_by") %in% names(ev))) stop("evidence.tsv lacks evidence_id / review_level / reviewed_by")
  ev$review_level <- rdr_normalize_review_level(ev$review_level)
  ev[, setdiff(names(ev), "reviewed_by"), drop = FALSE]
}

.rdr_put <- function(write_fun, path) {
  tmp <- tempfile(fileext = ".tsv"); write_fun(tmp)
  if (file.exists(path)) {
    if (!identical(.rdr_sha(path), .rdr_sha(tmp))) stop("refusing to overwrite a different existing file: ", basename(path))
    return(invisible(path))
  }
  if (!file.copy(tmp, path)) stop("copy failed: ", basename(path))
  invisible(path)
}

#' 원본 폴더에서 검토용 복사본과 SHA256SUMS.tsv 를 만든다.
rdr_copy_source <- function(source_dir, out_dir) {
  need <- c(RDR_COPY_BYTE_IDENTICAL, "evidence.tsv", RDR_COPY_WITHHELD)
  miss <- need[!file.exists(file.path(source_dir, need))]
  if (length(miss)) stop("source folder lacks: ", paste(miss, collapse = ", "))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (f in RDR_COPY_BYTE_IDENTICAL) {
    .rdr_put(function(tmp) file.copy(file.path(source_dir, f), tmp, overwrite = TRUE), file.path(out_dir, f))
    rows[[f]] <- data.frame(file = f, sha256 = .rdr_sha(file.path(out_dir, f)), source_file = f, source_sha256 = .rdr_sha(file.path(source_dir, f)),
                            copy_mode = "byte_identical", note_ko = "원본과 바이트 동일")
  }
  ev <- rdr_evidence_public_subset(.rdr_read_raw_tsv(file.path(source_dir, "evidence.tsv")))
  f <- "evidence_public_subset.tsv"
  .rdr_put(function(tmp) .rdr_write_raw_tsv(ev, tmp), file.path(out_dir, f))
  rows[[f]] <- data.frame(file = f, sha256 = .rdr_sha(file.path(out_dir, f)), source_file = "evidence.tsv", source_sha256 = .rdr_sha(file.path(source_dir, "evidence.tsv")),
                          copy_mode = "column_subset", note_ko = "reviewed_by 열 제거; review_level 은 문헌 검토 수준 값만 남김; 나머지 열 · 행 · 값은 원본 그대로")
  for (f in RDR_COPY_WITHHELD)
    rows[[f]] <- data.frame(file = f, sha256 = "not_copied", source_file = f, source_sha256 = .rdr_sha(file.path(source_dir, f)),
                            copy_mode = "withheld", note_ko = "자유 서술에 내부 작업 메모가 있어 공개 복사본에서 보류; 이 검토의 빌더는 쓰지 않음")
  sums <- do.call(rbind, rows); rownames(sums) <- NULL
  .rdr_put(function(tmp) .rdr_write_raw_tsv(sums, tmp), file.path(out_dir, "SHA256SUMS.tsv"))
  invisible(sums)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  arg <- function(k, default = NULL) { i <- match(k, a); if (is.na(i)) default else a[i + 1L] }
  src <- arg("--source-dir"); if (is.null(src)) stop("--source-dir is required")
  s <- rdr_copy_source(src, arg("--out-dir", file.path("data", "priors", "rna_dag_review", "source_v1_1_0")))
  print(s[, c("file", "copy_mode", "sha256")], row.names = FALSE)
}
