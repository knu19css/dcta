# =============================================================================
# 27_run_nrv.R — 성인 Whole BM HSC → Neutrophil 경로 검증: 단계형 실행기
# =============================================================================
#
# 단계 (attempt 폴더 아래 stage_<이름>/ 에 쓰고, stage_manifest.json 으로 파일 hash 를 고정한다)
#   freeze → structure → stream_global → freeze_controls → stream_sets → select → stream_route → evaluate → verdict → author → view
#   - 앞 단계 manifest 와 파일 hash 가 맞아야 다음 단계가 돈다. 앞의 어느 단계 파일을 읽든 그 단계 manifest 의 hash 를 다시 확인한다.
#   - 이미 쓴 단계 폴더는 덮어쓰지 않는다(exit 6). 실패한 단계는 새 attempt 폴더에서 다시 돌린다.
#   - 저자 라벨(obs/AuthorCellType)은 author 단계에서만 읽는다. UMAP 은 view 단계에서만 읽는다.
#   - 단계마다 소스를 source_snapshot_<단계>/ 로 복사하고, 그 사본을 적재해 실행한다(R18 과 같은 방식).
#   - 인자(분할, 배정, 재표본 수, null 수, block 수, chunk, seed)는 freeze 에서 고정한다. 뒤 단계 인자가 다르면 멈춘다.
#
# 실행
#   run-r dcta R/branch_biology/27_run_nrv.R --stage <단계> --attempt-dir <폴더> --data-root <MAIN 루트> --partition p92|p51 --allocation all
# exit code: 0 정상, 2 인자, 4 입력 불일치·검사 실패·계산 오류, 6 단계 폴더가 비어 있지 않음
# =============================================================================

suppressPackageStartupMessages(library(Matrix))

NRV_STAGES <- c("freeze", "structure", "stream_global", "freeze_controls", "stream_sets", "select", "stream_route", "evaluate", "verdict", "author", "view")
NRV_SOURCES <- c("05_cell_trajectory_io.R", "10_expression_edge_graph.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                 "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "28_view_nrv.R")
# 입력 (모두 읽기 전용, --data-root 아래 상대 경로)
NRV_INPUTS <- c(graph = "results/projection_cluster_context_2026-09-25/attempt_01/BoneMarrowMap/graph_base.npz",
                graph_ids = "results/projection_cluster_context_2026-09-25/attempt_01/BoneMarrowMap/graph_ids.txt",
                p92 = "results/branch_resolution_2026-10-02/attempt_02/partitions/res8.0_seed0/membership.tsv.gz",
                p51 = "results/branch_resolution_2026-10-02/attempt_02/partitions/res4.0_seed0/membership.tsv.gz",
                h5ad = "data/external/hay_hca_source.h5ad",
                eval_mapping = "results/path_context_2026-09-26/attempt_01/BoneMarrowMap/evaluation_mapping_used.tsv",
                readout51 = "results/simple_joint_annotation_2026-10-01/attempt_03/BoneMarrowMap/unit_readouts.tsv.gz")
NRV_ROUTE_KINDS <- c(full_route = "T_terminal", partial_path_early_granule = "T_early")   # 경로 종류 → 끝 anchor
NRV_EARLY_EDGES <- sprintf("E%02d", 1:4)                                                     # 초기 과립 부분 경로의 필수 프로그램 범위
NRV_LAYERS <- c("adult_confirmed", "human_age_unconfirmed")

#' 기준 상수 값 (R23 NRV_THRESHOLDS 에서 읽는다. 값을 여기서 다시 쓰지 않는다)
.nrv_thr <- function(name) {
  v <- NRV_THRESHOLDS$value[NRV_THRESHOLDS$name == name]
  if (length(v) != 1L) stop("기준 상수가 없다: ", name)
  as.numeric(v)
}


# ── 단계 폴더, manifest, hash 사슬 ──────────────────────────────────────────
nrv_stage_dir <- function(attempt, stage) file.path(attempt, paste0("stage_", stage))
.nrv_sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)

nrv_write_stage_manifest <- function(attempt, stage, files, upstream) {
  d <- nrv_stage_dir(attempt, stage); m <- file.path(d, "stage_manifest.json")
  if (file.exists(m)) stop("이미 쓴 단계다(덮어쓰지 않는다): ", stage)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  nd <- normalizePath(d, mustWork = TRUE)                                   # 키 = 단계 폴더 기준 상대 경로 (하위 폴더 그림 포함, Codex 검토)
  rel <- vapply(files, function(p) { np <- normalizePath(p, mustWork = TRUE)
    if (!startsWith(np, paste0(nd, "/"))) stop("단계 폴더 밖 파일은 manifest 에 넣지 않는다: ", p)
    substring(np, nchar(nd) + 2L) }, character(1))
  man <- list(stage = stage, written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), upstream = upstream,
              files = as.list(setNames(vapply(files, .nrv_sha, character(1)), rel)))
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE, null = "null"), m)
  invisible(m)
}

#' 바로 앞 단계 manifest 와 그 파일 hash 를 확인한다. 돌려주는 값(앞 manifest 의 hash)이 이 단계 manifest 의 upstream 이 된다.
nrv_check_upstream <- function(attempt, stage) {
  k <- match(stage, NRV_STAGES); if (is.na(k)) stop("모르는 단계: ", stage)
  if (k == 1L) return(invisible(NULL))
  prev <- NRV_STAGES[k - 1L]; m <- file.path(nrv_stage_dir(attempt, prev), "stage_manifest.json")
  if (!file.exists(m)) stop("앞 단계(", prev, ") manifest 가 없다.")
  man <- jsonlite::fromJSON(m)
  for (f in names(man$files))
    if (!identical(.nrv_sha(file.path(nrv_stage_dir(attempt, prev), f)), man$files[[f]])) stop("앞 단계 파일 hash 가 다르다: ", f)
  # 사슬을 freeze 까지 거꾸로 확인한다: 각 manifest 의 upstream = 바로 앞 manifest 의 hash (Codex 검토)
  mf <- function(s) file.path(nrv_stage_dir(attempt, s), "stage_manifest.json")
  if (k >= 3L) for (j in (k - 1L):2L) {
    mj <- mf(NRV_STAGES[j]); mi <- mf(NRV_STAGES[j - 1L])
    if (!file.exists(mi)) stop("hash 사슬이 끊겼다(manifest 없음): ", NRV_STAGES[j - 1L])
    up <- jsonlite::fromJSON(mj)$upstream
    if (is.null(up) || !identical(up, .nrv_sha(mi))) stop("hash 사슬이 끊겼다: ", NRV_STAGES[j], " 의 upstream 이 ", NRV_STAGES[j - 1L], " manifest 와 다르다.")
  }
  invisible(.nrv_sha(m))
}

#' structure 단계에서 고정한 입력 hash 를 다시 확인한다 (reader 가 입력 파일 목록을 줄 때). 다르면 멈춘다 (Codex 검토)
.nrv_check_inputs <- function(a, reader) {
  if (!is.function(reader$input_files)) return(invisible(NULL))
  fx <- utils::read.delim(nrv_stage_file(a$attempt, "structure", "input_sha256.tsv"), stringsAsFactors = FALSE, colClasses = "character")
  f <- reader$input_files(); miss <- setdiff(fx$input, names(f)[file.exists(f)]); f <- f[names(f) %in% fx$input & file.exists(f)]
  bad <- names(f)[vapply(f, .nrv_sha, "") != fx$sha256[match(names(f), fx$input)]]
  if (length(bad) || length(miss)) stop("입력 hash 가 structure 단계와 다르다: ", paste(c(bad, miss), collapse = ", "))
  invisible(TRUE)
}

#' 앞의 어느 단계 파일이든, 그 단계 manifest 의 hash 와 같을 때만 경로를 돌려준다
nrv_stage_file <- function(attempt, stage, file) {
  d <- nrv_stage_dir(attempt, stage); m <- file.path(d, "stage_manifest.json")
  if (!file.exists(m)) stop("앞 단계(", stage, ") manifest 가 없다.")
  h <- jsonlite::fromJSON(m)$files[[file]]; p <- file.path(d, file)
  if (is.null(h) || !file.exists(p) || !identical(.nrv_sha(p), h)) stop("앞 단계 파일 hash 가 다르다: ", stage, "/", file)
  p
}
.nrv_rds <- function(a, stage, file) readRDS(nrv_stage_file(a$attempt, stage, file))

#' 저자 라벨은 author 단계에서만 읽는다
nrv_stage_reads_author <- function(stage, reader) if (identical(stage, "author")) reader$author() else invisible(NULL)

#' 단계 출력 기록기: 쓴 파일을 모아 manifest 에 넘긴다. 표는 mdu_write(%.17g), R 재검산용 사본은 RDS
.nrv_writer <- function(a, stage) {
  d <- nrv_stage_dir(a$attempt, stage); dir.create(d, recursive = TRUE, showWarnings = FALSE); files <- character(0)
  add <- function(p) { files <<- c(files, p); invisible(p) }
  list(dir = d, files = function() files, add = add,                                  # add: 다른 함수가 쓴 파일(그림, CSV)을 manifest 에 넣는다
       tsv = function(df, f) { p <- file.path(d, f); mdu_write(df, p); add(p) },
       rds = function(x, f) { p <- file.path(d, f); mdu_save_lossless(x, p); add(p) },
       json = function(x, f) { p <- file.path(d, f); writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), p); add(p) })
}

#' freeze 에서 고정한 인자와 지금 인자가 같은지 확인하고, 고정 내용을 돌려준다
.nrv_frozen <- function(a) {
  fz <- jsonlite::fromJSON(nrv_stage_file(a$attempt, "freeze", "plan_freeze.json"))
  keys <- c("partition", "allocation", "replicates", "n_null", "expected_blocks", "chunk", "seed")
  bad <- keys[!vapply(keys, function(k) identical(as.character(fz$args[[k]]), as.character(a[[k]])), logical(1))]
  if (length(bad)) stop("인자가 freeze 단계의 고정값과 다르다: ", paste(bad, collapse = ", "))
  fz
}


# ── 실제 입력 reader (읽기 전용) ─────────────────────────────────────────────
#' 테스트는 같은 모양의 mock 을 넣는다: ids, graph, partition, obs, symbols, raw_chunks, author, umap, eval_mapping (+ readout51, input_files)
nrv_h5_reader <- function(data_root) {
  P <- function(k) file.path(data_root, NRV_INPUTS[[k]]); h5 <- P("h5ad")
  with_h5 <- function(fun) { f <- hdf5r::H5File$new(h5, "r"); on.exit(f$close_all()); fun(f) }
  obs_ids <- function(f) as.character(f[[paste0("obs/", hdf5r::h5attr(f[["obs"]], "_index"))]]$read())
  list(ids = function() { g <- readLines(P("graph_ids")); g[nzchar(g)] },
       graph = function() ct_read_npz_sparse(P("graph")),
       partition = function(par) utils::read.delim(gzfile(P(par)), stringsAsFactors = FALSE, colClasses = c(cell = "character")),
       obs = function(name) with_h5(function(f) if (identical(name, "cell_id")) obs_ids(f) else ct_h5_obs(f, name)),
       symbols = function() with_h5(function(f) {
         gs <- f[["raw/var/gene_symbols"]]
         if (inherits(gs, "H5Group")) { cats <- as.character(gs[["categories"]]$read()); cats[as.integer(gs[["codes"]]$read()) + 1L] } else as.character(gs$read())
       }),
       raw_chunks = function(chunk) nrv_chunks_h5(h5, chunk),
       author = function() with_h5(function(f) {
         if (!f$exists("obs/AuthorCellType")) stop("obs/AuthorCellType 이 없다(cell_type 으로 대신하지 않는다).")
         ct_h5_obs(f, "AuthorCellType")
       }),
       umap = function() with_h5(function(f) ct_read_h5_matrix(f, "obsm/X_umap", length(obs_ids(f)))),
       eval_mapping = function() utils::read.delim(P("eval_mapping"), stringsAsFactors = FALSE),
       readout51 = function() utils::read.delim(gzfile(P("readout51")), stringsAsFactors = FALSE),
       input_files = function() vapply(names(NRV_INPUTS), P, character(1)))
}

#' raw chunk 검사를 끼운 반복자: 정수 count 가 아니거나 library 가 0 이면 멈춘다. 세포 수가 그래프와 다르면 멈춘다.
.nrv_checked_chunks <- function(each_chunk, n_cells) function(callback) {
  seen <- 0L
  each_chunk(function(X, rows) {
    X <- methods::as(X, "CsparseMatrix")
    if (length(X@x) && (any(!is.finite(X@x)) || any(X@x < 0) || any(X@x != round(X@x)))) stop("raw count 가 아니다(비유한, 음수, 정수 아님).")
    if (any(Matrix::rowSums(X) <= 0)) stop("library 합이 0 인 세포가 있다.")
    seen <<- seen + length(rows); callback(X, rows)
  })
  if (seen != n_cells) stop("raw 세포 수가 그래프 세포 수와 다르다: ", seen, " vs ", n_cells)
}

.nrv_set_matrix <- function(sets, n) {
  if (!length(sets)) return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, n)))
  nrv_set_matrix(sets, n)
}


# ── 실행 시작점 ──────────────────────────────────────────────────────────────
#' 인자 해석 → 단계 폴더 검사(6) → 소스 스냅샷 → 사본 적재 → 사본의 nrv_run_stage 실행. root / reader / after_snapshot 은 테스트용.
nrv_main <- function(argv, root = NULL, reader = NULL, after_snapshot = NULL) {
  a <- list(stage = NULL, attempt = NULL, data_root = NULL, partition = "p92", allocation = "all", replicates = 2000L, n_null = 1000L,
            expected_blocks = 42L, chunk = 5000L, seed = 20261013L)
  chr <- c("--stage" = "stage", "--attempt-dir" = "attempt", "--data-root" = "data_root", "--partition" = "partition", "--allocation" = "allocation")
  int <- c("--replicates" = "replicates", "--n-null" = "n_null", "--expected-blocks" = "expected_blocks", "--chunk" = "chunk", "--seed" = "seed")
  i <- 1L
  while (i <= length(argv)) {
    k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (k %in% names(chr)) a[[chr[[k]]]] <- v else if (k %in% names(int)) a[[int[[k]]]] <- suppressWarnings(as.integer(v)) else {
      cat("알 수 없는 인자:", k, "\n"); return(2L) }
    i <- i + 2L
  }
  if (is.null(a$stage) || is.na(a$stage) || !a$stage %in% NRV_STAGES || is.null(a$attempt) || is.null(a$data_root)) {
    cat("--stage / --attempt-dir / --data-root 가 필요하다.\n"); return(2L) }
  if (!a$partition %in% c("p92", "p51") || !a$allocation %in% c("all", names(NRV_ALLOCATIONS))) { cat("--partition 은 p92|p51, --allocation 은 all|A|B|C\n"); return(2L) }
  iv <- unlist(a[unname(int)]); if (anyNA(iv) || any(iv < 1L)) { cat("정수 인자는 양의 정수여야 한다.\n"); return(2L) }
  sd <- nrv_stage_dir(a$attempt, a$stage)
  if (dir.exists(sd) && length(list.files(sd, all.files = TRUE, no.. = TRUE))) { cat("단계 폴더가 비어 있지 않다:", sd, "\n"); return(6L) }
  dir.create(a$attempt, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    line <- paste(format(Sys.time(), "%H:%M:%S"), a$stage, ...); cat(line, "\n"); cat(line, "\n", file = file.path(a$attempt, "run.log"), append = TRUE)
  }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap_dir <- file.path(a$attempt, paste0("source_snapshot_", a$stage))
    if (dir.exists(snap_dir)) stop("소스 스냅샷 폴더가 이미 있다(새 attempt 폴더를 쓴다): ", snap_dir)
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", NRV_SOURCES), snap_dir)
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    if (is.function(after_snapshot)) after_snapshot()
    ld <- mdu_load_snapshot(snap, files = setdiff(NRV_SOURCES, "28_view_nrv.R"))
    log("loaded from snapshot:", paste(basename(ld$loaded$snapshot_path), collapse = ", "))
    body <- get("nrv_run_stage", envir = ld$env, inherits = FALSE)
    body(stage = a$stage, a = a, snap = snap, loaded = ld$loaded, reader = reader, log = log)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

#' 단계 본체 (스냅샷 사본에서 적재되어 실행된다)
nrv_run_stage <- function(stage, a, snap, loaded, reader, log) {
  up <- nrv_check_upstream(a$attempt, stage)
  if (is.null(reader)) reader <- nrv_h5_reader(a$data_root)
  if (stage %in% c("stream_global", "stream_sets", "stream_route", "author", "view")) .nrv_check_inputs(a, reader)   # 입력을 다시 여는 단계
  t0 <- Sys.time(); w <- .nrv_writer(a, stage)
  w$tsv(snap, "source_snapshot.tsv")
  fun <- switch(stage, freeze = nrv_stage_freeze, structure = nrv_stage_structure, stream_global = nrv_stage_stream_global,
                freeze_controls = nrv_stage_freeze_controls, stream_sets = nrv_stage_stream_sets, select = nrv_stage_select,
                stream_route = nrv_stage_stream_route, evaluate = nrv_stage_evaluate, verdict = nrv_stage_verdict,
                author = nrv_stage_author, view = nrv_stage_view)
  fun(a = a, reader = reader, log = log, w = w, snap = snap)
  nrv_write_stage_manifest(a$attempt, stage, files = w$files(), upstream = up)
  log("done:", length(w$files()), "files,", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
}

# ── 1. freeze: 역할 배정, 필수 프로그램, 선택 유전자, 기준 상수, 스트레스 점검 (자료를 읽지 않는다) ──────
nrv_stage_freeze <- function(a, reader, log, w, ...) {
  als <- if (a$allocation == "all") names(NRV_ALLOCATIONS) else a$allocation
  ref <- nrv_read_reference(NRV_REF_DIR); refg <- sort(unique(c(ref$gene_long$gene, ref$markers$gene)), method = "radix")
  by_al <- function(f) do.call(rbind, lapply(als, function(al) { x <- f(al); if (!nrow(x)) return(NULL); cbind(allocation = al, x, stringsAsFactors = FALSE) }))
  roles <- do.call(rbind, lapply(als, function(al) nrv_gene_roles(ref, al)))
  progs <- by_al(function(al) nrv_programs(nrv_gene_roles(ref, al)))
  req <- by_al(function(al) nrv_required_programs(nrv_programs(nrv_gene_roles(ref, al))))
  req$lit_rank <- as.integer(sub("^E", "", req$edge_id))                                     # 문헌 순서 = edge 순서 (같은 edge 안은 동률)
  req$in_early_path <- req$edge_id %in% NRV_EARLY_EDGES
  panels <- by_al(function(al) nrv_masked_panels(ref, al))
  selg <- by_al(function(al) do.call(rbind, lapply(c("anchor_route", "combined", "marker_only_full"), function(cd)
    data.frame(condition = cd, gene = nrv_selection_genes(ref, al, cd), stringsAsFactors = FALSE))))
  anchors <- by_al(function(al) { an <- nrv_anchor_genes(al)
    data.frame(anchor = rep(names(an), lengths(an)), gene = unlist(an, use.names = FALSE), stringsAsFactors = FALSE) })
  for (al in als) nrv_check_no_leak(roles[roles$allocation == al, ], progs$gene[progs$allocation == al],
                                    selg$gene[selg$allocation == al & selg$condition == "combined"])    # 실행 전 누설 검사
  ref_files <- sort(list.files(NRV_REF_DIR, full.names = TRUE), method = "radix")
  w$tsv(roles, "gene_roles.tsv"); w$tsv(progs, "programs.tsv"); w$tsv(req, "required_programs.tsv"); w$tsv(panels, "masked_panels.tsv")
  w$tsv(selg, "selection_genes.tsv"); w$tsv(anchors, "anchor_genes.tsv"); w$tsv(NRV_THRESHOLDS, "thresholds.tsv")
  w$tsv(nrv_stress_audit(), "stress_inflammation_audit.tsv")
  w$rds(list(ref = ref, refg = refg, als = als, roles = roles, progs = progs, req = req, panels = panels, selg = selg, anchors = anchors), "freeze_lossless.rds")
  w$json(list(args = a[c("partition", "allocation", "replicates", "n_null", "expected_blocks", "chunk", "seed")], allocations = als,
              reference_dir = NRV_REF_DIR, reference_sha256 = as.list(setNames(vapply(ref_files, .nrv_sha, ""), basename(ref_files))),
              route_kinds = as.list(NRV_ROUTE_KINDS), early_path_edges = NRV_EARLY_EDGES,
              rules = list(holm_families = "F1 adult_confirmed, F2 human_age_unconfirmed, F4 internal (키·배정·경로 종류 안); 유전자 F3 BY 는 서술",
                           gene_null_exclusion = "참고 37 유전자와 그 배정의 anchor 대조 유전자는 null 후보에서 뺀다",
                           internal_coordinate = "t_pooled (고정 전체 그래프), 보류 block 구조는 block_structure 로 따로",
                           single_gene = "유전자 1 개인 필수 프로그램(E07_core_up, E08_core_up = MMP9)은 단일 유전자 근거, 라벨 _single_gene_component",
                           stress_inflammation = "not_evaluated (승인된 근거 목록 없음)",
                           author_labels = "author 단계에서만 읽는다"),
              created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), "plan_freeze.json")
  log("frozen allocations:", paste(als, collapse = ","), "| required programs:", nrow(req))
}


# ── 2. structure: ID, 양수 edge, block, fold, 군집 쌍 표 (전체 + fold 학습) ────────────────────────────
#' block 정보: study 는 block 마다 하나여야 한다. sorting 구분은 Bulk → bulk, CD34+ 로 시작 → cd34_enriched (fold 층 나누기용, 운영 기준)
nrv_block_meta <- function(cells) {
  sc <- ifelse(cells$sorting == "Bulk", "bulk", ifelse(startsWith(cells$sorting, "CD34+"), "cd34_enriched", cells$sorting))
  sp <- split(seq_len(nrow(cells)), cells$block)
  do.call(rbind, lapply(sort(names(sp), method = "radix"), function(b) {
    i <- sp[[b]]; st <- unique(cells$study[i]); if (length(st) != 1L) stop("block 하나에 study 가 둘 이상이다: ", b)
    data.frame(block = b, study = st, sorting_class = paste(sort(unique(sc[i]), method = "radix"), collapse = "|"), n_cells = length(i),
               stringsAsFactors = FALSE)
  }))
}

#' study 하나 빼기: study 마다 fold 하나 (study 이름 순서)
nrv_loso_folds <- function(bm) {
  st <- sort(unique(bm$study), method = "radix")
  data.frame(block = bm$block, fold = match(bm$study, st), held_out_study = bm$study, stringsAsFactors = FALSE)
}

nrv_stage_structure <- function(a, reader, log, w, ...) {
  .nrv_frozen(a)
  ids <- as.character(reader$ids()); ct_check_ids(ids, reader$obs("cell_id"))
  part <- reader$partition(a$partition); ct_check_ids(ids, part$cell)
  G <- reader$graph(); if (nrow(G) != length(ids) || ncol(G) != length(ids)) stop("그래프 크기가 세포 수와 다르다: ", nrow(G), " vs ", length(ids))
  e <- nrv_positive_upper_edges(G); dropped <- attr(e, "dropped"); rm(G); log("positive upper edges:", nrow(e))
  block <- nrv_blocks(reader$obs("donor_id"), expected_n = a$expected_blocks)
  cells <- data.frame(cell = seq_along(ids), cell_id = ids, cluster = as.character(part$unit), block = block,
                      study = as.character(reader$obs("Study")), sorting = as.character(reader$obs("Sorting")), stringsAsFactors = FALSE)
  if (anyNA(cells$cluster) || anyNA(cells$study) || anyNA(cells$sorting)) stop("군집, Study, Sorting 에 결측이 있다.")
  bm <- nrv_block_meta(cells)
  folds <- rbind(cbind(holdout_scheme = "block5", nrv_folds(bm, k = .nrv_thr("n_folds"), seed = a$seed), held_out_study = NA_character_,
                       stringsAsFactors = FALSE),
                 cbind(holdout_scheme = "loso", nrv_loso_folds(bm), stringsAsFactors = FALSE))
  keys <- unique(folds[, c("holdout_scheme", "fold")]); rownames(keys) <- NULL
  pt_full <- nrv_pair_table(e, cells$cluster, cells$block, cells$study)
  pt_fold <- lapply(seq_len(nrow(keys)), function(r) {
    ho <- folds$block[folds$holdout_scheme == keys$holdout_scheme[r] & folds$fold == keys$fold[r]]
    nrv_pair_table(e, cells$cluster, cells$block, cells$study, cell_mask = !(cells$block %in% ho))      # 양끝이 모두 학습 세포인 edge 만
  })
  names(pt_fold) <- paste(keys$holdout_scheme, keys$fold, sep = "|")
  log("blocks:", nrow(bm), "| clusters:", length(unique(cells$cluster)), "| folds:", nrow(keys), "| observed pairs:", nrow(pt_full))
  if (is.function(reader$input_files)) {
    f <- reader$input_files(); f <- f[file.exists(f)]
    w$tsv(data.frame(input = names(f), path = unname(f), sha256 = vapply(f, .nrv_sha, "", USE.NAMES = FALSE), stringsAsFactors = FALSE), "input_sha256.tsv")
  }
  w$tsv(data.frame(n_cells = nrow(cells), n_clusters = length(unique(cells$cluster)), n_blocks = nrow(bm), n_positive_upper_edges = nrow(e),
                   n_dropped_diagonal = unname(dropped["diagonal"]), n_dropped_stored_zero = unname(dropped["stored_zero"]),
                   n_cross_cluster_edges = sum(cells$cluster[e$i] != cells$cluster[e$j]), n_observed_pairs = nrow(pt_full)), "edge_audit.tsv")
  w$tsv(cells, "cells.tsv.gz"); w$tsv(bm, "block_meta.tsv"); w$tsv(folds, "folds.tsv"); w$tsv(pt_full, "pairs_full.tsv.gz")
  w$tsv(do.call(rbind, lapply(names(pt_fold), function(k) cbind(key = k, pt_fold[[k]], stringsAsFactors = FALSE))), "pairs_folds.tsv.gz")
  w$rds(e, "edges.rds"); w$rds(list(cells = cells, bm = bm, folds = folds, keys = keys), "cells_lossless.rds")
  w$rds(list(full = pt_full, folds = pt_fold), "pairs.rds")
}


# ── 3. stream_global: 전체 유전자 전역 평균·검출률 (raw 1 회) ──────────────────────────────────────────
nrv_stage_stream_global <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  sym <- as.character(reader$symbols()); if (!length(sym) || anyNA(sym)) stop("raw 유전자 symbol 이 없거나 결측이 있다.")
  gs <- nrv_gene_global_stats(.nrv_checked_chunks(reader$raw_chunks(a$chunk), nrow(cl$cells)), sym)
  log("genes:", nrow(gs), "| expressed (mean_log > 0):", sum(gs$mean_log > 0))
  w$tsv(gs, "gene_stats.tsv.gz"); w$rds(list(gs = gs, symbols = sym), "gene_stats_lossless.rds")
}


# ── 4. freeze_controls: anchor·node 대조 유전자, gene null 집합 (발현 결과를 보기 전) ──────────────────
#' 대조를 붙일 module 목록: 배정별 anchor(S, T_early, T_terminal), 배정별 masked node panel(비지 않은 것), 공통 full node panel
.nrv_modules <- function(fr) {
  an <- fr$anchors; u <- unique(an[, c("allocation", "anchor")])
  m1 <- data.frame(module_id = paste(u$allocation, "anchor", u$anchor, sep = ":"), allocation = u$allocation, kind = "anchor", name = u$anchor,
                   genes = vapply(seq_len(nrow(u)), function(r) paste(an$gene[an$allocation == u$allocation[r] & an$anchor == u$anchor[r]], collapse = "|"), ""),
                   stringsAsFactors = FALSE)
  pn <- fr$panels; ok <- pn$status == "ok"
  m2 <- data.frame(module_id = paste(pn$allocation[ok], "masked", pn$node_id[ok], sep = ":"), allocation = pn$allocation[ok], kind = "masked",
                   name = pn$node_id[ok], genes = pn$masked_panel[ok], stringsAsFactors = FALSE)
  fp <- unique(pn[, c("node_id", "full_panel")])
  m3 <- data.frame(module_id = paste("all", "full", fp$node_id, sep = ":"), allocation = "all", kind = "full", name = fp$node_id, genes = fp$full_panel,
                   stringsAsFactors = FALSE)
  rbind(m1, m2, m3)
}

nrv_stage_freeze_controls <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); gs <- .nrv_rds(a, "stream_global", "gene_stats_lossless.rds")$gs
  ref_cols <- gs$col[gs$symbol %in% fr$refg]; cols_of <- function(g) gs$col[gs$symbol %in% g]
  mods <- .nrv_modules(fr)
  ctl <- do.call(rbind, lapply(seq_len(nrow(mods)), function(r) {
    Q <- cols_of(strsplit(mods$genes[r], "|", fixed = TRUE)[[1]]); if (!length(Q)) return(NULL)
    cbind(module_id = mods$module_id[r], nrv_select_controls(gs, Q, exclude_cols = ref_cols, seed = a$seed), stringsAsFactors = FALSE)
  }))
  mods$n_panel_cols <- vapply(mods$genes, function(g) length(cols_of(strsplit(g, "|", fixed = TRUE)[[1]])), integer(1), USE.NAMES = FALSE)
  mods$n_controls <- as.integer(table(factor(ctl$module_id, levels = mods$module_id)))
  null <- lapply(fr$als, function(al) {
    pg <- fr$progs[fr$progs$allocation == al & fr$progs$program_id %in% fr$req$program_id[fr$req$allocation == al], ]
    pg$col <- vapply(pg$gene, function(g) { k <- cols_of(g); if (length(k)) as.integer(min(k)) else NA_integer_ }, integer(1), USE.NAMES = FALSE)
    pg$n_cols <- vapply(pg$gene, function(g) length(cols_of(g)), integer(1), USE.NAMES = FALSE)
    tg <- pg[!is.na(pg$col), ]; rownames(tg) <- NULL
    excl <- sort(unique(c(ref_cols, ctl$col[ctl$module_id %in% mods$module_id[mods$allocation == al & mods$kind == "anchor"]])))
    N <- nrv_null_sets(gs, tg$col, excl, R = a$n_null, n_bins = 10L, seed = a$seed)
    list(targets = tg, missing = pg[is.na(pg$col), ], N = N, exclude_cols = excl)
  })
  names(null) <- fr$als
  pool <- sort(unique(unlist(lapply(null, function(x) as.integer(x$N)))))
  log("modules:", nrow(mods), "| control columns:", length(unique(ctl$col)), "| null pool columns:", length(pool))
  w$tsv(mods, "modules.tsv"); w$tsv(ctl, "controls.tsv")
  w$tsv(do.call(rbind, lapply(names(null), function(al) null[[al]]$targets)), "null_targets.tsv")
  w$tsv(data.frame(col = pool, symbol = gs$symbol[match(pool, gs$col)], stringsAsFactors = FALSE), "null_pool.tsv")
  w$rds(list(modules = mods, controls = ctl, null = null, pool = pool), "controls_lossless.rds")
}


# ── 5. stream_sets: 군집 x block 합(참고 + 대조), 같은 block 경계 합과 전체 경계 합(참고) (raw 1 회) ─────
nrv_stage_stream_sets <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  e <- .nrv_rds(a, "structure", "edges.rds"); ct <- .nrv_rds(a, "freeze_controls", "controls_lossless.rds")
  sym <- as.character(reader$symbols()); n <- nrow(cl$cells)
  map_ref <- nrv_gene_map(sym, fr$refg); ctl_cols <- sort(unique(ct$controls$col)); map_ctl <- setNames(as.list(ctl_cols), paste0("col", ctl_cols))
  sd <- nrv_set_definitions(e, cl$cells); H <- lapply(sd, function(x) .nrv_set_matrix(x$sets, n))
  log("sets: cb", nrow(H$cb), "| sb", nrow(H$sb), "| pool", nrow(H$pool), "| genes ref", length(map_ref), "ctl", length(map_ctl))
  out <- nrv_stream_multi(.nrv_checked_chunks(reader$raw_chunks(a$chunk), n), maps = list(ref = map_ref, ctl = map_ctl),
                          targets = list(cb_ref = list(H = H$cb, map = "ref", det = FALSE), cb_ctl = list(H = H$cb, map = "ctl", det = FALSE),
                                         sb = list(H = H$sb, map = "ref", det = TRUE), pool = list(H = H$pool, map = "ref", det = TRUE)))
  w$tsv(data.frame(kind = names(sd), n_sets = vapply(sd, function(x) length(x$sets), integer(1)),
                   n_memberships = vapply(sd, function(x) sum(lengths(x$sets)), integer(1)), stringsAsFactors = FALSE), "set_summary.tsv")
  w$tsv(data.frame(gene = fr$refg, in_data = !fr$refg %in% attr(map_ref, "missing"), n_cols = unname(lengths(map_ref)[fr$refg]),
                   stringsAsFactors = FALSE), "reference_gene_presence.tsv")
  w$rds(list(meta = lapply(sd, `[[`, "meta"), sums = out, missing_ref_genes = attr(map_ref, "missing")), "set_sums.rds")
}


# ── 공통: 경로 구간 배열, 군집 x block 보기, module 점수, 경로 종류의 필수 프로그램 ──────────────────────

#' 경로 구간 배열 (구간 x block x 유전자). 경로 방향 i → j 의 Δμ = 평균 x(B^b_{j→i}) − 평균 x(B^b_{i→j}) (하류 − 상류, E3).
#' det = 두 경계의 검출 세포 합. 같은 block 경계 집합이 없거나 한쪽 세포가 min_cells 보다 적은 block 은 Δμ = NA (det 는 있는 만큼 센다).
#' sums 의 행 순서는 meta$set_id 와 같아야 한다.
nrv_segment_arrays <- function(sums, meta, path, blocks, genes, min_cells = 5) {
  seg <- paste(path[-length(path)], path[-1], sep = ">")
  D <- array(NA_real_, c(length(seg), length(blocks), length(genes)), dimnames = list(seg, blocks, genes)); E <- D; E[] <- 0
  row_of <- function(f, t) match(paste("sb", paste(f, t, sep = ">"), blocks, sep = "|"), meta$set_id)
  for (s in seq_along(seg)) {
    rf <- row_of(path[s], path[s + 1L]); rt <- row_of(path[s + 1L], path[s])            # from 쪽 경계 B_{i→j}, to 쪽 경계 B_{j→i}
    nf <- ifelse(is.na(rf), 0, sums$n[rf]); nt <- ifelse(is.na(rt), 0, sums$n[rt]); ok <- nf >= min_cells & nt >= min_cells
    if (any(ok)) D[s, ok, ] <- sums$sum[rt[ok], genes, drop = FALSE] / nt[ok] - sums$sum[rf[ok], genes, drop = FALSE] / nf[ok]
    df <- matrix(0, length(blocks), length(genes)); dt <- df
    if (any(!is.na(rf))) df[!is.na(rf), ] <- sums$det[rf[!is.na(rf)], genes, drop = FALSE]
    if (any(!is.na(rt))) dt[!is.na(rt), ] <- sums$det[rt[!is.na(rt)], genes, drop = FALSE]
    E[s, , ] <- df + dt
  }
  list(delta = D, det = E)
}

#' 군집 x block 보기: 세포 수 행렬과, 유전자(참고 symbol 또는 대조 "col<번호>") 의 군집 x block 평균 행렬
.nrv_cb_view <- function(ss) {
  m <- ss$meta$cb; cl <- .nrv_lv(m$cluster); bl <- .nrv_lv(m$block); pos <- cbind(match(m$cluster, cl), match(m$block, bl))
  n <- matrix(0, length(cl), length(bl), dimnames = list(cl, bl)); n[pos] <- ss$sums$cb_ref$n
  have <- c(colnames(ss$sums$cb_ref$sum), colnames(ss$sums$cb_ctl$sum))
  mean_of <- function(g) {
    S <- if (g %in% colnames(ss$sums$cb_ref$sum)) ss$sums$cb_ref$sum else ss$sums$cb_ctl$sum
    M <- matrix(NA_real_, length(cl), length(bl), dimnames = dimnames(n)); M[pos] <- S[, g] / ss$sums$cb_ref$n; M
  }
  list(n = n, mean_of = mean_of, has = function(g) g %in% have)
}

#' module 점수 (군집 x block) = 패널 유전자 평균의 평균 − 대조 유전자 평균의 평균 (세포 점수 "패널 평균 − 대조 평균" 의 평균과 같다)
.nrv_module_score <- function(cbv, genes, ctl_cols) {
  g <- genes[vapply(genes, cbv$has, logical(1))]; k <- paste0("col", ctl_cols)
  if (!length(g) || !length(k)) return(NULL)
  Reduce(`+`, lapply(g, cbv$mean_of)) / length(g) - Reduce(`+`, lapply(k, cbv$mean_of)) / length(k)
}

#' 경로 종류의 필수 프로그램 정의: 자료에 있는 유전자만 쓰고, gene null 대상 행(null_rows)과 같은 순서로 둔다
.nrv_route_programs <- function(fr, ct, al, kind) {
  rq <- fr$req[fr$req$allocation == al, ]; if (kind == "partial_path_early_granule") rq <- rq[rq$in_early_path, ]
  tg <- ct$null[[al]]$targets
  out <- lapply(seq_len(nrow(rq)), function(k) {
    r <- which(tg$program_id == rq$program_id[k]); by_rule <- isTRUE(rq$program_evaluable_by_rule[k])
    list(program_id = rq$program_id[k], edge_id = rq$edge_id[k], layer = rq$layer[k], lit_rank = rq$lit_rank[k], genes = tg$gene[r], signs = tg$sign[r],
         null_rows = r, single_gene = !by_rule, min_genes = if (by_rule) .nrv_thr("program_min_genes") else 1,
         conditional = grepl("_conditional_", rq$program_id[k], fixed = TRUE))
  })
  setNames(out, rq$program_id)
}

.nrv_key_cols <- function(a, en) data.frame(partition = a$partition, holdout_scheme = en$holdout_scheme, fold = en$fold, allocation = en$allocation,
                                            stringsAsFactors = FALSE)

.nrv_route_row <- function(rt) {
  s <- function(x) if (is.null(x) || !length(x)) "" else paste(x, collapse = ">")
  data.frame(structural_path = rt$status, primary_path = s(rt$primary), primary_cost = if (is.null(rt$primary_cost)) NA_real_ else rt$primary_cost,
             n_tied_primary = if (is.null(rt$n_tied_primary)) NA_integer_ else as.integer(rt$n_tied_primary),
             overlap = paste(rt$overlap, collapse = "|"), nontrivial_status = if (is.null(rt$nontrivial)) "" else rt$nontrivial$status,
             nontrivial_path = s(rt$nontrivial$primary), n_alternatives = if (is.null(rt$alternatives)) 0L else nrow(rt$alternatives),
             stringsAsFactors = FALSE)
}


# ── 6. select: anchor 후보, 경로(유전자 없음), marker 라벨, 구간 σ (모두 학습 block 만), 내부 좌표 ─────────
nrv_stage_select <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  e <- .nrv_rds(a, "structure", "edges.rds"); pr <- .nrv_rds(a, "structure", "pairs.rds")
  ct <- .nrv_rds(a, "freeze_controls", "controls_lossless.rds"); ss <- .nrv_rds(a, "stream_sets", "set_sums.rds")
  cbv <- .nrv_cb_view(ss); blocks <- cl$bm$block; mods <- ct$modules
  scores <- setNames(lapply(mods$module_id, function(id) .nrv_module_score(cbv, strsplit(mods$genes[mods$module_id == id], "|", fixed = TRUE)[[1]],
                                                                         ct$controls$col[ct$controls$module_id == id])), mods$module_id)
  freq <- .nrv_thr("anchor_candidate_freq"); minc <- .nrv_thr("anchor_block_min_cells"); minb <- .nrv_thr("anchor_min_blocks")
  marker <- function(ids, nodes, train) {
    sc <- setNames(lapply(ids, function(id) { s <- scores[[id]]
      if (is.null(s)) matrix(NA_real_, nrow(cbv$n), length(train), dimnames = list(rownames(cbv$n), train)) else s[, train, drop = FALSE] }), nodes)
    nrv_marker_only(sc, cbv$n[, train, drop = FALSE], R = a$replicates, seed = a$seed, min_cells = minc)
  }
  entries <- list()
  for (r in seq_len(nrow(cl$keys))) for (al in fr$als) {
    sch <- cl$keys$holdout_scheme[r]; fo <- cl$keys$fold[r]; key <- paste(sch, fo, sep = "|")
    held <- cl$folds$block[cl$folds$holdout_scheme == sch & cl$folds$fold == fo]; train <- setdiff(blocks, held)
    atab <- do.call(rbind, lapply(c("S", "T_early", "T_terminal"), function(an) {
      s <- scores[[paste(al, "anchor", an, sep = ":")]]
      if (is.null(s)) return(data.frame(anchor = an, cluster = NA_character_, score = NA_real_, n_eligible_blocks = 0L, top_freq = NA_real_,
                                        n_draws_no_candidate = NA_integer_, anchor_status = "anchor_unavailable", candidate = FALSE, stringsAsFactors = FALSE))
      x <- nrv_anchor_scores(s[, train, drop = FALSE], cbv$n[, train, drop = FALSE], min_cells = minc, min_blocks = minb, R = a$replicates, seed = a$seed)
      x$candidate <- x$anchor_status == "candidate_only" & !is.na(x$top_freq) & x$top_freq >= freq - MDU_TIE_TOL   # 빈도 0.10 (반올림 잡음만 허용)
      cbind(anchor = an, x, stringsAsFactors = FALSE)
    }))
    cand <- function(an) atab$cluster[atab$anchor == an & atab$candidate]
    routes <- lapply(NRV_ROUTE_KINDS, function(tk) nrv_route(pr$folds[[key]], cand("S"), cand(tk), k_alt = .nrv_thr("yen_k")))
    pn <- fr$panels[fr$panels$allocation == al, ]
    mo_masked <- marker(paste(al, "masked", pn$node_id, sep = ":"), pn$node_id, train)
    mo_full <- marker(paste("all", "full", pn$node_id, sep = ":"), pn$node_id, train)
    sig <- do.call(rbind, lapply(names(NRV_ROUTE_KINDS), function(kind) {
      rt <- routes[[kind]]; if (!identical(rt$status, "route_ok")) return(NULL)
      P <- .nrv_route_programs(fr, ct, al, kind)
      arr <- nrv_segment_arrays(ss$sums$sb, ss$meta$sb, rt$primary, blocks, unique(unlist(lapply(P, `[[`, "genes"))), .nrv_thr("set_min_cells"))
      do.call(rbind, lapply(P, function(p) {
        s <- nrv_route_program_stats(arr, p$genes, p$signs, train, held, p$min_genes, .nrv_thr("gene_min_detected_cells"))
        tv <- if (is.na(s$sigma)) NA_real_ else s$train_seg[[s$sigma]]
        data.frame(route_kind = kind, program_id = p$program_id, sigma = s$sigma, sigma_index = match(s$sigma, names(s$train_seg)), train_value = tv,
                   n_genes_at_sigma = if (is.na(s$sigma)) 0L else s$n_genes[[s$sigma]],
                   n_tied_at_max = if (is.na(tv)) 0L else sum(abs(s$train_seg - tv) <= MDU_TIE_TOL, na.rm = TRUE), stringsAsFactors = FALSE)
      }))
    }))
    tri <- do.call(rbind, lapply(names(NRV_ROUTE_KINDS), function(kind) {
      rt <- routes[[kind]]; p <- rt$primary; m <- length(p)
      if (!identical(rt$status, "route_ok") || m < 3L) return(NULL)
      data.frame(route_kind = kind, position = 2:(m - 1L), B = p[2:(m - 1L)], A = p[1:(m - 2L)], C = p[3:m], stringsAsFactors = FALSE)
    }))
    entries[[paste(key, al, sep = "|")]] <- list(holdout_scheme = sch, fold = fo, key = key, allocation = al, train = train, held = held, anchors = atab,
                                                routes = routes, marker_masked = mo_masked, marker_full = mo_full, sigma = sig, interior = tri)
  }
  # 내부 좌표 (그래프만, 고정 전체 그래프): A, B, C 사이 edge 만 남겨도 결과가 같다
  tri_all <- unique(do.call(rbind, c(list(data.frame(B = character(0), A = character(0), C = character(0), stringsAsFactors = FALSE)),
                                     lapply(entries, function(x) x$interior[, c("B", "A", "C")]))))
  ci <- cl$cells$cluster[e$i]; cj <- cl$cells$cluster[e$j]; icoord <- list()
  for (r in seq_len(nrow(tri_all))) {
    tr <- tri_all[r, ]; abc <- c(tr$A, tr$B, tr$C); es <- e[ci %in% abc & cj %in% abc, , drop = FALSE]
    ic <- nrv_internal_coordinate(es, cl$cells$cluster, tr$B, tr$A, tr$C)
    icoord[[paste(tr$B, tr$A, tr$C, sep = "|")]] <- list(B = tr$B, A = tr$A, C = tr$C, status = ic$status, cells = ic$cells, n_shared = ic$n_shared,
      n_unreachable = ic$n_unreachable, block_structure = nrv_block_structure(es, cl$cells$cluster, cl$cells$block, tr$B, tr$A, tr$C))
  }
  log("entries:", length(entries), "| route_ok full:", sum(vapply(entries, function(x) identical(x$routes$full_route$status, "route_ok"), logical(1))),
      "| interior triples:", length(icoord))
  kc <- function(en) .nrv_key_cols(a, en)
  w$tsv(do.call(rbind, lapply(entries, function(en) cbind(kc(en), en$anchors, stringsAsFactors = FALSE))), "anchors.tsv")
  w$tsv(do.call(rbind, lapply(entries, function(en) do.call(rbind, lapply(names(en$routes), function(k) cbind(kc(en), route_kind = k,
    S_candidates = paste(en$anchors$cluster[en$anchors$anchor == "S" & en$anchors$candidate], collapse = "|"),
    T_candidates = paste(en$anchors$cluster[en$anchors$anchor == NRV_ROUTE_KINDS[[k]] & en$anchors$candidate], collapse = "|"),
    .nrv_route_row(en$routes[[k]]), stringsAsFactors = FALSE))))), "routes.tsv")
  alt <- do.call(rbind, lapply(entries, function(en) do.call(rbind, lapply(names(en$routes), function(k) {
    x <- en$routes[[k]]$alternatives; if (is.null(x) || !nrow(x)) NULL else cbind(kc(en), route_kind = k, x, stringsAsFactors = FALSE) }))))
  if (!is.null(alt)) w$tsv(alt, "route_alternatives.tsv")
  w$tsv(do.call(rbind, lapply(entries, function(en) rbind(cbind(kc(en), panel = "masked", en$marker_masked, stringsAsFactors = FALSE),
                                                         cbind(kc(en), panel = "full_circular", en$marker_full, stringsAsFactors = FALSE)))), "marker_labels.tsv")
  sg <- do.call(rbind, lapply(entries, function(en) if (is.null(en$sigma)) NULL else cbind(kc(en), en$sigma, stringsAsFactors = FALSE)))
  if (!is.null(sg)) w$tsv(sg, "sigma.tsv")
  if (length(icoord)) {
    w$tsv(do.call(rbind, lapply(names(icoord), function(id) { x <- icoord[[id]]
      data.frame(triple = id, B = x$B, A = x$A, C = x$C, pooled_structure = x$status, n_cells_B = if (is.null(x$cells)) 0L else nrow(x$cells),
                 n_t_finite = if (is.null(x$cells)) 0L else sum(is.finite(x$cells$t)), n_shared = x$n_shared, n_unreachable = x$n_unreachable,
                 stringsAsFactors = FALSE) })), "internal_structure.tsv")
    w$tsv(do.call(rbind, lapply(names(icoord), function(id) cbind(triple = id, icoord[[id]]$block_structure, stringsAsFactors = FALSE))), "block_structure.tsv")
  }
  w$rds(list(entries = entries, icoord = icoord), "selection.rds")
}


# ── 7. stream_route: 경로 구간 경계의 null 풀 합, 내부 군집 가중 집합 합 (raw 1 회) ─────────────────────
#' 내부 군집 가중 집합: (B|A|C, block) 마다 t 가 유한한 B 세포. H1 = 지시(Σx), Ht = 가중치 t(Σ t·x). 세 무리(서술용) 집합도 만든다.
#' 세포 수준 행렬을 저장하지 않고 기울기를 nrv_slope_from_sums 로 계산하기 위한 것이다(세포 수준 OLS 와 같다, 테스트로 고정).
.nrv_interior_sets <- function(icoord, cells, n) {
  empty <- function() Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0L, n))
  df <- do.call(rbind, lapply(names(icoord), function(id) {
    cc <- icoord[[id]]$cells; if (is.null(cc)) return(NULL); cc <- cc[is.finite(cc$t), c("cell", "t")]; if (!nrow(cc)) return(NULL)
    data.frame(triple = id, cell = cc$cell, t = cc$t, block = cells$block[cc$cell], stringsAsFactors = FALSE)
  }))
  if (is.null(df)) return(list(H1 = empty(), Ht = empty(), Hg = empty(),
                               meta = data.frame(set_id = character(0), triple = character(0), block = character(0), n = integer(0), t1 = numeric(0), t2 = numeric(0)),
                               gmeta = data.frame(set_id = character(0), triple = character(0), block = character(0), group = character(0), n = integer(0))))
  sid <- paste("int", df$triple, df$block, sep = "|"); ids <- sort(unique(sid), method = "radix"); r <- match(sid, ids); f1 <- match(ids, sid)
  H1 <- Matrix::sparseMatrix(i = r, j = df$cell, x = 1, dims = c(length(ids), n), dimnames = list(ids, NULL))
  Ht <- Matrix::sparseMatrix(i = r, j = df$cell, x = df$t, dims = c(length(ids), n), dimnames = list(ids, NULL))
  meta <- data.frame(set_id = ids, triple = df$triple[f1], block = df$block[f1], n = tabulate(r, length(ids)),
                     t1 = as.numeric(rowsum(df$t, r)), t2 = as.numeric(rowsum(df$t^2, r)), stringsAsFactors = FALSE)
  grp <- as.character(cut(df$t, c(-Inf, 1 / 3, 2 / 3, Inf), labels = c("A_side", "middle", "C_side"), right = FALSE))
  gid <- paste("grp", df$triple, df$block, grp, sep = "|"); gids <- sort(unique(gid), method = "radix"); rg <- match(gid, gids); fg <- match(gids, gid)
  Hg <- Matrix::sparseMatrix(i = rg, j = df$cell, x = 1, dims = c(length(gids), n), dimnames = list(gids, NULL))
  gmeta <- data.frame(set_id = gids, triple = df$triple[fg], block = df$block[fg], group = grp[fg], n = tabulate(rg, length(gids)), stringsAsFactors = FALSE)
  list(H1 = H1, Ht = Ht, Hg = Hg, meta = meta, gmeta = gmeta)
}

nrv_stage_stream_route <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  e <- .nrv_rds(a, "structure", "edges.rds"); ct <- .nrv_rds(a, "freeze_controls", "controls_lossless.rds"); sel <- .nrv_rds(a, "select", "selection.rds")
  n <- nrow(cl$cells); sym <- as.character(reader$symbols()); sd <- nrv_set_definitions(e, cl$cells)
  seg <- do.call(rbind, lapply(sel$entries, function(x) do.call(rbind, lapply(x$routes, function(rt) {
    if (!identical(rt$status, "route_ok")) return(NULL); p <- rt$primary; k <- length(p)
    data.frame(from = c(p[-k], p[-1]), to = c(p[-1], p[-k]), stringsAsFactors = FALSE) }))))
  m <- sd$sb$meta; keep <- if (is.null(seg)) logical(nrow(m)) else paste(m$from, m$to) %in% paste(seg$from, seg$to)
  route_meta <- m[keep, , drop = FALSE]; rownames(route_meta) <- NULL; H_route <- .nrv_set_matrix(sd$sb$sets[route_meta$set_id], n)
  int <- .nrv_interior_sets(sel$icoord, cl$cells, n)
  map_pool <- setNames(as.list(ct$pool), paste0("col", ct$pool)); map_ref <- nrv_gene_map(sym, fr$refg)
  log("route sets:", nrow(H_route), "| interior sets:", nrow(int$H1), "| group sets:", nrow(int$Hg), "| null pool:", length(map_pool))
  out <- nrv_stream_multi(.nrv_checked_chunks(reader$raw_chunks(a$chunk), n), maps = list(ref = map_ref, pool = map_pool),
                          targets = list(route_pool = list(H = H_route, map = "pool", det = TRUE), int1_pool = list(H = int$H1, map = "pool", det = FALSE),
                                         intt_pool = list(H = int$Ht, map = "pool", det = FALSE), int1_ref = list(H = int$H1, map = "ref", det = FALSE),
                                         intt_ref = list(H = int$Ht, map = "ref", det = FALSE), grp_ref = list(H = int$Hg, map = "ref", det = FALSE)))
  w$tsv(route_meta, "route_sets.tsv"); w$tsv(int$meta, "interior_sets.tsv"); w$tsv(int$gmeta, "interior_group_sets.tsv")
  w$rds(list(route_meta = route_meta, int_meta = int$meta, grp_meta = int$gmeta, pool = ct$pool, sums = out), "route_sums.rds")
}


# ── 8. evaluate: 보류 block 검정, 선택을 반복한 gene null, 내부 연속성, 파동 서술, 누설 검사 ─────────────
#' 구간 칸의 관측 서술 (보류 block, 보정 없는 부호 검정과 CI). 값은 문헌 방향 부호를 곱한 Δ 이므로 increase = 문헌 방향 변화다.
.nrv_desc <- function(h, alpha) {
  if (is.null(h) || h$n_eligible < 2L || is.na(h$sign_test$p_exact)) return("insufficient")
  if (h$sign_test$p_exact <= alpha && !is.na(h$ci_low) && h$ci_low > 0) return("observed_increase")
  if (h$sign_test$p_exact <= alpha && !is.na(h$ci_high) && h$ci_high < 0) return("observed_decrease")
  "plateau_or_not_significant"
}

#' 프로그램 판정 (E6): Holm 보정 부호 검정 p ≤ α, CI 가 문헌 방향에서 0 을 벗어남, 선택을 반복한 gene null p ≤ α 를 모두 만족해야 한다
.nrv_verdict <- function(sigma, p_holm, ci_low, ci_high, p_up, p_down, conditional, alpha) {
  if (is.na(sigma)) return(c("insufficient_evidence", "no_evaluable_genes"))
  if (is.na(p_holm)) return(c("insufficient_evidence", "holdout_unavailable"))
  up <- p_holm <= alpha && !is.na(ci_low) && ci_low > 0; dn <- p_holm <= alpha && !is.na(ci_high) && ci_high < 0
  if (up && !is.na(p_up) && p_up <= alpha) return(c(if (isTRUE(conditional)) "supported_conditional" else "supported", ""))
  if (dn && !is.na(p_down) && p_down <= alpha) return(c("conflicting", ""))
  c("insufficient_evidence", if (p_holm > alpha) "not_significant_after_holm" else if (!(up || dn)) "ci_includes_zero" else "gene_null_not_significant")
}

#' 내부 판정 (E4): pooled 구조가 span_ok 일 때만 평가한다. 기준은 프로그램 판정과 같은 세 가지다(Holm F4).
.nrv_internal_status <- function(pooled, selected, p_holm, ci_low, ci_high, p_up, p_down, alpha) {
  if (!identical(pooled, "span_ok")) return(c("insufficient_internal", paste0("pooled_", pooled)))
  if (is.na(selected)) return(c("insufficient_internal", "no_eligible_training_blocks"))
  if (is.na(p_holm)) return(c("insufficient_internal", "holdout_unavailable"))
  up <- p_holm <= alpha && !is.na(ci_low) && ci_low > 0; dn <- p_holm <= alpha && !is.na(ci_high) && ci_high < 0
  if (up && !is.na(p_up) && p_up <= alpha) return(c("internally_continuous", ""))
  if (dn && !is.na(p_down) && p_down <= alpha) return(c("internally_reversed", ""))
  c("not_supported_internal", if (p_holm > alpha) "not_significant_after_holm" else if (!(up || dn)) "ci_includes_zero" else "gene_null_not_significant")
}

#' 프로그램마다: 관측 통계, 보류 block 추론, 선택을 반복한 gene null, 구간 칸 서술, σ 의 유전자별 근거
.nrv_eval_programs <- function(P, arr, arrp, N, en, a, mind, alpha) {
  seg <- dimnames(arr$delta)[[1]]; pr <- list(); sr <- list(); gr <- list(); nulls <- list()
  inf <- function(x) nrv_block_inference(x, R = a$replicates, seed = a$seed)
  for (p in P) {
    so <- nrv_route_program_stats(arr, p$genes, p$signs, en$train, en$held, p$min_genes, mind); h <- inf(so$heldout); st <- h$sign_test
    ne <- vapply(seq_len(nrow(N)), function(r) nrv_route_program_stats(arrp, paste0("col", N[r, p$null_rows]), p$signs, en$train, en$held,
                                                                      p$min_genes, mind)$estimate, numeric(1))      # 같은 평가 가능성·σ 선택 반복
    pu <- nrv_gene_null_p(so$estimate, ne); pd <- nrv_gene_null_p(-so$estimate, -ne); nulls[[p$program_id]] <- ne
    pr[[p$program_id]] <- data.frame(program_id = p$program_id, edge_id = p$edge_id, layer = p$layer, lit_rank = p$lit_rank, conditional = p$conditional,
      single_gene = p$single_gene, n_genes_program = length(p$genes), sigma = so$sigma, sigma_index = match(so$sigma, seg),
      train_value = if (is.na(so$sigma)) NA_real_ else so$train_seg[[so$sigma]], n_genes_at_sigma = if (is.na(so$sigma)) 0L else so$n_genes[[so$sigma]],
      n_held_eligible = h$n_eligible, estimate = h$estimate, ci_low = h$ci_low, ci_high = h$ci_high, n_pos = st$n_pos, n_neg = st$n_neg, n_zero = st$n_zero,
      sign_test_status = st$status, p_sign = st$p_exact, p_null_up = as.numeric(pu), p_null_down = as.numeric(pd), n_null_valid = sum(!is.na(ne)),
      null_reason = if (is.na(pu)) attr(pu, "reason") else "", stringsAsFactors = FALSE)
    for (k in seq_along(seg)) {
      ev <- so$n_genes[[seg[k]]] > 0L && so$n_genes[[seg[k]]] >= p$min_genes; hs <- if (ev) inf(so$dP[seg[k], en$held]) else NULL
      sr[[length(sr) + 1L]] <- data.frame(program_id = p$program_id, segment_index = k, segment = seg[k], n_genes = so$n_genes[[seg[k]]],
        n_held_eligible = if (ev) hs$n_eligible else 0L, estimate = if (ev) hs$estimate else NA_real_, ci_low = if (ev) hs$ci_low else NA_real_,
        ci_high = if (ev) hs$ci_high else NA_real_, p_sign = if (ev) hs$sign_test$p_exact else NA_real_, description = .nrv_desc(hs, alpha),
        phase_assigned = identical(seg[k], so$sigma), value_definition = "signed_by_reference_direction", stringsAsFactors = FALSE)
    }
    if (!is.na(so$sigma) && length(p$genes)) {                                  # σ 의 유전자별 근거 (F3, BY 는 키 안에서)
      dtr <- colSums(matrix(arr$det[so$sigma, en$train, p$genes], length(en$train), length(p$genes)))
      for (j in which(dtr >= mind)) {
        x <- p$signs[j] * arr$delta[so$sigma, en$held, p$genes[j]]; tt <- mdu_sign_test(x)
        gr[[length(gr) + 1L]] <- data.frame(program_id = p$program_id, gene = p$genes[j], sigma = so$sigma, n_held = tt$n_nonmissing,
                                            mean_signed_delta = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE), p_sign = tt$p_exact, stringsAsFactors = FALSE)
      }
    }
  }
  list(programs = do.call(rbind, pr), segments = do.call(rbind, sr), genes = do.call(rbind, gr), nulls = nulls)
}

#' Holm (F1 = 성인 확인 층, F2 = 연령 미확인 층) → 판정 → phase 상태
.nrv_assign_verdicts <- function(prog, alpha) {
  prog$family <- ifelse(prog$layer == "adult_confirmed", "F1_adult_confirmed", "F2_human_age_unconfirmed")
  prog$p_holm <- NA_real_; prog$family_size <- 0L
  for (f in unique(prog$family)) {
    i <- which(prog$family == f & !is.na(prog$p_sign)); if (length(i)) prog$p_holm[i] <- stats::p.adjust(prog$p_sign[i], "holm")
    prog$family_size[prog$family == f] <- length(i)
  }
  prog$n_for_strictest_holm_threshold <- ifelse(prog$family_size > 0, ceiling(log2(2 * prog$family_size / alpha)), NA_real_)   # 설명 열 (거르지 않음)
  v <- vapply(seq_len(nrow(prog)), function(i) .nrv_verdict(prog$sigma[i], prog$p_holm[i], prog$ci_low[i], prog$ci_high[i], prog$p_null_up[i],
                                                            prog$p_null_down[i], prog$conditional[i], alpha), character(2))
  prog$verdict <- v[1, ]; prog$reason <- v[2, ]
  prog$phase_status <- ifelse(startsWith(prog$verdict, "supported"), ifelse(prog$single_gene, "phase_supported_single_gene", "phase_supported"),
                              ifelse(prog$verdict == "conflicting", "phase_reversed", "phase_insufficient"))
  prog
}

#' 내부 군집 하나의 집합 합 (block 순서 = blocks). 집합이 없는 block 은 0 (세포 0 → 기울기 NA)
.nrv_isum <- function(rs, triple, blocks, part) {
  m <- rs$int_meta; r <- which(m$triple == triple); b <- match(blocks, m$block[r]); ok <- !is.na(b)
  S1 <- rs$sums[[paste0("int1_", part)]]$sum; St <- rs$sums[[paste0("intt_", part)]]$sum
  s1 <- matrix(0, length(blocks), ncol(S1), dimnames = list(blocks, colnames(S1))); st <- s1
  if (any(ok)) { s1[ok, ] <- S1[r[b[ok]], , drop = FALSE]; st[ok, ] <- St[r[b[ok]], , drop = FALSE] }
  v <- function(x) { o <- setNames(numeric(length(blocks)), blocks); o[ok] <- x[r[b[ok]]]; o }
  list(n = v(m$n), t1 = v(m$t1), t2 = v(m$t2), S1 = s1, St = st)
}

#' 세 무리(A 쪽, 가운데, C 쪽) 프로그램 점수의 보류 block 평균 (p 없음, 서술)
.nrv_between <- function(rs, triple, held, genes, signs) {
  gm <- rs$grp_meta; S <- rs$sums$grp_ref$sum; w <- signs / length(genes)
  vapply(c("A_side", "middle", "C_side"), function(g) {
    r <- which(gm$triple == triple & gm$group == g & gm$block %in% held); if (!length(r)) return(NA_real_)
    mean(as.numeric(S[r, genes, drop = FALSE] %*% w) / gm$n[r])
  }, numeric(1))
}

#' 경로 내부 군집마다: 학습 block 기울기로 프로그램 선택 → 보류 block 검정, 같은 선택을 반복한 gene null, 보류 block 구조
.nrv_eval_internal <- function(path, P, en, sel, rs, N, a, minc, alpha) {
  m <- length(path); if (m < 3L) return(NULL)
  blocks <- c(en$train, en$held); Pi <- lapply(P, function(p) list(genes = p$genes, signs = p$signs, min_genes = p$min_genes)); rows <- list()
  for (k in 2:(m - 1L)) {
    id <- paste(path[k], path[k - 1L], path[k + 1L], sep = "|"); ic <- sel$icoord[[id]]
    io <- nrv_internal_program_stats(.nrv_isum(rs, id, blocks, "ref"), Pi, en$train, en$held, minc)
    h <- nrv_block_inference(io$heldout, R = a$replicates, seed = a$seed)
    ne <- rep(NA_real_, nrow(N))
    if (identical(ic$status, "span_ok")) {
      isp <- .nrv_isum(rs, id, blocks, "pool")
      ne <- vapply(seq_len(nrow(N)), function(r) nrv_internal_program_stats(isp, lapply(P, function(p) list(genes = paste0("col", N[r, p$null_rows]),
                     signs = p$signs, min_genes = p$min_genes)), en$train, en$held, minc)$estimate, numeric(1))
    }
    pu <- nrv_gene_null_p(io$estimate, ne); pd <- nrv_gene_null_p(-io$estimate, -ne)
    used <- names(io$heldout)[!is.na(io$heldout)]; bs <- ic$block_structure; span <- bs$block_structure[match(used, bs$block)]
    btw <- if (is.na(io$selected)) c(NA_real_, NA_real_, NA_real_) else .nrv_between(rs, id, en$held, P[[io$selected]]$genes, P[[io$selected]]$signs)
    rows[[length(rows) + 1L]] <- data.frame(position = k, B = path[k], A = path[k - 1L], C = path[k + 1L], pooled_structure = ic$status,
      selected_program = io$selected, selected_single_gene = if (is.na(io$selected)) NA else P[[io$selected]]$single_gene,
      train_stat = if (is.na(io$selected)) NA_real_ else io$train_stat[[io$selected]], n_held_eligible = h$n_eligible, estimate = h$estimate,
      ci_low = h$ci_low, ci_high = h$ci_high, p_sign = h$sign_test$p_exact, p_null_up = as.numeric(pu), p_null_down = as.numeric(pd),
      n_null_valid = sum(!is.na(ne)), n_held_blocks_used = length(used), n_held_span_ok = sum(span == "span_ok", na.rm = TRUE),
      block_structure_all_span_ok = length(used) > 0 && all(span == "span_ok", na.rm = FALSE) && !anyNA(span),
      between_A_side = btw[[1]], between_middle = btw[[2]], between_C_side = btw[[3]], geometry = "t_pooled", stringsAsFactors = FALSE)
  }
  it <- do.call(rbind, rows); i <- !is.na(it$p_sign); it$p_holm <- NA_real_
  if (any(i)) it$p_holm[i] <- stats::p.adjust(it$p_sign[i], "holm")                                  # F4: 경로의 내부 군집들
  v <- vapply(seq_len(nrow(it)), function(r) .nrv_internal_status(it$pooled_structure[r], it$selected_program[r], it$p_holm[r], it$ci_low[r],
                                                                  it$ci_high[r], it$p_null_up[r], it$p_null_down[r], alpha), character(2))
  it$status <- v[1, ]; it$reason <- v[2, ]; it
}

#' 파동 site: S 출구, k1 입구, k1 출구, …, T 입구 (그래프로만 정한 같은 block 경계). 보류 block x site 의 프로그램 점수 평균
.nrv_wave_Y <- function(ss, path, held, genes, signs, min_cells) {
  k <- length(path); sites <- character(0); sid <- character(0)
  for (s in seq_len(k - 1L)) {
    sites <- c(sites, sprintf("%02d_%s_exit", 2L * s - 1L, path[s]), sprintf("%02d_%s_entry", 2L * s, path[s + 1L]))
    sid <- c(sid, paste(path[s], path[s + 1L], sep = ">"), paste(path[s + 1L], path[s], sep = ">"))
  }
  w <- signs / length(genes); Y <- matrix(NA_real_, length(held), length(sites), dimnames = list(held, sites))
  for (j in seq_along(sites)) {
    r <- match(paste("sb", sid[j], held, sep = "|"), ss$meta$sb$set_id); nn <- ifelse(is.na(r), 0, ss$sums$sb$n[r]); ok <- nn >= min_cells
    if (any(ok)) Y[ok, j] <- as.numeric(ss$sums$sb$sum[r[ok], genes, drop = FALSE] %*% w) / nn[ok]
  }
  Y
}

.nrv_eval_waves <- function(ss, path, en, P, a, minc) {
  do.call(rbind, lapply(P, function(p) {
    if (!length(p$genes)) return(NULL)
    Y <- .nrv_wave_Y(ss, path, en$held, p$genes, p$signs, minc); wp <- nrv_wave_profile(Y, 1, R = a$replicates, seed = a$seed)
    data.frame(program_id = p$program_id, site_index = seq_len(ncol(Y)), site = colnames(Y), site_mean = unname(wp$site_mean),
               site_ci_low = unname(wp$site_ci_low), site_ci_high = unname(wp$site_ci_high), n_draws_site_missing = unname(wp$n_draws_site_missing),
               n_blocks_site = unname(wp$n_blocks_site), peak_freq = unname(wp$peak_freq), n_draws_no_peak = wp$n_draws_no_peak, R = wp$R,
               entry_site = wp$entry_site, exit_site = wp$exit_site, status = wp$status, stringsAsFactors = FALSE)
  }))
}

#' pooled 경계 대비 (block 을 가리지 않음): pooled_not_replicated 서술만
.nrv_eval_pooled <- function(ss, path, P, minc) {
  seg <- paste(path[-length(path)], path[-1], sep = ">"); m <- ss$meta$pool; S <- ss$sums$pool
  do.call(rbind, lapply(seq_along(seg), function(k) {
    rf <- match(paste0("pool|", path[k], ">", path[k + 1L]), m$set_id); rt <- match(paste0("pool|", path[k + 1L], ">", path[k]), m$set_id)
    nf <- if (is.na(rf)) 0 else S$n[rf]; nt <- if (is.na(rt)) 0 else S$n[rt]
    do.call(rbind, lapply(P, function(p) {
      ok <- nf >= minc && nt >= minc && length(p$genes) > 0
      d <- if (ok) sum(p$signs / length(p$genes) * (S$sum[rt, p$genes] / nt - S$sum[rf, p$genes] / nf)) else NA_real_
      data.frame(segment_index = k, segment = seg[k], program_id = p$program_id, pooled_signed_delta = d, n_from = nf, n_to = nt,
                 value_status = "pooled_not_replicated", stringsAsFactors = FALSE)
    }))
  }))
}

nrv_stage_evaluate <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  ct <- .nrv_rds(a, "freeze_controls", "controls_lossless.rds"); gs <- .nrv_rds(a, "stream_global", "gene_stats_lossless.rds")$gs
  ss <- .nrv_rds(a, "stream_sets", "set_sums.rds"); sel <- .nrv_rds(a, "select", "selection.rds"); rs <- .nrv_rds(a, "stream_route", "route_sums.rds")
  alpha <- .nrv_thr("alpha"); minc <- .nrv_thr("set_min_cells"); mind <- .nrv_thr("gene_min_detected_cells"); blocks <- cl$bm$block
  # 누설 검사 (배정마다): 판정 유전자 ∩ (anchor ∪ masked panel ∪ anchor·masked 대조 유전자) = ∅, null 후보 ∩ 제외 목록 = ∅
  leak <- do.call(rbind, lapply(fr$als, function(al) {
    tg <- ct$null[[al]]$targets; mids <- ct$modules$module_id[ct$modules$allocation == al & ct$modules$kind %in% c("anchor", "masked")]
    ctl_sym <- unique(gs$symbol[match(ct$controls$col[ct$controls$module_id %in% mids], gs$col)])
    selg <- union(fr$selg$gene[fr$selg$allocation == al & fr$selg$condition == "combined"], ctl_sym)
    nrv_check_no_leak(fr$roles[fr$roles$allocation == al, ], unique(tg$gene), selg)
    bad <- intersect(as.integer(ct$null[[al]]$N), ct$null[[al]]$exclude_cols)
    if (length(bad)) stop("gene null 후보가 제외 목록(참고·anchor 대조 유전자)과 겹친다: ", length(bad))
    data.frame(allocation = al, status = "ok", n_evaluation_genes = length(unique(tg$gene)), n_selection_genes = length(selg),
               n_null_pool = length(unique(as.integer(ct$null[[al]]$N))), n_null_excluded_overlap = 0L, stringsAsFactors = FALSE)
  }))
  keep <- list(program_results = list(), segment_descriptions = list(), gene_results = list(), internal_results = list(),
               wave_profiles = list(), pooled_segment_contrasts = list(), routes_not_evaluated = list())
  nulls <- list()
  for (id in names(sel$entries)) {
    en <- sel$entries[[id]]; kc <- .nrv_key_cols(a, en); N <- ct$null[[en$allocation]]$N; pool_cols <- paste0("col", sort(unique(as.integer(N))))
    for (kind in names(NRV_ROUTE_KINDS)) {
      rt <- en$routes[[kind]]; tag <- function(df) if (is.null(df) || !nrow(df)) NULL else cbind(kc, route_kind = kind, df, stringsAsFactors = FALSE)
      if (!identical(rt$status, "route_ok")) { keep$routes_not_evaluated[[length(keep$routes_not_evaluated) + 1L]] <- tag(data.frame(structural_path = rt$status)); next }
      P <- .nrv_route_programs(fr, ct, en$allocation, kind); path <- rt$primary
      arr <- nrv_segment_arrays(ss$sums$sb, ss$meta$sb, path, blocks, unique(unlist(lapply(P, `[[`, "genes"))), minc)
      arrp <- nrv_segment_arrays(rs$sums$route_pool, rs$route_meta, path, blocks, pool_cols, minc)
      ep <- .nrv_eval_programs(P, arr, arrp, N, en, a, mind, alpha)
      fz <- en$sigma[en$sigma$route_kind == kind, ]
      if (!identical(as.character(ep$programs$sigma), as.character(fz$sigma[match(ep$programs$program_id, fz$program_id)])))
        stop("σ 가 select 단계의 고정값과 다르다: ", id, " / ", kind)
      gr <- ep$genes; if (!is.null(gr)) { gr$p_by <- NA_real_; i <- !is.na(gr$p_sign); if (any(i)) gr$p_by[i] <- stats::p.adjust(gr$p_sign[i], "BY") }
      keep$program_results[[length(keep$program_results) + 1L]] <- tag(.nrv_assign_verdicts(ep$programs, alpha))
      keep$segment_descriptions[[length(keep$segment_descriptions) + 1L]] <- tag(ep$segments)
      keep$gene_results[[length(keep$gene_results) + 1L]] <- tag(gr)
      keep$internal_results[[length(keep$internal_results) + 1L]] <- tag(.nrv_eval_internal(path, P, en, sel, rs, N, a, minc, alpha))
      keep$wave_profiles[[length(keep$wave_profiles) + 1L]] <- tag(.nrv_eval_waves(ss, path, en, P, a, minc))
      keep$pooled_segment_contrasts[[length(keep$pooled_segment_contrasts) + 1L]] <- tag(.nrv_eval_pooled(ss, path, P, minc))
      nulls[[paste(id, kind, sep = "|")]] <- ep$nulls
    }
    log("evaluated:", id)
  }
  tabs <- lapply(keep, function(x) { x <- Filter(Negate(is.null), x); if (length(x)) do.call(rbind, x) else NULL })
  w$tsv(leak, "leak_check.tsv")
  for (nm in names(tabs)) if (!is.null(tabs[[nm]])) w$tsv(tabs[[nm]], paste0(nm, ".tsv"))
  w$rds(c(tabs, list(leak = leak, nulls = nulls)), "evaluation.rds")
}


# ── 9. verdict: 경로 다섯 근거와 라벨, 10 문헌 edge, 세 조건 표, 전환 집계 (저자 라벨 없음) ─────────────
.nrv_internal_summary <- function(st) {
  if (!length(st)) return("internal_not_applicable")
  if (all(st == "internally_continuous")) return("all_continuous_pooled_geometry")
  if (any(st == "insufficient_internal")) "some_insufficient" else "some_not_continuous"
}

nrv_stage_verdict <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); pr <- .nrv_rds(a, "structure", "pairs.rds")
  sel <- .nrv_rds(a, "select", "selection.rds"); ev <- .nrv_rds(a, "evaluate", "evaluation.rds"); ref <- fr$ref
  pick <- function(df, kc, kind) if (is.null(df)) NULL else
    df[df$holdout_scheme == kc$holdout_scheme & df$fold == kc$fold & df$allocation == kc$allocation & df$route_kind == kind, , drop = FALSE]
  none <- function(df) is.null(df) || !nrow(df)
  rs_rows <- list(); ea_rows <- list(); cc_rows <- list(); tr_rows <- list(); mh_rows <- list(); cm_rows <- list(); cmap <- list()
  for (id in names(sel$entries)) {
    en <- sel$entries[[id]]; kc <- .nrv_key_cols(a, en); al <- en$allocation; ro <- fr$roles[fr$roles$allocation == al, ]; rq <- fr$req[fr$req$allocation == al, ]
    key_rows <- list()
    for (kind in names(NRV_ROUTE_KINDS)) {
      rt <- en$routes[[kind]]; pg <- pick(ev$program_results, kc, kind); sg <- pick(ev$segment_descriptions, kc, kind); it <- pick(ev$internal_results, kc, kind)
      req_ids <- if (kind == "partial_path_early_granule") rq$program_id[rq$in_early_path] else rq$program_id
      phase <- setNames(rep("phase_insufficient", length(req_ids)), req_ids); pos <- setNames(rep(NA_real_, length(req_ids)), req_ids)
      if (!none(pg)) { phase[pg$program_id] <- pg$phase_status; pos[pg$program_id] <- pg$sigma_index }
      rank <- setNames(rq$lit_rank[match(req_ids, rq$program_id)], req_ids); lay <- rq$layer[match(req_ids, rq$program_id)]
      layers <- intersect(NRV_LAYERS, unique(lay))                                # 경로 종류에 필수 프로그램이 있는 층만
      cov <- setNames(lapply(layers, function(L) nrv_coverage(req_ids[lay == L], phase, pos, rank)), layers)
      ist <- if (none(it)) character(0) else setNames(it$status, it$B)
      seg_obs <- if (none(sg)) logical(0) else as.logical(tapply(sg$description != "insufficient", sg$segment_index, any))
      bs_all <- if (none(it)) NA else all(it$block_structure_all_span_ok %in% TRUE)
      T_ok <- any(en$anchors$anchor == NRV_ROUTE_KINDS[[kind]] & en$anchors$candidate)
      lab <- nrv_full_route_label(rt$status, T_ok, ist, seg_obs, phase, cov, bs_all)
      cs <- function(L) if (L %in% names(cov)) cov[[L]]$status else "not_applicable"
      key_rows[[kind]] <- cbind(kc, route_kind = kind, .nrv_route_row(rt), anchor_identity = "candidate_only",
        S_candidates = paste(en$anchors$cluster[en$anchors$anchor == "S" & en$anchors$candidate], collapse = "|"),
        T_candidates = paste(en$anchors$cluster[en$anchors$anchor == NRV_ROUTE_KINDS[[kind]] & en$anchors$candidate], collapse = "|"),
        internal_continuity = .nrv_internal_summary(ist), n_interior = length(ist), n_internally_continuous = sum(ist == "internally_continuous"),
        block_structure_all_span_ok = bs_all, n_segments = length(seg_obs), n_segments_observed = sum(seg_obs),
        coverage_adult_confirmed = cs("adult_confirmed"), coverage_human_age_unconfirmed = cs("human_age_unconfirmed"),
        full_route_label = lab$label, reason_codes = paste(lab$reason_codes, collapse = ";"), single_gene_programs = paste(lab$single_gene_programs, collapse = "|"),
        n_programs_supported = if (none(pg)) 0L else sum(startsWith(pg$verdict, "supported")), n_programs_conflicting = if (none(pg)) 0L else sum(pg$verdict == "conflicting"),
        partial_supported_segments = if (none(pg)) "" else paste(unique(pg$sigma[startsWith(pg$verdict, "supported")]), collapse = "|"),
        layer_note = "층별로 따로 보고; 두 층을 모두 만족해도 성인 문헌만으로 검증된 전체 경로라고 부르지 않는다", stringsAsFactors = FALSE)
    }
    rsk <- do.call(rbind, key_rows); nrv_assert_single_key(rsk); rs_rows[[id]] <- rsk
    # 10 문헌 edge: 경로 근거(anchor_route)와 결합(combined). 둘 다 full route 의 프로그램 결과를 쓴다.
    pgf <- pick(ev$program_results, kc, "full_route"); itf <- pick(ev$internal_results, kc, "full_route")
    prf <- if (none(pgf)) data.frame(program_id = character(0), status = character(0), stringsAsFactors = FALSE) else
      data.frame(program_id = pgf$program_id, status = pgf$verdict, stringsAsFactors = FALSE)
    ea_route <- nrv_edge_assessment(ref, ro, prf, mapping = NULL)
    rtf <- en$routes$full_route; route_full <- if (identical(rtf$status, "route_ok")) rtf$primary else character(0)
    cand_m <- setNames(strsplit(en$marker_masked$candidate_set, "|", fixed = TRUE), en$marker_masked$cluster)
    cand_f <- setNames(strsplit(en$marker_full$candidate_set, "|", fixed = TRUE), en$marker_full$cluster)
    ph <- if (none(pgf)) data.frame(edge_id = character(0), program_id = character(0), segment_index = integer(0), status = character(0),
                                    conditional_only = logical(0), stringsAsFactors = FALSE) else
      data.frame(edge_id = pgf$edge_id, program_id = pgf$program_id, segment_index = pgf$sigma_index, status = pgf$phase_status,
                 conditional_only = pgf$conditional, stringsAsFactors = FALSE)
    cm <- nrv_combined_mapping(route_full, cand_m, ref$edges, ph, if (none(itf)) NULL else setNames(itf$status, itf$B))
    ea_comb <- nrv_edge_assessment(ref, ro, prf, mapping = cm); cmap[[id]] <- cm
    mh <- nrv_marker_edge_hypotheses(pr$folds[[en$key]], cand_m, ref$edges); mhf <- nrv_marker_edge_hypotheses(pr$folds[[en$key]], cand_f, ref$edges)
    ea_rows[[id]] <- rbind(cbind(kc, condition = "anchor_route", ea_route, stringsAsFactors = FALSE), cbind(kc, condition = "combined", ea_comb, stringsAsFactors = FALSE))
    cc_rows[[id]] <- nrv_condition_table(as.list(kc), marker_only = data.frame(edge_id = mh$edge_id, status = mh$status),
                                         route = data.frame(edge_id = ea_route$edge_id, status = ea_route$verdict),
                                         combined = data.frame(edge_id = ea_comb$edge_id, status = ea_comb$verdict))
    tr_rows[[id]] <- cbind(kc, nrv_transitions(data.frame(edge_id = mh$edge_id, status = mh$status), data.frame(edge_id = ea_comb$edge_id, status = ea_comb$verdict)),
                           stringsAsFactors = FALSE)
    mh_rows[[id]] <- rbind(cbind(kc, panel = "masked", mh, stringsAsFactors = FALSE), cbind(kc, panel = "full_circular", mhf, stringsAsFactors = FALSE))
    cm_rows[[id]] <- cbind(kc, route = paste(route_full, collapse = ">"), cm, stringsAsFactors = FALSE)
  }
  bind <- function(x) do.call(rbind, unname(x))
  rs_all <- bind(rs_rows); ea_all <- bind(ea_rows)
  log("route labels:", paste(names(table(rs_all$full_route_label)), table(rs_all$full_route_label), sep = "=", collapse = ", "))
  w$tsv(ea_all, "reference_edge_assessment.tsv"); w$tsv(rs_all, "route_summary.tsv"); w$tsv(bind(cc_rows), "condition_comparison.tsv")
  w$tsv(bind(tr_rows), "transitions.tsv"); w$tsv(bind(mh_rows), "marker_hypotheses.tsv"); w$tsv(bind(cm_rows), "combined_mapping.tsv")
  w$rds(list(route_summary = rs_all, reference_edge_assessment = ea_all, condition_comparison = bind(cc_rows), transitions = bind(tr_rows),
             marker_hypotheses = bind(mh_rows), combined_mapping = cmap), "verdict_lossless.rds")
}


# ── 10. author: 저자 라벨을 처음 읽는다. 세 조건을 같은 세포·군집 분모로 센다 ─────────────────────────
#' node → 저자 라벨: 기존 evaluation_mapping_used.tsv 에서 이름이 같은 node 만 (해결된 것만). 나머지는 unmapped_reference.
.nrv_node_author_map <- function(em, nodes) {
  res <- if ("resolved" %in% names(em)) em$resolved == "yes" else rep(TRUE, nrow(em))
  em <- em[res & em$state %in% nodes & !is.na(em$exact_fine) & nzchar(em$exact_fine), , drop = FALSE]
  setNames(strsplit(em$exact_fine, "|", fixed = TRUE), em$state)
}

#' 한 군집의 조건 라벨(node 후보)과 저자 라벨 비교: 후보가 없거나 여럿·unassigned·unavailable 이면 held, 대응 없는 node 는 unmapped
.nrv_compare <- function(nodes, labels, nm) {
  if (!length(nodes) || length(nodes) > 1L || nodes %in% c("unassigned", "unavailable")) return(rep("held", length(labels)))
  if (!nodes %in% names(nm)) return(rep("unmapped", length(labels)))
  ifelse(labels %in% nm[[nodes]], "agree", "disagree")
}

#' combined 조건의 군집 라벨: 대응이 single 이고 결합 근거가 supported 인 문헌 edge 의 부분경로 양끝에만 u, v 를 붙인다. 그 밖은 held.
.nrv_combined_labels <- function(cm, ref_edges) {
  out <- list(); if (is.null(cm) || !nrow(cm)) return(out)
  for (r in which(cm$mapping_status == "single" & startsWith(cm$combined_status, "supported"))) {
    sr <- strsplit(cm$subroute[r], ">", fixed = TRUE)[[1]]; ed <- ref_edges[ref_edges$edge_id == cm$edge_id[r], ]
    out[[sr[1]]] <- unique(c(out[[sr[1]]], ed$from)); out[[sr[length(sr)]]] <- unique(c(out[[sr[length(sr)]]], ed$to))
  }
  out
}

nrv_stage_author <- function(a, reader, log, w, ...) {
  .nrv_frozen(a); fr <- .nrv_rds(a, "freeze", "freeze_lossless.rds"); cl <- .nrv_rds(a, "structure", "cells_lossless.rds")
  sel <- .nrv_rds(a, "select", "selection.rds"); vd <- .nrv_rds(a, "verdict", "verdict_lossless.rds")
  au <- as.character(nrv_stage_reads_author("author", reader))                    # 저자 라벨은 여기서 처음 읽는다
  if (length(au) != nrow(cl$cells) || anyNA(au)) stop("저자 라벨 수가 세포 수와 다르거나 결측이 있다.")
  nm <- .nrv_node_author_map(reader$eval_mapping(), fr$ref$nodes$node_id)
  clu <- cl$cells$cluster; ucl <- .nrv_lv(clu); tab <- table(factor(clu, levels = ucl), au)
  maj <- setNames(colnames(tab)[apply(tab, 1, which.max)], ucl)                  # 군집 다수 저자 라벨 (동률이면 표 순서상 앞, 서술용)
  lev <- c("agree", "disagree", "unmapped", "held"); rows <- list(); comp <- list()
  for (id in names(sel$entries)) {
    en <- sel$entries[[id]]; kc <- .nrv_key_cols(a, en)
    labs <- list(marker_only = setNames(strsplit(en$marker_masked$candidate_set, "|", fixed = TRUE), en$marker_masked$cluster),
                 combined = .nrv_combined_labels(vd$combined_mapping[[id]], fr$ref$edges))
    for (cond in c("marker_only", "anchor_route", "combined")) {
      if (cond == "anchor_route") {                                                # 상태 라벨이 없다 → 비교하지 않는다
        rows[[length(rows) + 1L]] <- cbind(kc, condition = cond, unit = c("cluster", "cell"), comparison = "non_comparable", n = c(length(ucl), length(au)),
                                           stringsAsFactors = FALSE)
        next
      }
      L <- labs[[cond]]; ccl <- vapply(ucl, function(k) .nrv_compare(L[[k]], maj[[k]], nm), character(1)); ccell <- character(length(au))
      for (k in ucl) { i <- which(clu == k); ccell[i] <- .nrv_compare(L[[k]], au[i], nm) }
      cnt <- function(x, unit) cbind(kc, condition = cond, unit = unit, comparison = lev, n = as.integer(table(factor(x, levels = lev))), stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- rbind(cnt(ccl, "cluster"), cnt(ccell, "cell"))
    }
    for (kind in names(en$routes)) {                                               # 경로 군집의 저자 라벨 구성 (서술)
      rt <- en$routes[[kind]]; if (!identical(rt$status, "route_ok")) next
      for (k in seq_along(rt$primary)) {
        v <- sort(tab[rt$primary[k], ], decreasing = TRUE); v <- v[v > 0]; v <- v[seq_len(min(5L, length(v)))]
        comp[[length(comp) + 1L]] <- cbind(kc, route_kind = kind, position = k, cluster = rt$primary[k], author_label = names(v), n = as.integer(v),
                                           frac = as.numeric(v) / sum(tab[rt$primary[k], ]), stringsAsFactors = FALSE)
      }
    }
  }
  w$tsv(do.call(rbind, rows), "author_condition_comparison.tsv")
  if (length(comp)) w$tsv(do.call(rbind, comp), "author_route_composition.tsv")
  lc <- as.data.frame(tab, stringsAsFactors = FALSE); names(lc) <- c("cluster", "author_label", "n"); lc <- lc[lc$n > 0, ]
  lc$majority <- lc$author_label == maj[lc$cluster]; w$tsv(lc, "author_cluster_label_counts.tsv.gz")
  w$tsv(data.frame(node_id = names(nm), author_labels = vapply(nm, paste, "", collapse = "|"), stringsAsFactors = FALSE), "node_author_mapping_used.tsv")
  if (identical(a$partition, "p51") && is.function(reader$readout51)) {         # 51 의 기존 marker-only 판정 (92 의 baseline 으로 복사하지 않는다)
    ro <- reader$readout51(); ro <- ro[ro$fit == "descriptive" & ro$readout == "M_new_simple", c("unit", "call", "status")]
    names(ro) <- c("cluster", "existing_call_51", "existing_status_51"); ro$cluster <- as.character(ro$cluster)
    w$tsv(merge(data.frame(cluster = ucl, author_majority = unname(maj[ucl]), stringsAsFactors = FALSE), ro, by = "cluster", all.x = TRUE), "existing_call_51.tsv")
  }
  log("author comparison written for", length(sel$entries), "entries")
}


# ── 11. view: 28 을 스냅샷 사본에서 적재해 그림, Cytoscape 표, 그림 출처를 쓴다 ──────────────────────────
nrv_stage_view <- function(a, reader, log, w, snap, ...) {
  .nrv_frozen(a); k <- which(basename(snap$snapshot_path) == "28_view_nrv.R")
  if (length(k) != 1L || !identical(.nrv_sha(snap$snapshot_path[k]), snap$sha256[k])) stop("28 스냅샷 사본 hash 가 고정값과 다르다.")
  env <- new.env(parent = environment(nrv_stage_view)); sys.source(snap$snapshot_path[k], envir = env)
  log("loaded from snapshot: 28_view_nrv.R")
  env$nrv_view(a = a, reader = reader, log = log, w = w, snap = snap)
}


# 스크립트로 실행할 때: bootstrap 에 필요한 R18(스냅샷)·R23(배정 이름)만 먼저 읽는다. 단계 계산은 스냅샷 사본에서 한다.
if (!interactive() && sys.nframe() == 0L) {
  .nrv_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("18_run_marker_donor_uncertainty.R", "23_nrv_reference.R")) source(file.path(.nrv_bb, .f))
  quit(status = nrv_main(commandArgs(trailingOnly = TRUE)))
}
