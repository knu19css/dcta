# =============================================================================
# 60_run_nso.R — 호중구 성숙 단계 관측(NSO) 실행기: 입력 고정 → 원 RNA 근거 → 관측 연결 → 계산 표 고정 → (그 뒤) 원 UMAP 표시 · 사후 원저자 표
# =============================================================================
#
# 계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md (2026-10-19), 계획 Task 3. Ruling: 같은 폴더 claude_progress.md
#   - 입력 hash: 참고표 2 개 · graph_base.npz · h5ad 는 codex_input_baseline / 계획 값과 같아야 하고, graph · graph_ids · h5ad 는 O 의 NRV structure 기록과도 같아야 한다.
#   - O 의 NRV p92 · p51 stage_structure (cells_lossless · edges · pairs) 를 stage manifest hash 로 확인해 재사용한다. 군집 쌍 표는 다시 계산해 저장본과 대조한다.
#   - 계산 표를 모두 쓰고 calc_manifest.json 으로 hash 를 고정한 뒤에야 원 UMAP(표시)과 원저자 라벨(posthoc_ 표)을 읽는다.
#   - 예측 stage_label 은 만들지 않는다. 결과 폴더가 비어 있지 않으면 쓰지 않는다(exit 6).
#
# 실행 (N 에서, dcta 환경만)
#   run-r dcta R/branch_biology/60_run_nso.R --data-root <M> --o-root <O> --reference-dir <R> \
#     --out-dir <N>/results/neut_stage_observation_2026-10-19/attempt_01 [--cma-dir <R57 attempt_02>] [--fdt-dir <R54 attempt_01>] [--chunk 5000]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

NSO_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R")
NSO_GRAPH <- "results/projection_cluster_context_2026-09-25/attempt_01/BoneMarrowMap/graph_base.npz"       # M 기준 (R27 NRV_INPUTS 와 같은 경로)
NSO_GRAPH_IDS <- "results/projection_cluster_context_2026-09-25/attempt_01/BoneMarrowMap/graph_ids.txt"
NSO_H5AD <- "data/external/hay_hca_source.h5ad"
NSO_NRV <- c(p92 = "results/neutrophil_route_validation_2026-10-13/attempt_01_p92", p51 = "results/neutrophil_route_validation_2026-10-13/attempt_01_p51")   # O 기준
NSO_CMA_DIR <- file.path("results", "cluster_marker_annotation_2026-10-16", "attempt_02")    # N 기준: R57 군집 주석 (참조 열)
NSO_FDT_DIR <- file.path("results", "full_dag_top1_2026-10-16", "attempt_01")                # N 기준: R54 군집 x 39 패널 M (참조 열)
NSO_COMPETING <- c("GMP", "Neutrophil", "Eosinophil", "Basophil", "Monocyte", "CD14_Monocyte", "CD16_Monocyte")   # Ruling 6: 결과를 보기 전에 고정한 경쟁 후보 참조 열
NSO_PINNED <- list(reference = c(nodes = "8059caa90f7bf67561f2fd89603b493032435a295a34c2e327209adec83368f8",
                                 edges = "73f9f61fe35cdac0c9162bca930c6294737712c5e2973685f76c8b44eeadb58b"),
                   graph = "e63392de28aee7e4b629221bc8a0ba6118e5048c27b55071019b1d247ffa62de", h5ad = "4149b357b7d0027027988b8c78a92ae8894f36f3eae667762500dbc976f3c454")


# ── 1. 표 조립 (순수) ────────────────────────────────────────────────────────
#' 다시 계산한 군집 쌍 표와 NRV 저장본이 같은지 (키 순서 그대로, 수치는 1e-12 안)
.nso_same_pairs <- function(x, y) {
  if (nrow(x) != nrow(y) || !setequal(names(x), names(y))) return(FALSE)
  if (!nrow(x)) return(TRUE)
  num <- setdiff(names(x), c("cluster_a", "cluster_b"))
  identical(as.character(x$cluster_a), as.character(y$cluster_a)) && identical(as.character(x$cluster_b), as.character(y$cluster_b)) &&
    all(vapply(num, function(k) isTRUE(all.equal(as.numeric(x[[k]]), as.numeric(y[[k]]), tolerance = 1e-12)), TRUE))
}

#' 참조 열: R57 군집 주석(원저자 열은 읽지 않는다)과 R54 M 의 고정 경쟁 타입. 둘 다 그 결과 manifest hash 로 확인한다.
nso_reference_comparison <- function(cma_dir, fdt_dir, md) {
  fl <- .mbvq_manifest_file(fdt_dir, "manifest.json", "lossless.rds"); sc <- readRDS(fl)$scores          # 점수만 쓴다 (같은 파일의 사후 저자 구성표는 쓰지 않는다)
  tabs <- list(); inputs <- c(r54_lossless = fl)
  for (p in c("p92", "p51")) { cl <- if (p == "p92") md$cluster92 else md$cluster51; uc <- fdt_cluster_order_safe(cl); nc <- table(factor(cl, levels = uc))
    f <- .mbvq_manifest_file(cma_dir, "manifest.json", file.path(p, "cluster_marker_annotation.tsv")); inputs[[paste0("r57_cluster_", p)]] <- f
    hdr <- strsplit(readLines(f, n = 1L), "\t", fixed = TRUE)[[1]]
    t57 <- utils::read.delim(f, colClasses = stats::setNames(ifelse(grepl("author", hdr, ignore.case = TRUE), "NULL", "character"), hdr), quote = "", check.names = FALSE)
    if (!all(c("cluster", "predicted_type", "status", "warnings") %in% names(t57))) stop("R57 군집 주석에 필요한 열이 없다: ", p)
    if (anyDuplicated(t57$cluster) || !setequal(t57$cluster, uc)) stop("R57 군집 주석의 군집이 세포표와 다르다: ", p)
    M <- sc[[p]]$M; nn <- sc[[p]]$n
    if (is.null(M) || !setequal(names(nn), uc) || any(as.numeric(nn[uc]) != as.numeric(nc))) stop("R54 군집 크기가 세포표와 다르다: ", p)
    miss <- setdiff(NSO_COMPETING, colnames(M)); if (length(miss)) stop("R54 M 에 참조 타입이 없다: ", paste(miss, collapse = ", "))
    j <- match(uc, t57$cluster); Mx <- M[uc, NSO_COMPETING, drop = FALSE]; colnames(Mx) <- paste0("r54_M_", NSO_COMPETING)
    tabs[[p]] <- data.frame(partition = p, cluster = uc, n_cells = as.integer(nc), r57_predicted_type = t57$predicted_type[j], r57_status = t57$status[j],
                            r57_warnings = t57$warnings[j], Mx, row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE) }
  list(tables = tabs, inputs = inputs)
}

#' 단계 · 변화 참고 유전자의 입력 존재 · 전체 세포 검출 (관측 가능성을 확정한 단계표가 아니다)
nso_availability <- function(ref, ge) {
  n <- ge$n_cells; have <- colnames(ge$cell$log)
  nd <- stats::setNames(as.numeric(Matrix::colSums(ge$cell$det)), have); ml <- stats::setNames(as.numeric(Matrix::colMeans(ge$cell$log)), have)
  mk <- ref$stage_markers; ch <- ref$changes
  a <- data.frame(item_kind = "stage_marker", stage_or_edge = mk$stage_id, from = "", to = "", gene = mk$gene, direction = "", evidence = "", shared_with = mk$shared_with,
                  pb_shared = mk$pb_shared, evidence_scope = mk$evidence_scope, stringsAsFactors = FALSE)
  b <- if (nrow(ch)) data.frame(item_kind = "transition_change", stage_or_edge = ch$edge_id, from = ch$from, to = ch$to, gene = ch$gene, direction = ch$direction,
                                evidence = ch$evidence, shared_with = "", pb_shared = FALSE, evidence_scope = ch$reference_scope, stringsAsFactors = FALSE) else NULL
  d <- rbind(a, b); ii <- d$gene %in% have
  d$in_input <- ii; d$raw_cols <- ge$genes$raw_cols[match(d$gene, ge$genes$gene)]
  d$n_cells_detected <- ifelse(ii, unname(nd[d$gene]), NA_real_); d$det_rate_all <- d$n_cells_detected / n; d$mean_log_all <- ifelse(ii, unname(ml[d$gene]), NA_real_)
  d$availability <- ifelse(!ii, "missing_gene", ifelse(d$n_cells_detected == 0, "observed_zero", "observed"))
  d$interpretation <- "gene presence and detection over all cells of this dataset; not a stage call and not proof that a stage is observable or absent"
  d
}

#' 세포 표: 전체 원 cell_id 순서, 단계 marker 6 개 근거(검출 유전자 수 · logRNA 평균), 참고 유전자별 logRNA, unavailable 이유, 기존 군집 · R57 주석 참조. 예측 단계 없음.
nso_cell_table <- function(md, ge, pe, cmp) {
  ids <- paste0("marker_", NSO_STAGES); cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  d <- data.frame(cell_id = md$cell_id, cluster92 = md$cluster92, cluster51 = md$cluster51, block = md$block, study = md$study, library_size = ge$cell$library,
                  stringsAsFactors = FALSE)
  for (k in ids) { d[[paste0(k, "_n_detected")]] <- cn[, k]; d[[paste0(k, "_mean_log")]] <- cm[, k] }
  L <- ge$cell$log; for (g in colnames(L)) d[[paste0("log_", g)]] <- as.numeric(L[, g])
  d$unavailable_reason <- pe$cells$unavailable
  t92 <- cmp$p92; t51 <- cmp$p51; j92 <- match(md$cluster92, t92$cluster)
  d$r57_type_p92 <- t92$r57_predicted_type[j92]; d$r57_status_p92 <- t92$r57_status[j92]; d$r57_type_p51 <- t51$r57_predicted_type[match(md$cluster51, t51$cluster)]
  d
}

#' spot-check 표: p92 군집마다 원 행 순서 첫 세포(라벨 · 값을 보지 않는 고정 규칙)의 원 count · library · logRNA
nso_spot_table <- function(md, ge) {
  k <- ge$kept; if (is.null(k$counts) || !length(k$cell)) return(data.frame(cell_id = character(0), stringsAsFactors = FALSE))
  d <- data.frame(cell_id = md$cell_id[k$cell], cell_index = k$cell, cluster92 = md$cluster92[k$cell], block = md$block[k$cell], library_size = k$library, stringsAsFactors = FALSE)
  for (g in colnames(k$counts)) { d[[paste0("count_", g)]] <- as.numeric(k$counts[, g]); d[[paste0("log_", g)]] <- as.numeric(ge$cell$log[k$cell, g]) }
  d$rule <- "first cell in raw order of each p92 cluster; no labels or values used to pick; log = log1p(1e4 * count / library_size)"
  d
}

#' 한계 표: reference 구조 · 자료 가용성 · 연결 해석의 한계를 전체 분모와 함께 (원저자 라벨을 읽지 않고 만든다)
nso_limitations <- function(ref, ge, pe, oc, icn) {
  rows <- list(); add <- function(item, scope, evidence, consequence) rows[[length(rows) + 1L]] <<- data.frame(item = item, scope = scope, evidence = evidence,
                                                                                                                consequence = consequence, stringsAsFactors = FALSE)
  sm <- ref$shared_markers
  for (i in which(sm$n_stages >= 2L)) add("shared_marker", sm$stages[i], paste0(sm$gene[i], " is listed for ", gsub("|", " and ", sm$stages[i], fixed = TRUE)),
                                          "this gene alone cannot separate these stages; their marker-program values overlap by construction")
  for (i in which(sm$pb_shared)) add("pb_shared_marker", paste0(sm$stages[i], "|PB_PMN"), paste0(sm$gene[i], " is also listed for peripheral-blood PMN"),
                                     "late bone-marrow and blood-neutrophil programs share this gene; marrow_to_blood is kept separate")
  d <- pe$programs
  for (i in which(d$single_gene_program)) add("single_gene_program", d$program_id[i], paste0("only ", d$genes[i]), "single-gene evidence, reported as such")
  for (i in which(d$evidence == "conditional")) add("conditional_change", d$program_id[i], d$genes[i], "auxiliary candidate trend from the reference; not a required gate")
  for (i in which(d$status != "ok")) add("program_unavailable", d$program_id[i], paste0(d$status[i], ": ", d$missing_genes[i]), "no value is reported (not replaced by 0)")
  g <- ge$genes; for (i in which(!g$in_input)) add("missing_gene", g$gene[i], "symbol not found in raw/var/gene_symbols", "gene evidence unavailable; not biological absence")
  nd <- Matrix::colSums(ge$cell$det); for (k in names(nd)[nd == 0]) add("observed_zero_gene", k, "present in input, zero counts in every cell", "observed zero is not proof that a stage is absent")
  tr <- ref$transitions; for (i in which(tr$status != "change_rules_listed")) add("transition_without_rules", paste(tr$from[i], tr$to[i], sep = ">"), tr$status[i],
                                                                                 "no reference change rule for this adjacent pair")
  for (i in which(ref$scope$role == "marrow_to_blood")) add("marrow_to_blood_out_of_scope", ref$scope$id[i], paste(ref$scope$from[i], ref$scope$to[i], sep = ">"),
                                                            "tissue move from bone marrow to blood, not a marrow morphological stage; separate reference table, no data computed")
  add("evidence_scope", paste(unique(ref$changes$edge_id), collapse = "|"), paste(unique(ref$changes$reference_scope), collapse = "|"),
      "human marrow evidence with donor age unconfirmed in the review; not adult-validated adjacent-stage DEGs")
  add("evaluation_reference", "author labels", "spec 2026-10-18 section 14: the 55 author labels contain no exact Myeloblast/PM/MY/MM/BN/SN label",
      "stage accuracy cannot be computed here; author labels were not used in any calculation (post hoc tables only)")
  for (p in names(icn)) { tb <- table(factor(icn[[p]]$triples$status, levels = c("span_ok", "internal_disconnected", "no_internal_span", "no_boundary")))
    add("internal_continuity", p, paste(names(tb), as.integer(tb), sep = "=", collapse = ", "),
        "graph adjacency of A-B and B-C does not imply that the two boundaries of B are connected inside B") }
  for (p in names(oc)) { f <- oc[[p]]$features
    tb <- table(factor(f$status[f$feature_kind == "gene"], levels = c("ok", "single_block", "no_eligible_block", "no_same_block_boundary", "missing_gene")))
    add("same_block_contrast_availability", p, paste(names(tb), as.integer(tb), sep = "=", collapse = ", "),
        "equal-block contrast exists only where both boundary sides have enough cells in the same block; the pooled contrast is descriptive") }
  add("orientation", "observed_edge_evidence", NSO_ORIENTATION, "no time, velocity or differentiation direction is inferred")
  add("fixed_graph", "graph and partitions", "graph, latent space and partitions were built on all cells", "block-level contrasts are conditional descriptions, not external validation")
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

#' 사후 (고정 뒤): 군집 x 원저자 라벨 세포 수 · 군집 안 비율 (정확도 아님)
nso_posthoc_composition <- function(md, au) {
  lab <- ifelse(is.na(au) | !nzchar(au), "(no author label)", as.character(au))
  out <- do.call(rbind, lapply(c("p92", "p51"), function(p) { cl <- if (p == "p92") md$cluster92 else md$cluster51
    tb <- as.data.frame(table(cluster = cl, author_label = lab), stringsAsFactors = FALSE); tb <- tb[tb$Freq > 0, , drop = FALSE]; names(tb)[3] <- "n_cells"
    tot <- tapply(tb$n_cells, tb$cluster, sum); tb$frac_of_cluster <- tb$n_cells / as.numeric(tot[tb$cluster])
    tb <- tb[order(match(tb$cluster, fdt_cluster_order_safe(cl)), -tb$n_cells, tb$author_label, method = "radix"), , drop = FALSE]
    data.frame(partition = p, tb, row.names = NULL, stringsAsFactors = FALSE) }))
  out$note <- "post hoc after results were fixed; not used for any selection; not accuracy"
  out
}

#' 사후 (고정 뒤): 원저자 라벨마다 단계 marker 프로그램의 세포 평균과 검출 세포 비율 (서술)
nso_posthoc_programs <- function(au, pe) {
  lab <- ifelse(is.na(au) | !nzchar(au), "(no author label)", as.character(au)); ul <- sort(unique(lab), method = "radix"); fl <- factor(lab, levels = ul)
  cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  d <- data.frame(author_label = ul, n_cells = as.integer(table(fl)), stringsAsFactors = FALSE)
  for (k in paste0("marker_", NSO_STAGES)) {
    d[[paste0(k, "_mean_log")]] <- as.numeric(tapply(cm[, k], fl, function(v) if (all(is.na(v))) NA_real_ else mean(v)))
    d[[paste0(k, "_frac_detected")]] <- as.numeric(tapply(cn[, k] > 0, fl, function(v) if (all(is.na(v))) NA_real_ else mean(v))) }
  d$note <- "post hoc description after results were fixed; not used for any selection; not accuracy"
  d
}


# ── 2. 실행기 ────────────────────────────────────────────────────────────────
#' 인자 해석 → 결과 폴더 검사(6) → 소스 스냅샷 → 사본 적재 → 사본의 nso_build. root / reader / pins 는 테스트용.
nso_main <- function(argv, root = NULL, reader = NULL, pins = NSO_PINNED) {
  a <- list(data_root = NULL, o_root = NULL, ref_dir = NULL, out = NULL, cma = NULL, fdt = NULL, chunk = "5000")
  keys <- c("--data-root" = "data_root", "--o-root" = "o_root", "--reference-dir" = "ref_dir", "--out-dir" = "out", "--cma-dir" = "cma", "--fdt-dir" = "fdt", "--chunk" = "chunk")
  i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  a$chunk <- suppressWarnings(as.integer(a$chunk))
  bad <- vapply(a[c("data_root", "o_root", "ref_dir", "out")], function(x) is.null(x) || is.na(x) || !nzchar(x), TRUE)
  if (any(bad) || is.na(a$chunk) || a$chunk < 1L) { cat("--data-root, --o-root, --reference-dir, --out-dir 와 양의 정수 --chunk 가 필요하다.\n"); return(2L) }
  out <- a$out
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  if (is.null(a$cma)) a$cma <- file.path(root, NSO_CMA_DIR)
  if (is.null(a$fdt)) a$fdt <- file.path(root, NSO_FDT_DIR)
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", NSO_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(NSO_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = NSO_SOURCES)
    get("nso_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, reader = reader, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

nso_build <- function(a, out, snap, log, argv, reader, pins) {
  ck <- .mbvq_checks(); if (is.null(reader)) reader <- nrv_h5_reader(a$data_root)
  # ── 입력 고정: 계약 hash · NRV structure 기록 ──
  rf <- c(reference_nodes = file.path(a$ref_dir, "02_nodes.csv"), reference_edges = file.path(a$ref_dir, "03_edges.csv"))
  rs <- vapply(rf, .nrv_sha, "")
  if (!is.null(pins) && !identical(unname(rs), unname(pins$reference[c("nodes", "edges")]))) stop("참고표 hash 가 고정값(codex_input_baseline)과 다르다: ", paste(basename(rf), rs, collapse = ", "))
  st <- lapply(NSO_NRV, function(p) { att <- file.path(a$o_root, p)
    c(cells = nrv_stage_file(att, "structure", "cells_lossless.rds"), edges = nrv_stage_file(att, "structure", "edges.rds"),
      pairs = nrv_stage_file(att, "structure", "pairs.rds"), input = nrv_stage_file(att, "structure", "input_sha256.tsv")) })
  fk <- c("graph", "graph_ids", "h5ad"); f <- reader$input_files()[fk]
  if (anyNA(f) || !all(file.exists(f))) stop("입력 파일이 없다: ", paste(fk[is.na(f) | !file.exists(f)], collapse = ", "))
  cur <- vapply(f, .nrv_sha, "")
  for (p in names(st)) { rec <- utils::read.delim(st[[p]][["input"]], stringsAsFactors = FALSE, colClasses = "character")
    r <- rec$sha256[match(fk, rec$input)]
    if (anyNA(r) || !identical(unname(cur), r)) stop("현재 입력 hash 가 NRV structure 단계 기록과 다르다 (", p, "): ", paste(fk[is.na(r) | r != cur], collapse = ", ")) }
  if (!is.null(pins) && (!identical(unname(cur[["graph"]]), pins$graph) || !identical(unname(cur[["h5ad"]]), pins$h5ad))) stop("graph_base.npz 또는 h5ad hash 가 계약값과 다르다.")
  ck$add("inputs_pinned", TRUE, sprintf("reference nodes/edges, graph_base.npz and h5ad sha256 == codex_input_baseline and plan values%s; graph, graph_ids, h5ad == NRV structure input_sha256 (p92, p51)",
                                        if (is.null(pins)) " (pins disabled)" else ""))
  log("inputs pinned")
  # ── reference · 세포 · 계산용 메타데이터 ──
  ref <- nso_reference(rf[["reference_nodes"]], rf[["reference_edges"]])
  c92 <- readRDS(st$p92[["cells"]])$cells; c51 <- readRDS(st$p51[["cells"]])$cells
  chk <- nso_check_inputs(c92, c51, raw_ids = reader$obs("cell_id"), graph_ids = reader$ids(), reference = ref)
  for (r in seq_len(nrow(chk$checks))) ck$add(chk$checks$check[r], chk$checks$pass[r], chk$checks$detail[r])
  md <- nso_calculation_metadata(data.frame(chk$cells, donor = as.character(reader$obs("donor_id")), stringsAsFactors = FALSE))
  if (!identical(as.character(reader$obs("Study")), md$study)) stop("NRV 세포표 study 가 h5 obs/Study 와 다르다.")
  ck$add("metadata_matches_h5", TRUE, "block == R24 nrv_blocks(obs/donor_id) and study == obs/Study for every cell")
  n <- nrow(md); rm(c92, c51)
  he <- vapply(st, function(s) .nrv_sha(s[["edges"]]), ""); e <- readRDS(st$p92[["edges"]])
  if (he[["p92"]] != he[["p51"]] && !identical(e, readRDS(st$p51[["edges"]]))) stop("p92 와 p51 의 edges.rds 가 다르다 (같은 그래프여야 한다).")
  if (nrow(e) && any(e$i < 1L | e$j > n | e$i >= e$j)) stop("edge 가 1..n 상삼각 세포 번호가 아니다.")
  ck$add("edges_shared_by_partitions", TRUE, sprintf("%d positive upper-triangle edges, identical for p92 and p51 (%s)", nrow(e), if (he[["p92"]] == he[["p51"]]) "same file hash" else "same content"))
  # ── 원 RNA 근거 (raw 1 회) ──
  genes <- unique(c(ref$stage_markers$gene, ref$changes$gene)); gm <- nrv_gene_map(reader$symbols(), genes)
  spot <- vapply(fdt_cluster_order_safe(md$cluster92), function(k) which(md$cluster92 == k)[1L], 0L)
  groups <- list(p92 = data.frame(cluster = md$cluster92, block = md$block, stringsAsFactors = FALSE), p51 = data.frame(cluster = md$cluster51, block = md$block, stringsAsFactors = FALSE))
  ge <- nso_gene_evidence(.nrv_checked_chunks(reader$raw_chunks(a$chunk), n), gm, groups, keep_cells = spot)
  log("RNA streamed:", n, "cells;", sum(ge$genes$in_input), "/", length(genes), "reference genes in input")
  pe <- nso_program_evidence(ge, ref)
  den <- tapply(pe$groups$n_cells, paste(pe$groups$partition, pe$groups$program_id), sum)
  ck$add("denominators_preserved", all(den == n), sprintf("every partition x program sums cluster x block cells to %d (%d sums)", n, length(den)))
  # ── 관측 연결 · B 내부 연속성 ──
  oc <- list(); icn <- list()
  for (p in names(NSO_NRV)) { cl <- if (p == "p92") md$cluster92 else md$cluster51
    oc[[p]] <- nso_observed_connections(e, data.frame(cluster = cl, block = md$block, study = md$study, stringsAsFactors = FALSE), ge, programs = pe$programs)
    same <- .nso_same_pairs(oc[[p]]$pairs, readRDS(st[[p]][["pairs"]])$full)
    ck$add(paste0("pairs_match_stored_nrv_", p), same, sprintf("recomputed R24 pair table (%d observed pairs) == NRV structure pairs.rds$full", nrow(oc[[p]]$pairs)))
    if (!same) stop("다시 계산한 군집 쌍 표가 NRV 저장본과 다르다: ", p)
    icn[[p]] <- nso_internal_continuity(e, cl)
    log(p, "observed pairs:", nrow(oc[[p]]$pairs), "| triples:", nrow(icn[[p]]$triples), "| internal_disconnected:", sum(icn[[p]]$triples$status == "internal_disconnected")) }
  bnd <- nso_boundary_cells(e, md$cluster92); rm(e)
  # ── 참조 열 (R57 · R54; 원저자 열은 읽지 않는다) ──
  cmp <- nso_reference_comparison(a$cma, a$fdt, md)
  ck$add("comparison_inputs_checked", TRUE, "R57 cluster tables (author columns never loaded) and R54 lossless.rds checked against their manifest.json; R54 cluster sizes == cell table")
  # ── 계산 표 → calc_manifest.json (여기서 결과를 고정한다) ──
  files <- character(0); rtp <- logical(0)
  wt <- function(df, rel) { p0 <- file.path(out, rel); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); rtp[[rel]] <<- mbvt_write_tsv(df, p0)
    files <<- c(files, p0); invisible(p0) }
  bind_p <- function(x, what, empty) { z <- lapply(names(x), function(p) { d <- what(x[[p]]); if (nrow(d)) cbind(partition = p, d, stringsAsFactors = FALSE) else NULL })
    z <- do.call(rbind, z); if (is.null(z)) empty else { rownames(z) <- NULL; z } }
  wt(nso_availability(ref, ge), "stage_reference_availability.tsv")
  wt(nso_cell_table(md, ge, pe, cmp$tables), "cell_stage_evidence.tsv.gz")
  wt(pe$groups, "cluster_block_programs.tsv.gz"); wt(ge$groups, "cluster_block_genes.tsv.gz"); wt(pe$clusters, "cluster_programs.tsv"); wt(pe$programs, "program_definitions.tsv")
  e0 <- function(d) data.frame(partition = character(0), d[0, , drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)   # 모든 분할이 비었을 때도 열은 남긴다
  wt(bind_p(oc, function(z) z$pairs, e0(oc$p92$pairs)), "observed_pairs.tsv")
  wt(bind_p(oc, function(z) z$features, e0(.nso_feature_template())), "observed_edge_evidence.tsv.gz")
  wt(bind_p(icn, function(z) z$triples, e0(icn$p92$triples)), "internal_continuity.tsv.gz")
  wt(bind_p(icn, function(z) z$clusters, e0(icn$p92$clusters)), "internal_cluster_summary.tsv")
  wt(nso_limitations(ref, ge, pe, oc, icn), "stage_mapping_limitations.tsv")
  wt(do.call(rbind, unname(cmp$tables)), "cluster_reference_comparison.tsv"); wt(nso_spot_table(md, ge), "spot_check_cells.tsv")
  wt(ref$stages, "reference_stages.tsv"); wt(ref$scope, "reference_scope.tsv"); wt(ref$shared_markers, "reference_shared_markers.tsv"); wt(ref$transitions, "reference_transitions.tsv")
  p0 <- file.path(out, "lossless.rds")
  saveRDS(list(contract = "neut_stage_observation_2026-10-19/codex_instruction.md", settings = list(min_cells = NSO_MIN_CELLS, chunk = a$chunk, competing = NSO_COMPETING),
               reference = ref, metadata = md, gene_evidence = ge, programs = pe, connections = oc, internal = icn, comparison = cmp$tables), p0, version = 3)
  files <- c(files, p0)
  rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  writeLines(jsonlite::toJSON(list(phase = "calculation tables fixed before the original UMAP and author labels are read", written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), file.path(out, "calc_manifest.json"))
  calc_files <- files; log("calculation fixed:", length(files), "files")
  # ── 표시 (고정 뒤): 원 UMAP ──
  U <- reader$umap()
  if (!is.matrix(U) || nrow(U) != n || ncol(U) < 2L || any(!is.finite(U[, 1:2]))) stop("UMAP 세포 수 · 모양이 세포 표와 다르다 (", if (is.matrix(U)) paste(dim(U), collapse = " x ") else class(U)[1], " vs ", n, ").")
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  f1 <- file.path(fig, "fig01_neut_stage_evidence_umap.png"); f2 <- file.path(fig, "fig02_observed_connections_umap.png")
  t1 <- nso_draw_stage_evidence(f1, U[, 1:2], pe, ref, ko, fam); t2 <- nso_draw_connections(f2, U[, 1:2], md$cluster92, oc$p92$pairs, icn$p92$clusters, bnd, ko, fam)
  files <- c(files, f1, f2); wt(rbind(t1, t2), "figure_text.tsv"); log("figures written")
  # ── 사후 (고정 뒤): 원저자 라벨 ──
  au <- reader$author(); if (length(au) != n) stop("원저자 라벨 수가 세포 수와 다르다: ", length(au), " vs ", n)
  wt(nso_posthoc_composition(md, au), "posthoc_author_cluster_composition.tsv"); wt(nso_posthoc_programs(au, pe), "posthoc_author_program_summary.tsv")
  calc_now <- vapply(calc_files, .nrv_sha, ""); cm0 <- jsonlite::fromJSON(file.path(out, "calc_manifest.json"))$files
  ck$add("author_and_umap_read_after_fix", identical(unname(calc_now), unname(unlist(cm0[rel(calc_files)]))),
         "UMAP and author labels were opened only after calc_manifest.json was written; fixed calculation files unchanged afterwards")
  # ── provenance ──
  allinp <- c(rf, f, unlist(lapply(names(st), function(p) stats::setNames(st[[p]], paste0("nrv_", p, "_", names(st[[p]]))))), cmp$inputs)
  sh <- vapply(allinp, function(x) if (x %in% f) unname(cur[match(x, f)]) else .nrv_sha(x), "")
  wt(data.frame(input = names(allinp), path = unname(allinp), sha256 = unname(sh), stringsAsFactors = FALSE), "input_sha256.tsv")
  wt(snap, "source_snapshot.tsv")
  pk <- c("Matrix", "hdf5r", "igraph", "digest", "jsonlite")
  wt(data.frame(item = c("R_executable", "R", "platform", "libPaths", "packages", "title_font", "chunk", "min_cells", "random_seed", "command", "working_dir", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pk, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), if (ko) MBVR_TITLE_FAMILY else "default",
                          as.character(a$chunk), as.character(NSO_MIN_CELLS), "none (no random numbers are drawn)", paste(c("60_run_nso.R", argv), collapse = " "), getwd(),
                          format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), stringsAsFactors = FALSE), "environment.tsv")
  ck$add("tables_roundtrip_ok", all(rtp), sprintf("%d TSVs re-read with the same rows, field counts and text columns (mbvt_write_tsv)", length(rtp)))
  cks <- ck$table(); wt(cks, "checks.tsv")
  sm <- function(p) { z <- pe$clusters[pe$clusters$partition == p & pe$clusters$program_id %in% paste0("marker_", NSO_STAGES), ]; length(unique(z$cluster)) }
  nd <- colSums(pe$cells$n_detected > 0, na.rm = TRUE)
  rd <- c("# 호중구 성숙 단계 관측 — 첫 구현 묶음 결과", "",
          "- 계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md · 계획 2026-10-19-neutrophil-stage-observation.md. 설명서: N 루트 README_neut_stage_observation_ko.md.",
          sprintf("- 세포 %s (p92 군집 %d · p51 군집 %d, 전체 분모 유지). 참고 유전자 %d / %d 입력에 있음.", .nso_fmt_n(n), sm("p92"), sm("p51"), sum(ge$genes$in_input), length(genes)),
          paste0("- 단계 marker 프로그램 ≥1 검출 세포: ", paste(sprintf("%s %s", NSO_STAGES, .nso_fmt_n(nd[paste0("marker_", NSO_STAGES)])), collapse = " · "), "."),
          sprintf("- 관측 군집 쌍: p92 %d · p51 %d (NRV 저장본과 같음). B 내부 단절 삼중쌍: p92 %d / %d · p51 %d / %d.", nrow(oc$p92$pairs), nrow(oc$p51$pairs),
                  sum(icn$p92$triples$status == "internal_disconnected"), nrow(icn$p92$triples), sum(icn$p51$triples$status == "internal_disconnected"), nrow(icn$p51$triples)),
          "- 예측 stage_label 은 만들지 않았다. 값은 서술값이고 단계 정답 확률이 아니다. 원저자 라벨과 UMAP 은 calc_manifest.json 으로 계산 표를 고정한 뒤에만 읽었다.",
          "- 그림: figures/fig01_neut_stage_evidence_umap.png (6 단계 marker 프로그램), figures/fig02_observed_connections_umap.png (관측 연결 · 경계 세포 · 내부 단절).")
  p0 <- file.path(out, "README_ko.md"); writeLines(rd, p0); files <- c(files, p0)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .nso_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.nso_bb, .f))
  quit(status = nso_main(commandArgs(trailingOnly = TRUE)))
}
