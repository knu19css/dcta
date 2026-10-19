# =============================================================================
# 72_run_gse253355.R — GSE253355 성숙 호중구 포함 골수에 고정 atlas 6x10 패널을 군집 단위로 다시 적용하는 실행기
# =============================================================================
#
# 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_gse253355_instruction.md + codex_gse253355_local_addendum.md + codex_gse253355_review.md.
#   Ruling: 같은 폴더 claude_gse253355_progress.md
#   순서: (1) 입력 sha256 고정값 · 패널 = attempt_02 source 고정본(source_freeze.json sha256 고정) → (2) 객체 점검 기록(class · 슬롯 · assay · meta.data 열 · reduction · graph)
#         → (3) 원 그래프 · 군집 · 좌표가 없으면(NEW) 새 표현 그래프(HVG/PCA/kNN/SNN/Louvain/DrL, 고정 파라미터)를 점수 전에 만든다
#         → (4) raw count 점수 · 군집 점수 · 후보 · 단계별 최고 군집 · 전체 그래프 연결 감사 → calc_manifest.json
#         → (5) 그 뒤에만 저자 라벨: 구성표 · 성숙 호중구 확인 · 과립구 라벨 유도 그래프 경로 · 그림(PNG + PDF) · 요약(25 줄 이하) · 공개 사본(집계만).
#   결과 폴더 · 공개 폴더가 비어 있지 않으면 쓰지 않는다(exit 6). .part 입력은 쓰지 않는다.
#
# 실행 (W 에서, dcta 환경만)
#   run-r dcta R/branch_biology/72_run_gse253355.R --out-dir results/neutrophil_atlas_rna_2026-10-19/gse253355_01 --cluster-col NEW --graph NEW --umap NEW \
#     --donor-col donor --author-col cluster_anno_l2 --granulocyte-labels "<라벨;…>" --neutrophil-labels "Neutrophil" [--input …] [--public-dir …]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

GSX_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                 "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                 "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R", "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R")
GSX_INPUT <- "data/external/gse253355/S2_steady_GSE253355.local_copy.rds"                       # W 기준 (Codex addendum)
GSX_PANELS <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/atlas_stage_markers.csv"
GSX_PANEL_ATTEMPT <- "results/neutrophil_atlas_rna_2026-10-19/attempt_02"
GSX_HGNC <- "data/external/granulopoiesis_atlas_2026-10-19/hgnc_complete_set.txt"
GSX_PUBLIC <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_01"
GSX_PINNED <- list(input = "83fd17576c48e710e5a03a60be485beedbc088f82c6a97fa79d7b80d5367186a",       # local_copy.rds (Codex addendum)
                   source_freeze = "cdb4874a40946e60b3183d0e0ee582a43d35e3373fa7350a487725bfcaec6a00")  # attempt_02 source_freeze.json (지시서)
GSX_PUBLIC_FILES <- c("summary_ko.md", "object_inspection.json", "object_metadata_columns.csv", "panel_gene_map.csv", "panel_gene_detection.csv", "new_graph_parameters.csv",
                      "graph_checks.csv", "cluster_scores.csv", "cluster_candidates.csv", "stage_peak_clusters.csv", "cluster_edges.csv", "peak_paths.csv",
                      "posthoc_author_composition.csv", "posthoc_donor_composition.csv", "posthoc_mature_neutrophil.csv", "posthoc_peak_paths_granulocyte.csv",
                      "posthoc_peak_path_identities.csv", "checks.csv", "environment.csv", "figure_text.csv")


#' 긴 군집 점수 표 → 군집 x 단계 행렬
.gsx_wide <- function(cl, col) { u <- unique(cl$cluster); M <- matrix(NA_real_, length(u), length(NSO_STAGES), dimnames = list(u, NSO_STAGES))
  M[cbind(match(cl$cluster, u), match(cl$stage, NSO_STAGES))] <- cl[[col]]; M }

#' 군집 후보 표: 주 점수(고정 분모, 모든 세포) 순위 + donor 동일 가중 · 관측 유전자 분모 진단 순위, 작은 군집 경고, 불완전 패널 표시
gsx_candidates <- function(clus, tol, min_cells, incomplete_genes = "") {
  S <- .gsx_wide(clus, "mean_all_cells"); E <- .gsx_wide(clus, "mean_donor_equal"); A <- .gsx_wide(clus, "mean_available_genes")
  rp <- atlas_cl_rank(S, tol); re <- atlas_cl_rank(E, tol); ra <- atlas_cl_rank(A, tol); z <- clus[clus$stage == NSO_STAGES[1], ]; n <- z$n_cells[match(rownames(S), z$cluster)]
  d <- data.frame(cluster = rownames(S), n_cells = as.integer(n), n_donors = z$n_donors[match(rownames(S), z$cluster)], small_cluster = n < min_cells, stringsAsFactors = FALSE)
  d$warning <- ifelse(d$small_cluster, sprintf("n_cells %d < %d: kept, not merged or hidden", d$n_cells, as.integer(min_cells)), "")
  for (s in NSO_STAGES) d[[paste0("score_", s)]] <- unname(S[, s])
  d <- cbind(d, rp[, setdiff(names(rp), "cluster")])
  if (nzchar(incomplete_genes)) d$comparison_scope[d$comparison_scope != "none"] <- "incomplete_panel_comparison"
  d$missing_panel_genes <- rep(incomplete_genes, nrow(d))
  for (s in NSO_STAGES) d[[paste0("donor_equal_", s)]] <- unname(E[, s])
  d$donor_equal_top_status <- re$score_top_status; d$donor_equal_top_candidate <- re$score_top_candidate
  pos <- c("score_top_candidate", "tied")
  d$donor_equal_agreement <- ifelse(re$score_top_status == "unavailable", "donor_equal_unavailable",
                                    ifelse(d$score_top_status %in% pos & re$score_top_status %in% pos, ifelse(d$score_top_candidate == re$score_top_candidate, "same_top", "different_top"),
                                           ifelse(d$score_top_status == "unassigned_zero" & re$score_top_status == "unassigned_zero", "both_zero", "not_comparable")))
  d$available_gene_top_candidate <- ra$score_top_candidate; d$available_gene_agreement <- ifelse(ra$score_top_candidate == d$score_top_candidate, "same_top", "different_top")
  d$note <- rep("score_top_candidate = highest cluster-mean panel score (fixed 10-gene denominator); not a probability, confirmed cell type or differentiation call", nrow(d)); d
}

gsx_main <- function(argv, root = NULL, pins = GSX_PINNED) {
  keys <- c("--input" = "input", "--panels" = "panels", "--panel-attempt" = "panel_attempt", "--hgnc" = "hgnc", "--out-dir" = "out", "--public-dir" = "public",
            "--assay" = "assay", "--cluster-col" = "cluster_col", "--donor-col" = "donor_col", "--author-col" = "author_col", "--umap" = "umap", "--graph" = "graph",
            "--granulocyte-labels" = "gran", "--neutrophil-labels" = "neut", "--n-hvg" = "n_hvg", "--n-pcs" = "n_pcs", "--k" = "k", "--resolution" = "resolution", "--seed" = "seed")
  a <- list(); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  req <- c("out", "cluster_col", "donor_col", "author_col", "umap", "graph", "gran", "neut")
  if (any(vapply(req, function(k) is.null(a[[k]]) || is.na(a[[k]]) || !nzchar(a[[k]]), TRUE))) {
    cat("--out-dir, --cluster-col, --donor-col, --author-col, --umap, --graph, --granulocyte-labels, --neutrophil-labels 가 필요하다.\n"); return(2L) }
  if (is.null(root)) root <- mdu_find_root()
  dflt <- list(input = GSX_INPUT, panels = GSX_PANELS, panel_attempt = GSX_PANEL_ATTEMPT, hgnc = GSX_HGNC, public = GSX_PUBLIC)
  for (k in names(dflt)) if (is.null(a[[k]])) a[[k]] <- file.path(root, dflt[[k]])
  if (is.null(a$assay)) a$assay <- "RNA"
  gp <- list(n_hvg = 2000L, n_pcs = 30L, k = 20L, prune = 1 / 15, resolution = 0.5, seed = 42L, scale_max = 10, loess_span = 0.3)   # = R70 GSX_GRAPH_PARAMS (R70 은 아직 적재 전)
  for (k in c("n_hvg", "n_pcs", "k", "seed")) if (!is.null(a[[k]])) gp[[k]] <- suppressWarnings(as.integer(a[[k]]))
  if (!is.null(a$resolution)) gp$resolution <- suppressWarnings(as.numeric(a$resolution))
  if (any(vapply(gp, function(x) is.na(x) || x <= 0, TRUE))) { cat("그래프 파라미터는 양수여야 한다.\n"); return(2L) }
  a$graph_params <- gp; out <- a$out; pub <- a$public
  for (d in c(out, pub)) if (dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", d, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", GSX_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(GSX_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = GSX_SOURCES)
    get("gsx_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

gsx_build <- function(a, out, snap, log, argv, pins) {
  ck <- .mbvq_checks(); t0 <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"); rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  files <- character(0); crt <- logical(0)
  wc <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  wj <- function(x, f) { p0 <- file.path(out, f); writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), p0); files <<- c(files, p0); invisible(p0) }
  gp <- a$graph_params
  # ── (1) 입력 · 패널 고정 ──
  if (grepl("\\.part$", a$input)) stop("다운로드 중인 .part 입력은 쓰지 않는다: ", a$input)
  if (!file.exists(a$input)) stop("입력 파일이 없다: ", a$input)
  if (is.null(pins$input)) stop("입력 hash 고정값이 없다 (Codex 기록의 sha256 이 필요하다).")
  hi <- .nrv_sha(a$input); if (!identical(hi, pins$input)) stop("입력 hash 가 고정값과 다르다: ", hi)
  ck$add("input_pinned", TRUE, sprintf("%s sha256 %s == pinned (%.0f bytes)", basename(a$input), substr(hi, 1, 12), file.size(a$input)))
  rdc <- function(p) utils::read.csv(p, colClasses = "character", na.strings = character(0), check.names = FALSE)
  pan <- rdc(a$panels); src <- rdc(file.path(a$panel_attempt, "source", "source_panels.csv")); kc <- c("stage", "rank", "ensg", "symbol")
  fz <- file.path(a$panel_attempt, "source_freeze.json"); hf <- .nrv_sha(fz); fzf <- jsonlite::fromJSON(fz)$files
  okp <- all(kc %in% names(pan)) && identical(pan[kc], src[kc]) && identical(fzf[["source/source_panels.csv"]], .nrv_sha(file.path(a$panel_attempt, "source", "source_panels.csv"))) &&
    (is.null(pins$source_freeze) || identical(hf, pins$source_freeze)) && all(pan$stage %in% NSO_STAGES)
  if (!okp) stop("패널이 attempt_02 source 고정본과 다르다 (stage · rank · ENSG · symbol · source_freeze hash).")
  pan$rank <- as.integer(pan$rank); pan <- pan[order(match(pan$stage, NSO_STAGES), pan$rank), kc]
  ck$add("panels_frozen_source", TRUE, sprintf("%d panel genes (%s) identical to attempt_02 source_panels.csv; source_freeze.json %s lists that file", nrow(pan),
                                               paste(names(table(factor(pan$stage, levels = NSO_STAGES))), as.integer(table(factor(pan$stage, levels = NSO_STAGES))), sep = "=", collapse = " "), substr(hf, 1, 12)))
  hg <- utils::read.delim(a$hgnc, quote = "", colClasses = "character", stringsAsFactors = FALSE); miss <- setdiff(c("symbol", "prev_symbol", "alias_symbol"), names(hg))
  if (length(miss)) stop("HGNC 표에 열이 없다: ", paste(miss, collapse = ", ")); hg <- hg[, c("symbol", "prev_symbol", "alias_symbol")]
  # ── (2) 객체 점검 · 추출 ──
  ro <- gsx_open_rds(a$input); obj <- ro$object; insp <- gsx_inspect(obj); md <- gsx_meta(obj)
  wj(c(list(input = basename(a$input), sha256 = hi, read_format = ro$format), insp), "object_inspection.json"); wc(gsx_meta_columns(md), "object_metadata_columns.csv")
  log("object read:", ro$format, "|", insp$object_class, "|", nrow(md), "cells | assays:", paste(vapply(insp$assays, `[[`, "", "name"), collapse = ","),
      "| graphs:", length(insp$graphs), "| reductions:", length(insp$reductions))
  cn <- gsx_counts(obj, a$assay); X <- cn$counts
  ck$add("raw_counts_used", TRUE, sprintf("%s: %d genes x %d cells, non-negative integers; data / scale.data not used", paste(cn$source, collapse = "+"), nrow(X), ncol(X)))
  new_g <- identical(a$graph, "NEW"); new_c <- identical(a$cluster_col, "NEW"); new_u <- identical(a$umap, "NEW")
  if ((new_c || new_u) && !new_g) stop("새 군집 · 새 좌표는 새 그래프에서만 만든다 (--graph NEW).")
  G <- if (new_g) NULL else gsx_graph(obj, a$graph); U <- if (new_u) NULL else gsx_umap(obj, a$umap)[, 1:2, drop = FALSE]; rm(obj, ro); invisible(gc())
  cells <- colnames(X); n <- length(cells)
  for (k in c(a$donor_col, a$author_col, if (!new_c) a$cluster_col)) if (!k %in% names(md)) stop("meta.data 에 열이 없다: ", k)
  donor <- as.character(md[[a$donor_col]]); if (anyNA(donor) || any(!nzchar(donor))) stop("donor 열에 결측이 있다.")
  # ── (3) 새 표현 그래프 (원 그래프 · 군집 · 좌표가 없을 때만, 점수 전에) ──
  if (new_g) { t1 <- Sys.time()
    hv <- gsx_hvg(X, gp$n_hvg, gp$loess_span); log("new graph: HVG", length(hv$genes)); pc <- gsx_pca(X, hv$genes, gp$n_pcs, gp$scale_max); log("new graph: PCA", pc$n_pcs)
    nn <- gsx_knn(pc$scores, gp$k); log("new graph: kNN k =", ncol(nn)); G <- gsx_snn(nn, gp$prune); dimnames(G) <- list(cells, cells); log("new graph: SNN stored", length(G@x))
    clusters <- if (new_c) gsx_louvain(G, gp$resolution, gp$seed) else as.character(md[[a$cluster_col]])
    if (new_u) { U <- gsx_layout(G, gp$seed); rownames(U) <- cells }
    log("new graph done:", length(unique(clusters)), "clusters;", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min")
    wc(data.frame(name = c("graph_source", "hvg_method", "n_hvg", "pca_method", "n_pcs", "scale_max", "knn", "k", "snn", "prune", "clustering", "resolution", "seed", "display_layout",
                           "panels_or_labels_used"),
                  value = c("new_expression_snn (not original; original graph / reductions / seurat_clusters absent in the object)", "vst-equivalent on raw counts (Matrix re-implementation of Seurat vst formula, loess span 0.3, clip sqrt(N))",
                            as.character(length(hv$genes)), "log1p(1e4*count/library) of HVG, gene-wise centred and scaled (upper clip), covariance eigen decomposition, sign-fixed",
                            as.character(pc$n_pcs), as.character(gp$scale_max), "exact Euclidean on PCs, self first", as.character(ncol(nn)), "Jaccard of self-including neighbour sets",
                            format(gp$prune), if (new_c) "igraph cluster_louvain on SNN (diagonal excluded), clusters renamed 0.. by size" else paste("existing column", a$cluster_col),
                            format(gp$resolution), as.character(gp$seed), if (new_u) "igraph layout_with_drl on SNN (NOT UMAP; uwot absent in dcta); display only" else paste("existing reduction", a$umap),
                            "none (no stage panels, author labels or DAG used to build the graph)"), stringsAsFactors = FALSE), "new_graph_parameters.csv")
    wc(hv$table[hv$table$selected, c("gene", "mean", "variance", "variance_expected", "variance_standardized")], "new_graph_hvg.csv")
    ck$add("new_graph_built_before_scores", TRUE, sprintf("new expression-only graph built before any stage score: %d HVG, %d PCs, k %d, prune %.4f, Louvain resolution %s seed %d -> %d clusters; display = %s",
                                                         length(hv$genes), pc$n_pcs, ncol(nn), gp$prune, format(gp$resolution), gp$seed, length(unique(clusters)), if (new_u) "DrL layout (not UMAP)" else a$umap))
    gsrc <- "new_expression_snn (not original)"; csrc <- if (new_c) "new Louvain clusters (not author clusters)" else a$cluster_col
    usrc <- if (new_u) (if (exists("mbvr_hangul_font") && mbvr_hangul_font()) "새 표시 좌표 = 새 SNN 의 DrL 그래프 레이아웃 (UMAP 아님, 원 좌표 아님)" else "new display = DrL layout of the new SNN (not UMAP, not original)") else a$umap
    gbuild <- list(hvg = hv$genes, pca_sdev = pc$sdev, pca_scores = pc$scores, knn = nn); rm(hv, pc, nn); invisible(gc())
  } else { clusters <- as.character(md[[a$cluster_col]]); gsrc <- paste0("original object graph ", a$graph); csrc <- a$cluster_col; usrc <- paste0("object reduction ", a$umap); gbuild <- NULL
    ck$add("original_graph_used", TRUE, sprintf("graph %s, clusters %s, coordinates %s taken from the object as stored", a$graph, a$cluster_col, a$umap)) }
  if (anyNA(clusters) || any(!nzchar(clusters))) stop("군집이 비어 있는 세포가 있다.")
  gsx_check_ids(cells, meta = md, umap = U, graph = G)
  ck$add("cell_ids_aligned", TRUE, sprintf("%d cell IDs identical in order: counts, meta.data, coordinates, graph rows / columns", n))
  gc0 <- gsx_graph_checks(G); wc(data.frame(item = c("graph_source", "cluster_source", "coordinate_source", names(gc0)), value = c(gsrc, csrc, usrc, vapply(gc0, function(z) as.character(z)[1], "")),
                                            stringsAsFactors = FALSE), "graph_checks.csv")
  if (!gc0$square || !gc0$nonnegative) stop("그래프가 정사각 · 음수 없음이 아니다.")
  # ── (4) 점수 · 후보 · 최고 군집 · 전체 그래프 감사 (저자 라벨 없음) ──
  pm <- gsx_panel_map(rownames(X), pan, hg); cs <- gsx_cell_scores(X, pm)
  if ("nCount_RNA" %in% names(md)) { okl <- isTRUE(all.equal(cs$library, as.numeric(md$nCount_RNA), tolerance = 0))
    ck$add("library_matches_nCount_RNA", okl, "per-cell all-gene library (colSums of raw counts) equals meta.data nCount_RNA (same cells, same order)"); if (!okl) stop("library 가 nCount_RNA 와 다르다.") }
  sg <- cs$stage_genes; inc <- paste(pm$symbol[pm$status != "exact"], collapse = "|")
  ck$add("panel_genes_measured", TRUE, sprintf("%d / %d panel genes measured by exact symbol; missing: %s", sum(pm$status == "exact"), nrow(pm), if (nzchar(inc)) inc else "none"))
  clus <- gsx_cluster_scores(cs$fixed, cs$n_detected, clusters, donor, NSO_MIN_CELLS, A = cs$available)
  okd <- all(tapply(clus$n_cells, clus$stage, sum) == n); ck$add("denominators_preserved", okd, sprintf("every stage sums cluster sizes to %d cells; all clusters kept (small clusters warned)", n))
  tt <- as.numeric(NRV_THRESHOLDS$value[NRV_THRESHOLDS$name == "tie_tol"]); ck$add("tie_tol_matches_r23", identical(tt, ATLAS_CL_TIE_TOL), sprintf("tie tolerance %g; zero = exactly 0", ATLAS_CL_TIE_TOL))
  cand <- gsx_candidates(clus, ATLAS_CL_TIE_TOL, NSO_MIN_CELLS, inc)
  S1 <- .gsx_wide(clus, "mean_all_cells"); E1 <- .gsx_wide(clus, "mean_donor_equal"); nnc <- stats::setNames(cand$n_cells, cand$cluster)
  pk <- atlas_cl_stage_peaks(S1, cand[, c("cluster", "score_top_candidate")], ATLAS_CL_TIE_TOL, nnc)
  pe <- atlas_cl_stage_peaks(E1, data.frame(cluster = cand$cluster, score_top_candidate = cand$donor_equal_top_candidate, stringsAsFactors = FALSE), ATLAS_CL_TIE_TOL, nnc)
  peaks <- data.frame(pk, donor_equal_peak_status = pe$peak_status, donor_equal_peak_clusters = pe$peak_clusters,
                      comparison_note = rep("stage peak = highest cluster for this stage among all clusters; differs from the within-cluster top", nrow(pk)), stringsAsFactors = FALSE)
  e <- gsx_upper_edges(G); ce <- gsx_cluster_edges(e, clusters, donor); comp <- gsx_components(e, n)
  pp <- gsx_peak_paths(e, clusters, donor, peaks[, c("stage", "peak_clusters")], n); full_paths <- attr(pp, "cell_paths")
  log("scores:", nrow(cand), "clusters |", paste(names(table(cand$score_top_status)), table(cand$score_top_status), sep = "=", collapse = " "), "| peaks:",
      paste(peaks$stage, ifelse(nzchar(peaks$peak_clusters), peaks$peak_clusters, peaks$peak_status), collapse = " "), "| graph components:", comp$n_components)
  wc(pm, "panel_gene_map.csv"); wc(gsx_gene_detection(pm, cs), "panel_gene_detection.csv"); wc(sg, "panel_stage_genes.csv")
  wc(clus, "cluster_scores.csv"); wc(cand, "cluster_candidates.csv"); wc(peaks, "stage_peak_clusters.csv")
  wc(ce, "cluster_edges.csv"); wc(pp, "peak_paths.csv")
  wc(data.frame(graph = "full", n_cells = n, n_edges = nrow(e), n_components = comp$n_components, largest_component = comp$largest, n_isolated = comp$n_isolated, stringsAsFactors = FALSE),
     "graph_components.csv")
  p0 <- file.path(out, "lossless.rds")
  saveRDS(list(contract = "codex_gse253355_instruction.md + local_addendum + review", input_sha256 = hi, cells = cells, clusters = clusters, donor = donor, library = cs$library,
               fixed = cs$fixed, available = cs$available, n_detected = cs$n_detected, coordinates = U, graph_build = gbuild, graph_params = gp, peak_cell_paths_full = full_paths,
               scores = clus, candidates = cand, peaks = peaks, panel_map = pm), p0, version = 3); files <- c(files, p0)
  cmf <- file.path(out, "calc_manifest.json")
  writeLines(jsonlite::toJSON(list(phase = "scores, candidates, stage peaks and full-graph connectivity fixed before author labels are read", written_at = stamp(), input_sha256 = hi,
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  calc_files <- files; log("calculation fixed:", length(files), "files")
  # ── (5) 사후 (고정 뒤): 저자 라벨 ──
  au <- as.character(md[[a$author_col]]); au[is.na(au) | !nzchar(au)] <- "(no author label)"
  gl <- trimws(strsplit(a$gran, ";", fixed = TRUE)[[1]]); nl <- trimws(strsplit(a$neut, ";", fixed = TRUE)[[1]])
  ml <- setdiff(c(gl, nl), au); if (length(ml)) stop("지정한 저자 라벨이 자료에 없다 (철자 확인): ", paste(ml, collapse = ", "))
  ck$add("granulocyte_labels_found", TRUE, sprintf("granulocyte-associated labels (%s) and mature-neutrophil labels (%s) exist in %s", paste(gl, collapse = "; "), paste(nl, collapse = "; "), a$author_col))
  comp_tab <- do.call(rbind, lapply(cand$cluster, function(k) { w <- clusters == k; tb <- sort(table(au[w]), decreasing = TRUE)
    data.frame(cluster = k, author_label = names(tb), n_cells = as.integer(tb), frac_of_cluster = as.numeric(tb) / sum(w), stringsAsFactors = FALSE) }))
  comp_tab$note <- rep("post hoc after calc_manifest.json; not used for any score, candidate or threshold; not accuracy", nrow(comp_tab)); wc(comp_tab, "posthoc_author_composition.csv")
  atop <- vapply(cand$cluster, function(k) { z <- comp_tab[comp_tab$cluster == k, ]; z$author_label[1] }, ""); names(atop) <- cand$cluster
  dtab <- as.data.frame(table(cluster = clusters, donor = donor), stringsAsFactors = FALSE); dtab <- dtab[dtab$Freq > 0, ]; names(dtab)[3] <- "n_cells"
  dtab <- dtab[order(match(dtab$cluster, cand$cluster), dtab$donor), ]; rownames(dtab) <- NULL; wc(dtab, "posthoc_donor_composition.csv")
  mt <- do.call(rbind, lapply(nl, function(l) { w <- au == l; tb <- sort(table(clusters[w]), decreasing = TRUE); top <- utils::head(names(tb), 5)
    tc <- cand$score_top_candidate[match(clusters[w], cand$cluster)]; ft <- table(factor(tc, levels = c(NSO_STAGES, "")))
    data.frame(author_label = l, n_cells = sum(w), n_donors = length(unique(donor[w])), donors = paste(sort(unique(donor[w])), collapse = "|"),
               main_clusters = paste(sprintf("%s:%d(%.2f)", top, as.integer(tb[top]), as.numeric(tb[top]) / sum(w)), collapse = " "),
               main_cluster_top_candidates = paste(cand$score_top_candidate[match(top, cand$cluster)], collapse = " "),
               t(stats::setNames(as.numeric(ft[NSO_STAGES]) / sum(w), paste0("frac_cells_in_cluster_top_", NSO_STAGES))), stringsAsFactors = FALSE, check.names = FALSE) }))
  wc(mt, "posthoc_mature_neutrophil.csv")
  gmask <- au %in% gl; gp2 <- gsx_peak_paths(e, clusters, donor, peaks[, c("stage", "peak_clusters")], n, keep = gmask); gpaths <- attr(gp2, "cell_paths")
  gcomp <- gsx_components(e, n, keep = gmask); gp2$induced_by <- rep(paste0(a$author_col, " in {", paste(gl, collapse = "; "), "} (post hoc subset)"), nrow(gp2))
  wc(gp2, "posthoc_peak_paths_granulocyte.csv")
  ident <- function(p) if (!nzchar(p)) "" else paste(vapply(strsplit(p, ">", fixed = TRUE)[[1]], function(k) sprintf("%s(%s)", k, atop[[k]]), ""), collapse = ">")
  pid <- rbind(data.frame(graph = "full", pp[, c("from_stage", "to_stage", "relation", "quotient_status", "quotient_path", "cells_connected", "cell_hops", "cell_path_clusters", "path_agreement")],
                          stringsAsFactors = FALSE),
               data.frame(graph = "granulocyte_induced", gp2[, c("from_stage", "to_stage", "relation", "quotient_status", "quotient_path", "cells_connected", "cell_hops", "cell_path_clusters", "path_agreement")],
                          stringsAsFactors = FALSE))
  pid$quotient_path_identity <- vapply(pid$quotient_path, ident, ""); pid$cell_path_identity <- vapply(pid$cell_path_clusters, ident, "")
  pid$note <- rep("identities = most frequent author label of each cluster (post hoc); quotient path is not a cell trajectory", nrow(pid)); wc(pid, "posthoc_peak_path_identities.csv")
  gce <- gsx_cluster_edges(e[gmask[e$i] & gmask[e$j], , drop = FALSE], clusters, donor)
  # ── 그림 ──
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  clab <- if (new_u) (if (ko) "좌표 = 새 SNN 의 DrL 그래프 레이아웃 (UMAP 아님 · 원 좌표 아님 · 표시 전용)" else "coordinates = DrL layout of the new SNN (not UMAP, not original; display only)")
          else sprintf(if (ko) "좌표 = 객체의 %s (표시 전용)" else "coordinates = object reduction %s (display only)", a$umap)
  f1 <- gsx_draw_fig1(file.path(fig, "fig1_umap_author_vs_candidate"), U, clusters, au, cand, clab, ko, fam)
  f2 <- gsx_draw_fig2(file.path(fig, "fig2_granulocyte_zoom_heatmap"), U, clusters, gmask, S1, cand, peaks, gce, atop, gl, clab, ko, fam)
  files <- c(files, f1$files, f2$files); wc(rbind(f1$text, f2$text), "figure_text.csv"); log("figures written")
  calc_now <- vapply(calc_files, .nrv_sha, ""); cm0 <- jsonlite::fromJSON(cmf)$files
  ck$add("calc_fixed_before_author_labels", identical(unname(calc_now), unname(unlist(cm0[rel(calc_files)]))),
         "author labels were read only after calc_manifest.json; fixed calculation files unchanged afterwards")
  # ── 요약 · 기록 ──
  fmt <- function(x) .nca_fmt_n(x); cons <- pp[pp$consecutive, ]; gcons <- gp2[gp2$consecutive, ]
  ptxt <- function(z) paste(sprintf("%s>%s %s", z$from_stage, z$to_stage, ifelse(z$relation == "different_clusters", z$path_agreement, z$relation)), collapse = " · ")
  tbt <- table(factor(cand$score_top_candidate[cand$score_top_status == "score_top_candidate"], levels = NSO_STAGES))
  pkc <- peaks$peak_clusters; subtxt <- paste(sprintf("%s %s %s/%s", peaks$stage, pkc, vapply(pkc, function(z) if (nzchar(z)) fmt(sum(clusters %in% strsplit(z, "|", fixed = TRUE)[[1]] & gmask)) else "-", ""),
                                                      vapply(pkc, function(z) if (nzchar(z)) fmt(sum(clusters %in% strsplit(z, "|", fixed = TRUE)[[1]])) else "-", "")), collapse = " · ")
  topfr <- function(i) paste(sprintf("%s %.2f", NSO_STAGES, unlist(mt[i, paste0("frac_cells_in_cluster_top_", NSO_STAGES)])), collapse = " · ")
  sm <- c("# GSE253355 골수 — 고정 atlas 6 단계 패널의 군집 단위 재적용 (요약)", "",
          sprintf("- 입력: %s (sha256 %s…, Codex local copy, GEO 원본과 바이트 동일 아님). 세포 %s · 유전자 %s · donor %d (%s).", basename(a$input), substr(hi, 1, 12), fmt(n), fmt(nrow(X)),
                  length(unique(donor)), a$donor_col),
          sprintf("- 패널 유전자 측정: %d / %d (정확한 symbol). 없음: %s (점수는 고정 분모 10, 관측 유전자 분모는 진단).", sum(pm$status == "exact"), nrow(pm), if (nzchar(inc)) gsub("|", ", ", inc, fixed = TRUE) else "없음"),
          sprintf("- 그래프 · 군집 · 좌표: %s · %s · %s.", gsrc, csrc, if (new_u) "새 DrL 레이아웃(UMAP 아님)" else a$umap),
          sprintf("- 군집 %d 개 1 위 후보: %s · 동률 %d · 0 %d (작은 군집 경고 %d).", nrow(cand), paste(sprintf("%s %d", NSO_STAGES, as.integer(tbt)), collapse = " · "),
                  sum(cand$score_top_status == "tied"), sum(cand$score_top_status == "unassigned_zero"), sum(cand$small_cluster)),
          paste0("- 단계별 최고 군집 (사후 저자 최다 라벨): ", paste(sprintf("%s %s", peaks$stage, ifelse(nzchar(peaks$peak_clusters),
                  vapply(peaks$peak_clusters, function(z) paste(sprintf("%s(%s)", strsplit(z, "|", fixed = TRUE)[[1]], atop[strsplit(z, "|", fixed = TRUE)[[1]]]), collapse = "|"), ""), peaks$peak_status)), collapse = " · "), "."),
          paste0("- 저자 라벨: ", paste(sprintf("%s %s 세포 · donor %d · 주 군집 %s", mt$author_label, fmt(mt$n_cells), mt$n_donors, mt$main_clusters), collapse = "; "), " (성숙 호중구 = Neutrophil 라벨)."),
          paste0("- 그 세포들 군집의 1 위 후보 비율: ", paste(sprintf("%s [%s]", mt$author_label, vapply(seq_len(nrow(mt)), topfr, "")), collapse = "; "), "."),
          "- atlas PM 패널은 attempt_02 에서 기록한 후기 유전자(MMP9 · FCGR3B · FPR1 · AQP9 · CCL4L2)를 포함한다; Neutrophil 군집의 PM 1 위는 이 한계와 함께 읽는다.",
          paste0("- 연속 단계 최고 군집 연결 (전체 그래프): ", ptxt(cons), "."),
          paste0("- 같은 연결 (과립구 라벨 유도 그래프, 사후): ", ptxt(gcons), "."),
          paste0("- 유도 그래프 안 최고 군집 세포 (과립구 라벨/전체): ", subtxt, ". 소수 라벨 세포가 잇는 연결은 성숙 근거가 아니다."),
          "- cluster_path_only = 군집 요약에는 경로가 있으나 실제 세포 그래프에서는 끊김. 연결은 분화 증거가 아니고, 끊김만으로 생물학적 부재를 말하지 않는다.",
          "- 후보는 확률 · 확정 cell type · 분화 판정이 아니다. 저자 라벨은 계산 고정 뒤 사후 맥락으로만 썼다. 패널 유전자는 호중구 밖 계열에서도 발현된다.",
          "- 한계: 새 그래프 · 군집 · 좌표는 원 저자 그래프의 재현이 아니다. donor 동일 가중은 donor x 군집 >= 5 세포만(없으면 결측). 정확도 주장 없음.",
          "- 표: cluster_candidates · cluster_scores · stage_peak_clusters · panel_gene_map · peak_paths · posthoc_*.csv. 그림: figures/fig1 · fig2 (PNG · PDF).")
  p0 <- file.path(out, "summary_ko.md"); writeLines(sm, p0); files <- c(files, p0)
  allinp <- c(input = a$input, panels = a$panels, panel_source = file.path(a$panel_attempt, "source", "source_panels.csv"), panel_source_freeze = fz, hgnc = a$hgnc)
  wc(data.frame(input = names(allinp), path = unname(allinp), sha256 = vapply(seq_along(allinp), function(k) if (names(allinp)[k] == "input") hi else .nrv_sha(allinp[[k]]), ""),
                stringsAsFactors = FALSE), "input_sha256.csv")
  wc(snap, "source_snapshot.csv")
  pkg <- c("Matrix", "igraph", "digest", "jsonlite")
  wc(data.frame(item = c("R_executable", "R", "platform", "libPaths", "packages", "seurat_packages", "title_font", "graph_params", "min_cells", "random_seed", "command", "working_dir", "started_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pkg, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "),
                          paste(c("Seurat", "SeuratObject", "irlba", "RANN", "uwot"), vapply(c("Seurat", "SeuratObject", "irlba", "RANN", "uwot"), function(p) if (nzchar(system.file(package = p))) "present" else "absent", ""), sep = "=", collapse = " "),
                          if (ko) MBVR_TITLE_FAMILY else "default", paste(names(gp), unlist(gp), sep = "=", collapse = " "), as.character(NSO_MIN_CELLS), paste("set.seed", gp$seed, "before Louvain and DrL only"),
                          paste(c("72_run_gse253355.R", argv), collapse = " "), getwd(), format(t0, "%Y-%m-%dT%H:%M:%S%z"), stamp()), stringsAsFactors = FALSE), "environment.csv")
  ck$add("tables_roundtrip_ok", all(crt), sprintf("%d CSVs re-read with identical dimensions, names and cell text", length(crt)))
  pubf <- intersect(GSX_PUBLIC_FILES, rel(files))
  ck$add("public_copy_aggregate_only", !any(grepl("lossless|cell_level|hvg", pubf)), sprintf("%d aggregate files + figures copied to the public folder; no cell-level or matrix files", length(pubf)))
  cks <- ck$table(); wc(cks, "checks.csv"); allf <- c(files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), input_sha256 = hi, files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  pub <- a$public; dir.create(file.path(pub, "figures"), recursive = TRUE, showWarnings = FALSE)
  pf <- c(file.path(out, pubf), list.files(fig, full.names = TRUE)); dest <- ifelse(grepl("/figures/", pf), file.path(pub, "figures", basename(pf)), file.path(pub, basename(pf)))
  if (!all(file.copy(pf, dest, overwrite = FALSE))) stop("공개 사본 복사 실패")
  writeLines(jsonlite::toJSON(list(source_output = normalizePath(out), input_sha256 = hi, files = as.list(stats::setNames(vapply(dest, .nrv_sha, ""), substring(normalizePath(dest), nchar(normalizePath(pub)) + 2L)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(pub, "public_manifest.json"))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .gsx_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.gsx_bb, .f))
  quit(status = gsx_main(commandArgs(trailingOnly = TRUE)))
}
