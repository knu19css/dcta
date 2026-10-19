# =============================================================================
# 64_run_atlas_rna.R — atlas RNA 고정 패널의 target 적용 실행기: source 전용 고정 → target RNA 근거 → R62 형 원 edge 연결 → 기준선 비교 → 계산 고정
#                      → (그 뒤) 원 UMAP 그림 · 사후 원저자 표
# =============================================================================
#
# 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md, docs/experiments/neutrophil_atlas_rna_2026-10-19/design.md, plan Task 2.
#   Ruling: 같은 폴더 claude_progress.md
#   - 순서 (1) source 파일 hash = docs source_downloads.json (+ 고정값) → R63 reference → source/ CSV · RDS → source_freeze.json.
#          (2) 그 뒤에야 target 을 연다: h5ad · graph 계약 hash, O 의 NRV structure 기록(stage manifest), 기준선 NSO · connectivity manifest 고정값.
#          (3) 원 count 를 chunk 로 한 번 읽어 패널 · 변화 유전자 근거(R59 그대로) → R62 순서 도달 규칙 그대로의 원 edge 연결 → 기준선 비교 → calc_manifest.json.
#          (4) 고정 뒤에만 원 UMAP 그림과 원저자 posthoc 표. 예측 stage_label · 정확도는 만들지 않는다. 결과 폴더가 비어 있지 않으면 쓰지 않는다(exit 6).
#   - target 유전자는 ENSG 로 정확히 대응한다(raw/var/_index). 없는 유전자 = missing_gene (0 으로 바꾸지 않는다), 있으나 count 가 모두 0 = observed_zero.
#   - 판정불가 단계(빈 패널 · target 에 유전자 없음 · 기준보다 적은 관측 유전자)는 0 이 아니라 결측으로 두고, 그 단계 뒤의 순서 도달은 평가하지 않는다.
#   - 기준선 NSO 패널(단계당 2–4 유전자)과 atlas 패널(최대 10)은 크기가 달라, 검출 · 연결 증가를 단계 정확도 향상으로 읽지 않는다.
#   R05 · R17 · R18 · R23–R27 · R34 · R36 · R37 · R46 · R50 · R58–R63 을 수정 없이 재사용한다. NSO_PINNED · NCA_PINNED_MANIFEST 는 바꾸지 않는다.
#
# 실행 (W 에서, dcta 환경만)
#   run-r dcta R/branch_biology/64_run_atlas_rna.R --data-root <M> --o-root <O> --baseline-dir <N>/results/neut_stage_observation_2026-10-19/attempt_02 \
#     --baseline-connectivity-dir <N>/results/neut_stage_observation_2026-10-19/connectivity_audit_02 --out-dir <W>/results/neutrophil_atlas_rna_2026-10-19/attempt_01 \
#     [--source-dir <W>/data/external/granulopoiesis_atlas_2026-10-19] [--source-manifest <W>/docs/experiments/neutrophil_atlas_rna_2026-10-19/source_downloads.json] [--chunk 5000]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(igraph) })

ATLAS_RUN_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                       "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                       "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                       "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R")
ATLAS_SOURCE_DIR <- "data/external/granulopoiesis_atlas_2026-10-19"                                    # W 기준 (읽기 전용)
ATLAS_SOURCE_MANIFEST <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/source_downloads.json"      # Codex 소유 다운로드 기록 (읽기 전용)
ATLAS_PINNED <- list(source = c(tpm = "512ad89a2b60a6a535da9d830cb4741a73bc2d425c66d1ef659cd25c4e65631d", counts = "5a48bc2219e18d2ce1e640842b61ccc9154edb5afd84d1498f405709410cff8f",
                                soft = "1ab75c576c6f6baebccc64272291f512e4ef6472621f1bcb48e1050ee55d5e86", hgnc = "12363916be36f8f6f91879b0faade4eeaac973722de7136ec61329b77f2c2f91"),
                     h5ad = "4149b357b7d0027027988b8c78a92ae8894f36f3eae667762500dbc976f3c454", graph = "e63392de28aee7e4b629221bc8a0ba6118e5048c27b55071019b1d247ffa62de",
                     baseline = "4c31d681b0da2fc11b747dd865ce68431a10d161896f03515557b2d54d3c36a7",               # NSO attempt_02 manifest.json (지시서)
                     baseline_connectivity = "aae82d2eddb4b9faf761605de65e2474296cef718f1ce081b88e0d1bdee01fad",  # connectivity_audit_02 manifest.json (지시서)
                     n_cells = 263159L)
ATLAS_FULL_NAMES <- c("Myeloblast", "Promyelocyte", "Myelocyte", "Metamyelocyte", "Band neutrophil", "Segmented neutrophil")
ATLAS_NOT_ACCURACY <- "panel sizes differ from the baseline; more detected or connected cells is not higher stage accuracy and no stage label is made"


# ── 1. 조립 (순수) ───────────────────────────────────────────────────────────
#' source 파일 hash: docs 다운로드 기록(files = [{name, sha256}])과 같아야 하고, 고정값이 있으면 그것과도 같아야 한다
atlas_check_source_files <- function(source_dir, manifest, pinned = NULL) {
  if (!file.exists(manifest)) stop("source manifest 가 없다: ", manifest)
  m <- jsonlite::fromJSON(manifest)$files
  if (!is.data.frame(m) || !all(c("name", "sha256") %in% names(m))) stop("source manifest 형식이 다르다 (files = [{name, sha256}] 이어야 한다): ", manifest)
  p <- file.path(source_dir, unname(ATLAS_FILES)); miss <- unname(ATLAS_FILES)[!file.exists(p)]
  if (length(miss)) stop("source 파일이 없다: ", paste(miss, collapse = ", "))
  cur <- vapply(p, .nrv_sha, "", USE.NAMES = FALSE); rec <- m$sha256[match(unname(ATLAS_FILES), m$name)]
  if (anyNA(rec)) stop("source manifest 에 없는 source 파일: ", paste(unname(ATLAS_FILES)[is.na(rec)], collapse = ", "))
  bad <- unname(ATLAS_FILES)[cur != rec]; if (length(bad)) stop("source 파일 hash 가 source_downloads.json 과 다르다: ", paste(bad, collapse = ", "))
  if (!is.null(pinned)) { pv <- unname(pinned[names(ATLAS_FILES)])
    if (anyNA(pv) || any(cur != pv)) stop("source 파일 hash 가 고정값(ATLAS_PINNED$source)과 다르다: ", paste(unname(ATLAS_FILES)[is.na(pv) | cur != pv], collapse = ", ")) }
  data.frame(role = names(ATLAS_FILES), file = unname(ATLAS_FILES), path = p, sha256 = cur, manifest_sha256 = rec, bytes = as.numeric(file.size(p)), stringsAsFactors = FALSE)
}

#' CSV 저장 (design: source 표는 CSV). double 은 %.17g (mdu_write 와 같은 규칙), 글자 열만 따옴표. 다시 읽어 행 · 열 · 모든 칸 글자가 같은지 (TRUE / FALSE)
atlas_write_csv <- function(df, path) {
  df <- as.data.frame(df, stringsAsFactors = FALSE); qi <- which(vapply(df, function(x) is.character(x) || is.factor(x), logical(1)))
  for (j in seq_along(df)) { x <- df[[j]]
    df[[j]] <- if (is.double(x)) ifelse(is.na(x), "NA", sprintf("%.17g", x)) else if (is.logical(x)) ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))
               else if (is.integer(x)) ifelse(is.na(x), "NA", as.character(x))
               else { x <- as.character(x); if (anyNA(x)) stop("CSV 글자 열에 결측이 있다 (빈 글자로 쓰지 않고 멈춘다): ", names(df)[j], " (", basename(path), ")"); x } }
  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  tryCatch(utils::write.csv(df, con, row.names = FALSE, quote = if (length(qi)) qi else FALSE, na = "NA"), finally = close(con))
  back <- utils::read.csv(if (grepl("\\.gz$", path)) gzfile(path) else path, colClasses = "character", na.strings = character(0), check.names = FALSE, stringsAsFactors = FALSE)
  identical(dim(back), dim(df)) && identical(names(back), names(df)) && all(vapply(seq_along(df), function(j) identical(back[[j]], df[[j]]), logical(1)))
}

#' source 파라미터 · 규칙 표 (target 을 보기 전에 고정)
atlas_param_table <- function(ref) {
  p <- ref$params
  data.frame(name = c(names(p), "stages", "codes", "expression", "transition_rule", "panel_rule", "mapping_rule", "pb_pmn", "independence_note"),
             value = c(vapply(p, function(x) as.character(x), ""), paste(ref$stages, collapse = ">"), paste(ref$codes, collapse = "|"),
                       "median TPM per nominal code x stage (resequenced / sub-fraction columns of one code x stage pooled), then log2(TPM + 1)",
                       paste("adjacent pair (later minus earlier) over codes having both stages; candidate = n_pairs >= min_codes, high-side TPM >= min_tpm in >= min_codes_tpm codes,",
                             "|median| >= min_abs_median, sign agreement >= min_sign_agree; top_transition per sign by |median| (ties by ENSG); descriptive, not a DE test"),
                       paste("higher than every existing adjacent stage (same candidate rule on each contrast); end stages one neighbour; ranked by minimum adjacent effect,",
                             "top_panel (ties by ENSG); empty panel allowed (not filled); DCTA experimental rule, not the authors' official panel"),
                       "HGNC exact ensembl_gene_id one-to-one with locus_group protein-coding gene only; other candidates listed in source_exclusions.csv",
                       "peripheral-blood PMN columns kept separate (source_pb_pmn_*.csv); not used in contrasts or panels", ref$independence_note),
             stringsAsFactors = FALSE)
}

#' 혈액 PMN 별도 표: nominal code 묶음과 패널 유전자의 PMN log2(TPM + 1) (서술, 선택에 쓰지 않음)
atlas_pb_pmn_tables <- function(source, ref) {
  md <- source$metadata; k <- md$stage == "PB_PMN"
  gt <- data.frame(group = character(0), code = character(0), stage = character(0), n_columns = integer(0), columns = character(0), merge_notes = character(0),
                   subfractions = character(0), stringsAsFactors = FALSE)
  if (!any(k)) return(list(groups = gt, genes = data.frame(stage = character(0), rank = integer(0), ensg = character(0), symbol = character(0), stringsAsFactors = FALSE)))
  g <- atlas_code_stage_tpm(source$tpm[, k, drop = FALSE], md[k, , drop = FALSE])
  grp <- data.frame(g$groups, group_kind = rep("nominal_collection_group", nrow(g$groups)),
                    note = rep("peripheral-blood PMN; kept separate from the 6 marrow stages; not used in contrasts or panels", nrow(g$groups)), stringsAsFactors = FALSE)
  pn <- ref$panels; L <- log2(g$tpm_med[pn$ensg, , drop = FALSE] + 1)
  genes <- data.frame(stage = pn$stage, rank = pn$rank, ensg = pn$ensg, symbol = pn$symbol, stringsAsFactors = FALSE)
  genes$pb_pmn_log2_median <- if (nrow(pn)) apply(L, 1, stats::median) else numeric(0)
  for (cd in g$groups$code) genes[[paste0("pb_pmn_log2_", cd)]] <- unname(L[, paste(cd, "PB_PMN", sep = "|")])
  genes$note <- rep("descriptive only: blood PMN log2(TPM + 1) per nominal code; not used for any selection", nrow(genes))
  list(groups = grp, genes = genes)
}

#' target 열 이름: 패널 · 변화 ENSG 와 정확히 같은 열에는 HGNC symbol, 나머지 열은 ENSG 그대로 (symbol 이름 차이로 생기는 오대응 방지)
atlas_target_labels <- function(target_ensg, genes) {
  lab <- as.character(target_ensg); j <- match(lab, genes$ensg); lab[!is.na(j)] <- genes$symbol[j[!is.na(j)]]; lab
}

#' NSO(R58–R59) 호환 reference: 단계 marker = atlas 단계 양성 패널(공유 표시), 단계 사이 변화 = atlas 선택 변화(evidence "atlas"). 빈 패널 단계는 marker 0 개.
atlas_nso_reference <- function(ref) {
  pn <- ref$panels; tr <- ref$transitions
  mk <- data.frame(stage_id = pn$stage, stage_rank = match(pn$stage, NSO_STAGES), gene = pn$symbol, gene_order = as.integer(pn$rank), ensg = pn$ensg, stringsAsFactors = FALSE)
  if (anyNA(mk$stage_rank)) stop("패널 단계가 형태 6 단계가 아니다: ", paste(unique(mk$stage_id[is.na(mk$stage_rank)]), collapse = ", "))
  mk <- mk[order(mk$stage_rank, mk$gene_order, method = "radix"), , drop = FALSE]; rownames(mk) <- NULL
  mk$shared_with <- vapply(seq_len(nrow(mk)), function(i) paste(setdiff(mk$stage_id[mk$gene == mk$gene[i]], mk$stage_id[i]), collapse = "|"), "")
  mk$marker_kind <- rep("atlas_rna_stage_positive", nrow(mk))
  stages <- data.frame(stage_id = NSO_STAGES, stage_rank = seq_along(NSO_STAGES), source_stage_order = seq_along(NSO_STAGES), full_name = ATLAS_FULL_NAMES,
                       n_markers = as.integer(table(factor(mk$stage_id, levels = NSO_STAGES))),
                       marker_genes = vapply(NSO_STAGES, function(s) paste(mk$gene[mk$stage_id == s], collapse = "|"), ""), row.names = NULL, stringsAsFactors = FALSE)
  sp <- strsplit(tr$pair, ">", fixed = TRUE)
  ch <- data.frame(edge_id = tr$pair, from = vapply(sp, `[`, "", 1L), to = vapply(sp, `[`, "", 2L), gene = tr$symbol, direction = tr$direction,
                   evidence = rep("atlas", nrow(tr)), gene_order = as.integer(tr$rank), ensg = tr$ensg, stringsAsFactors = FALSE)
  bad <- ch$edge_id[match(ch$to, NSO_STAGES) != match(ch$from, NSO_STAGES) + 1L]
  if (length(bad)) stop("atlas 변화가 바로 다음 단계 쌍이 아니다: ", paste(unique(bad), collapse = ", "))
  ug <- unique(mk$gene); ns <- vapply(ug, function(g) sum(mk$gene == g), 0L)
  sh <- do.call(rbind, lapply(ug[ns >= 2L], function(g) { r <- mk$stage_rank[mk$gene == g]
    data.frame(gene = g, stages = paste(mk$stage_id[mk$gene == g], collapse = "|"), n_stages = length(r), adjacent_only = all(diff(r) == 1L), stringsAsFactors = FALSE) }))
  if (is.null(sh)) sh <- data.frame(gene = character(0), stages = character(0), n_stages = integer(0), adjacent_only = logical(0), stringsAsFactors = FALSE)
  list(stages = stages, stage_markers = mk, changes = ch, shared_markers = sh)
}

#' 단계 가용성: 빈 패널 · target 에 유전자 없음 · 검출 기준(min_markers)보다 적은 관측 유전자는 판정불가 (0 이 아니다)
atlas_stage_status <- function(n_panel, n_observed, min_markers) {
  ifelse(n_panel == 0L, "unavailable_empty_panel", ifelse(n_observed == 0L, "unavailable_missing_genes",
         ifelse(n_observed < min_markers, "unavailable_too_few_genes", "available")))
}

#' R62 순서 도달 표에 판정불가 단계를 표시한다: 그 단계 수는 결측, 첫 끊김이 판정불가 단계면 break_kind = unavailable_stage 이고 그 뒤 단계의 도달은 평가하지 않는다.
#' 실제 edge 부재로 먼저 끊겼으면(no_real_edge) 뒤 단계의 0 은 그대로 둔다. 첫 단계가 판정불가면 first_break = "start>Myeloblast".
atlas_mark_unavailable <- function(r, unavailable, reason = "unavailable_empty_panel") {
  t <- r$table; st <- t$stage; un <- st %in% unavailable; rs <- rep_len(reason, length(unavailable))
  t$stage_status <- ifelse(un, rs[match(st, unavailable)], "available")
  cnt <- c("n_cells", "n_components"); rch <- c("n_reachable_cells", "n_reachable_components", "n_reappearing_cells"); fb <- ""; kind <- ""; cut <- FALSE
  for (k in seq_along(st)) {
    if (un[k]) { t[k, c(cnt, rch)] <- NA; t$reached[k] <- NA
      if (!nzchar(fb)) { fb <- if (k == 1L) paste0("start>", st[1]) else paste(st[k - 1L], st[k], sep = ">"); kind <- "unavailable_stage"; cut <- TRUE }
    } else if (cut) { t[k, rch] <- NA; t$reached[k] <- NA; t$stage_status[k] <- "not_evaluable_after_unavailable"
    } else if (!nzchar(fb) && k > 1L && !isTRUE(t$reached[k])) { fb <- paste(st[k - 1L], st[k], sep = ">"); kind <- "no_real_edge" }
  }
  list(table = t, first_break = fb, break_kind = kind)
}

#' 패널 · 변화 유전자의 target 존재 · 전체 세포 검출 (단계 관측 가능성을 확정하는 표가 아니다)
atlas_availability <- function(ref, ge, gm, target_symbols) {
  n <- ge$n_cells; have <- colnames(ge$cell$log)
  nd <- stats::setNames(as.numeric(Matrix::colSums(ge$cell$det)), have); ml <- stats::setNames(as.numeric(Matrix::colMeans(ge$cell$log)), have)
  part <- function(kind, d, where, dir) if (!nrow(d)) NULL else data.frame(item_kind = kind, stage_or_pair = where, direction = dir, rank = as.integer(d$rank), ensg = d$ensg,
                                                                      symbol = d$symbol, stringsAsFactors = FALSE)
  d <- rbind(part("stage_panel", ref$panels, ref$panels$stage, rep("positive", nrow(ref$panels))), part("transition", ref$transitions, ref$transitions$pair, ref$transitions$direction))
  if (is.null(d)) stop("패널 · 변화 유전자가 하나도 없다 (target 에 적용할 것이 없다).")
  ii <- d$symbol %in% have
  d$in_input <- ii; d$target_columns <- vapply(d$symbol, function(g) paste(gm[[g]], collapse = "|"), "", USE.NAMES = FALSE)
  d$target_symbol <- vapply(d$symbol, function(g) paste(unique(as.character(target_symbols[gm[[g]]])), collapse = "|"), "", USE.NAMES = FALSE)
  d$n_cells_detected <- ifelse(ii, unname(nd[d$symbol]), NA_real_); d$det_rate_all <- d$n_cells_detected / n; d$mean_log_all <- ifelse(ii, unname(ml[d$symbol]), NA_real_)
  d$availability <- ifelse(!ii, "missing_gene", ifelse(d$n_cells_detected == 0, "observed_zero", "observed"))
  d$interpretation <- rep("exact-ENSG presence and detection over all cells of this dataset; missing_gene is not biological absence; not a stage call", nrow(d))
  d
}

#' 세포 표: 전체 원 cell_id 순서, atlas 단계 패널 6 개의 검출 유전자 수 · logRNA 평균(판정불가는 결측), unavailable 이유. 예측 단계 없음.
atlas_cell_table <- function(md, ge, pe) {
  cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  d <- data.frame(cell_id = md$cell_id, cluster92 = md$cluster92, cluster51 = md$cluster51, block = md$block, study = md$study, library_size = ge$cell$library,
                  stringsAsFactors = FALSE)
  for (k in paste0("marker_", NSO_STAGES)) { d[[paste0(k, "_n_detected")]] <- cn[, k]; d[[paste0(k, "_mean_log")]] <- cm[, k] }
  d$unavailable_reason <- pe$cells$unavailable
  d
}

#' R62 연결 규칙 그대로: 검출 집합(관측 패널 유전자 중 log > 0 인 것이 기준 이상) → 단계 유도 성분 · 인접 쌍 원 edge · 순서 도달 (기준 x 그래프).
#' 판정불가 단계는 결측 행으로 남긴다. markers = 단계 → target 에 있는 패널 유전자, n_panel = 단계 → source 패널 유전자 수
atlas_connectivity <- function(n, e, block, L, markers, n_panel, criteria = NCA_CRITERIA, graphs = NCA_GRAPHS) {
  st <- NSO_STAGES; markers <- lapply(markers[st], as.character); names(markers) <- st
  status <- lapply(criteria, function(k) stats::setNames(atlas_stage_status(as.integer(n_panel[st]), lengths(markers), k), st))
  sets <- lapply(names(criteria), function(cr) { m <- markers; m[status[[cr]] != "available"] <- list(character(0)); nca_detection_sets(L, m, criteria[[cr]]) })
  names(sets) <- names(criteria)
  sc <- list(); pe <- list(); rt <- list()
  for (cr in names(criteria)) for (gr in names(graphs)) { S <- sets[[cr]]; sb <- graphs[[gr]]; s0 <- status[[cr]]
    for (k in seq_along(st)) { z <- if (s0[k] == "available") nca_components(n, e, keep = S[, k], block = block, same_block = sb) else NULL
      sc[[length(sc) + 1L]] <- data.frame(criterion = cr, graph = gr, stage = st[k], stage_status = s0[[k]], n_cells = if (is.null(z)) NA_integer_ else z$n_cells,
                                          n_within_edges = if (is.null(z)) NA_integer_ else z$n_edges, n_components = if (is.null(z)) NA_integer_ else z$n_components,
                                          largest_component = if (is.null(z)) NA_integer_ else z$largest, n_isolated = if (is.null(z)) NA_integer_ else z$n_isolated,
                                          stringsAsFactors = FALSE) }
    for (k in seq_len(length(st) - 1L)) { ok <- s0[k] == "available" && s0[k + 1L] == "available"
      pz <- nca_pair_edges(e, S[, k], S[, k + 1L], block, same_block = sb); if (!ok) pz[1, ] <- NA
      pe[[length(pe) + 1L]] <- data.frame(criterion = cr, graph = gr, from = st[k], to = st[k + 1L], pair_status = if (ok) "available" else "unavailable_stage", pz,
                                          stringsAsFactors = FALSE) }
    r <- atlas_mark_unavailable(nca_ordered_reach(n, e, S, block = block, same_block = sb), unavailable = st[s0 != "available"], reason = unname(s0[s0 != "available"]))
    rt[[length(rt) + 1L]] <- data.frame(criterion = cr, graph = gr, r$table, first_break = r$first_break, break_kind = r$break_kind, stringsAsFactors = FALSE) }
  bind <- function(x) { d <- do.call(rbind, x); rownames(d) <- NULL; d }
  det <- data.frame(stage = st, n_panel_genes = as.integer(n_panel[st]), n_observed_genes = lengths(markers), observed_genes = vapply(markers, paste, "", collapse = "|"),
                    stringsAsFactors = FALSE, row.names = NULL)
  for (cr in names(criteria)) { det[[paste0("status_", cr)]] <- unname(status[[cr]])
    det[[paste0("n_cells_", cr)]] <- ifelse(status[[cr]] == "available", as.numeric(colSums(sets[[cr]])), NA_real_) }
  list(status = status, sets = sets, detection = det, stage_components = bind(sc), pair_edges = bind(pe), reach = bind(rt))
}

#' 순서 도달 표의 내부 일관성: 도달 ≤ 검출 (결측 제외), 첫 미도달 단계 뒤에는 도달 단계가 없다 (결측 제외)
atlas_reach_consistent <- function(rt) {
  ok1 <- all(is.na(rt$n_reachable_cells) | is.na(rt$n_cells) | rt$n_reachable_cells <= rt$n_cells)
  ok2 <- all(vapply(split(rt, paste(rt$criterion, rt$graph)), function(z) { r <- z$reached; b <- which(!is.na(r) & !r)
    !length(b) || !any(r[min(b):length(r)] %in% TRUE) }, TRUE))
  ok1 && ok2
}

#' 기준선 비교 (단계): 패널 크기 · 유전자 · 공유, ≥1/≥2 검출 세포 수 · 겹침, ge1 순서 도달 (전체 · 같은 block). 판정불가는 결측.
atlas_baseline_comparison <- function(base_mk, base_nd, base_reach, ref, con) {
  st <- NSO_STAGES; det <- con$detection
  rows <- lapply(seq_along(st), function(k) { s <- st[k]; bg <- as.character(base_mk$gene[base_mk$stage_id == s]); ag <- ref$panels$symbol[ref$panels$stage == s]
    d <- data.frame(stage = s, baseline_panel_size = length(bg), baseline_panel_genes = paste(bg, collapse = "|"), atlas_panel_size = length(ag),
                    atlas_panel_status = ref$panel_status$status[ref$panel_status$stage == s], atlas_panel_genes = paste(ag, collapse = "|"),
                    atlas_panel_in_target = det$n_observed_genes[k], shared_panel_genes = paste(intersect(bg, ag), collapse = "|"), stringsAsFactors = FALSE)
    for (cr in names(NCA_CRITERIA)) { v <- base_nd[, paste0("marker_", s)]; b <- !is.na(v) & v >= NCA_CRITERIA[[cr]]; av <- con$status[[cr]][[s]] == "available"; A <- con$sets[[cr]][, s]
      d[[paste0("baseline_n_", cr)]] <- sum(b); d[[paste0("atlas_n_", cr)]] <- if (av) sum(A) else NA_integer_
      d[[paste0("overlap_", cr)]] <- if (av) sum(A & b) else NA_integer_; d[[paste0("baseline_only_", cr)]] <- if (av) sum(b & !A) else NA_integer_
      d[[paste0("atlas_only_", cr)]] <- if (av) sum(A & !b) else NA_integer_; d[[paste0("atlas_status_", cr)]] <- con$status[[cr]][[s]] }
    for (gr in names(NCA_GRAPHS)) { zb <- base_reach[base_reach$criterion == "ge1" & base_reach$graph == gr & base_reach$stage == s, , drop = FALSE]
      za <- con$reach[con$reach$criterion == "ge1" & con$reach$graph == gr & con$reach$stage == s, , drop = FALSE]
      if (nrow(zb) != 1L || nrow(za) != 1L) stop("순서 도달 표에서 단계 행을 하나로 찾지 못했다: ", s, " ", gr)
      d[[paste0("baseline_reach_", gr)]] <- as.numeric(zb$n_reachable_cells); d[[paste0("atlas_reach_", gr)]] <- za$n_reachable_cells }
    d$note <- ATLAS_NOT_ACCURACY; d })
  do.call(rbind, rows)
}

#' 기준선 비교 (순서 도달 전체): 기준 x 그래프 x 단계의 검출 · 도달 · 첫 끊김
atlas_reach_comparison <- function(base_reach, con) {
  a <- con$reach; key <- function(d) paste(d$criterion, d$graph, d$stage); j <- match(key(a), key(base_reach))
  if (anyNA(j)) stop("기준선 순서 도달 표에 없는 기준 x 그래프 x 단계: ", paste(utils::head(key(a)[is.na(j)], 3), collapse = ", "))
  b <- base_reach[j, , drop = FALSE]
  data.frame(criterion = a$criterion, graph = a$graph, stage = a$stage, baseline_n_cells = as.numeric(b$n_cells), baseline_n_reachable_cells = as.numeric(b$n_reachable_cells),
             baseline_first_break = as.character(b$first_break), atlas_n_cells = a$n_cells, atlas_n_reachable_cells = a$n_reachable_cells, atlas_stage_status = a$stage_status,
             atlas_first_break = a$first_break, atlas_break_kind = a$break_kind, note = rep(ATLAS_NOT_ACCURACY, nrow(a)), stringsAsFactors = FALSE)
}

#' 기준선 비교 (p92 · p51 군집): 군집 안 ≥1/≥2 검출 세포 비율, 기준선 vs atlas (판정불가는 결측). 정확도 아님.
atlas_cluster_comparison <- function(md, base_nd, con) {
  out <- list()
  for (p in c("p92", "p51")) { cl <- if (p == "p92") md$cluster92 else md$cluster51; uc <- fdt_cluster_order_safe(cl); f <- factor(cl, levels = uc); nc <- as.integer(table(f))
    for (s in NSO_STAGES) { d <- data.frame(partition = p, cluster = uc, n_cells = nc, stage = s, stringsAsFactors = FALSE)
      for (cr in names(NCA_CRITERIA)) { v <- base_nd[, paste0("marker_", s)]; b <- !is.na(v) & v >= NCA_CRITERIA[[cr]]; av <- con$status[[cr]][[s]] == "available"
        d[[paste0("baseline_frac_", cr)]] <- as.numeric(tapply(b, f, mean))
        d[[paste0("atlas_frac_", cr)]] <- if (av) as.numeric(tapply(con$sets[[cr]][, s], f, mean)) else rep(NA_real_, length(uc))
        d[[paste0("atlas_status_", cr)]] <- rep(con$status[[cr]][[s]], length(uc)) }
      out[[length(out) + 1L]] <- d } }
  d <- do.call(rbind, out); rownames(d) <- NULL; d$note <- rep("descriptive cluster summary; not accuracy", nrow(d)); d
}

#' 실제 h5ad reader: R27 nrv_h5_reader + raw/var 색인(ENSG)
atlas_h5_reader <- function(data_root) {
  r <- nrv_h5_reader(data_root); h5 <- file.path(data_root, NRV_INPUTS[["h5ad"]])
  r$var_ids <- function() { f <- hdf5r::H5File$new(h5, "r"); on.exit(f$close_all()); g <- f[["raw/var"]]; x <- g[[as.character(hdf5r::h5attr(g, "_index"))]]
    if (inherits(x, "H5Group")) stop("raw/var 색인이 범주형이다 (ENSG 문자열 열이 아니다).")
    as.character(x$read()) }
  r
}


# ── 2. 실행기 ────────────────────────────────────────────────────────────────
#' 인자 해석 → 결과 폴더 검사(6) → 소스 스냅샷 → 사본 적재 → 사본의 atlas_build. root / reader / pins 는 테스트용.
atlas_main <- function(argv, root = NULL, reader = NULL, pins = ATLAS_PINNED) {
  a <- list(source_dir = NULL, source_manifest = NULL, data_root = NULL, o_root = NULL, base = NULL, base_con = NULL, out = NULL, chunk = "5000")
  keys <- c("--source-dir" = "source_dir", "--source-manifest" = "source_manifest", "--data-root" = "data_root", "--o-root" = "o_root", "--baseline-dir" = "base",
            "--baseline-connectivity-dir" = "base_con", "--out-dir" = "out", "--chunk" = "chunk")
  i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  a$chunk <- suppressWarnings(as.integer(a$chunk))
  bad <- vapply(a[c("data_root", "o_root", "base", "base_con", "out")], function(x) is.null(x) || is.na(x) || !nzchar(x), TRUE)
  if (any(bad) || is.na(a$chunk) || a$chunk < 1L) { cat("--data-root, --o-root, --baseline-dir, --baseline-connectivity-dir, --out-dir 와 양의 정수 --chunk 가 필요하다.\n"); return(2L) }
  out <- a$out
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", out, "\n"); return(6L) }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  if (is.null(a$source_dir)) a$source_dir <- file.path(root, ATLAS_SOURCE_DIR)
  if (is.null(a$source_manifest)) a$source_manifest <- file.path(root, ATLAS_SOURCE_MANIFEST)
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", ATLAS_RUN_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(ATLAS_RUN_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = ATLAS_RUN_SOURCES)
    get("atlas_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, reader = reader, pins = pins)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

atlas_build <- function(a, out, snap, log, argv, reader, pins) {
  ck <- .mbvq_checks(); t_start <- Sys.time()
  rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  # ── (1) source 전용 고정: target 은 아직 열지 않는다 ──
  sf <- atlas_check_source_files(a$source_dir, a$source_manifest, pins$source)
  ck$add("source_files_pinned", TRUE, sprintf("%d source files sha256 == source_downloads.json%s", nrow(sf), if (is.null(pins$source)) " (code pins disabled)" else " and ATLAS_PINNED$source"))
  s <- atlas_read_source(a$source_dir); ref <- atlas_build_reference(s); pmn <- atlas_pb_pmn_tables(s, ref)
  log("source reference built:", nrow(s$tpm), "genes,", ncol(s$tpm), "columns;", nrow(ref$panels), "panel genes;", nrow(ref$transitions), "transition genes")
  sdir <- file.path(out, "source"); dir.create(sdir, showWarnings = FALSE); sfiles <- character(0); crt <- logical(0)
  wc <- function(df, f) { p0 <- file.path(sdir, f); crt[[f]] <<- atlas_write_csv(df, p0); sfiles <<- c(sfiles, p0); invisible(p0) }
  md0 <- s$metadata; md0$collection_group <- paste(md0$code, md0$stage, sep = "|"); md0$group_kind <- rep("nominal_collection_group", nrow(md0))
  hs <- do.call(rbind, lapply(sort(unique(s$map$map_status)), function(m) data.frame(map_status = m, n_genes = sum(s$map$map_status == m),
                                                                                         n_use_for_panels = sum(s$map$use_for_panels & s$map$map_status == m), stringsAsFactors = FALSE)))
  wc(sf, "source_files.csv"); wc(atlas_param_table(ref), "source_params.csv"); wc(md0, "source_metadata.csv")
  wc(data.frame(ref$groups, group_kind = rep("nominal_collection_group", nrow(ref$groups)), stringsAsFactors = FALSE), "source_groups.csv")
  wc(ref$all_contrasts, "source_contrasts_all.csv.gz"); wc(ref$transitions, "source_transitions.csv"); wc(ref$panels, "source_panels.csv")
  wc(ref$panel_status, "source_panel_status.csv"); wc(ref$exclusions, "source_exclusions.csv"); wc(s$map, "source_hgnc_map.csv.gz"); wc(hs, "source_hgnc_summary.csv")
  wc(pmn$groups, "source_pb_pmn_groups.csv"); wc(pmn$genes, "source_pb_pmn_panel_genes.csv")
  p0 <- file.path(sdir, "source_reference.rds")
  saveRDS(list(contract = "neutrophil_atlas_rna_2026-10-19/codex_instruction.md", reference = ref, map = s$map, metadata = s$metadata, files = sf, pb_pmn = pmn), p0, version = 3)
  sfiles <- c(sfiles, p0)
  if (!all(crt)) stop("source CSV 왕복 확인 실패: ", paste(names(crt)[!crt], collapse = ", "))
  ck$add("source_tables_roundtrip_ok", TRUE, sprintf("%d source CSVs re-read with identical dimensions, names and every cell text", length(crt)))
  sfz <- file.path(out, "source_freeze.json")
  writeLines(jsonlite::toJSON(list(phase = "source-only reference fixed before any target RNA, graph, cell table, UMAP or author label was opened", written_at = stamp(),
                                   contract = "results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md", params = ref$params,
                                   source_files = as.list(stats::setNames(sf$sha256, sf$file)), files = as.list(stats::setNames(vapply(sfiles, .nrv_sha, ""), rel(sfiles)))),
                              auto_unbox = TRUE, pretty = TRUE, digits = NA), sfz)
  src_hash <- vapply(sfiles, .nrv_sha, ""); hz <- .nrv_sha(sfz)
  log("source frozen:", length(sfiles), "files; source_freeze.json", substr(hz, 1, 12)); rm(s)
  # ── (2) target 입력 고정 (source 고정 뒤) ──
  if (is.null(reader)) reader <- atlas_h5_reader(a$data_root)
  bm <- c(baseline = file.path(a$base, "manifest.json"), baseline_connectivity = file.path(a$base_con, "manifest.json"))
  if (!all(file.exists(bm))) stop("기준선 manifest 가 없다: ", paste(bm[!file.exists(bm)], collapse = ", "))
  hb <- vapply(bm, .nrv_sha, "")
  for (k in names(bm)) if (!is.null(pins[[k]]) && !identical(unname(hb[[k]]), pins[[k]])) stop(k, " manifest hash 가 지시서 고정값과 다르다: ", hb[[k]])
  st <- lapply(NSO_NRV, function(p) { att <- file.path(a$o_root, p)
    c(cells = nrv_stage_file(att, "structure", "cells_lossless.rds"), edges = nrv_stage_file(att, "structure", "edges.rds"), input = nrv_stage_file(att, "structure", "input_sha256.tsv")) })
  fk <- c("graph", "graph_ids", "h5ad"); f <- reader$input_files()[fk]
  if (anyNA(f) || !all(file.exists(f))) stop("입력 파일이 없다: ", paste(fk[is.na(f) | !file.exists(f)], collapse = ", "))
  cur <- vapply(f, .nrv_sha, "")
  for (p in names(st)) { rec <- utils::read.delim(st[[p]][["input"]], stringsAsFactors = FALSE, colClasses = "character"); r <- rec$sha256[match(fk, rec$input)]
    if (anyNA(r) || !identical(unname(cur), r)) stop("현재 입력 hash 가 NRV structure 단계 기록과 다르다 (", p, "): ", paste(fk[is.na(r) | r != cur], collapse = ", ")) }
  if (!is.null(pins$graph) && !identical(unname(cur[["graph"]]), pins$graph)) stop("graph_base.npz hash 가 계약값과 다르다: ", cur[["graph"]])
  if (!is.null(pins$h5ad) && !identical(unname(cur[["h5ad"]]), pins$h5ad)) stop("h5ad hash 가 계약값과 다르다: ", cur[["h5ad"]])
  ck$add("inputs_pinned", TRUE, sprintf("h5ad, graph_base.npz == contract%s; graph, graph_ids, h5ad == NRV structure input_sha256 (p92, p51); baseline manifests %s / %s == instruction%s",
                                        if (is.null(pins$h5ad)) " (pins disabled)" else "", substr(hb[["baseline"]], 1, 12), substr(hb[["baseline_connectivity"]], 1, 12),
                                        if (is.null(pins$baseline)) " (pins disabled)" else ""))
  log("inputs pinned")
  # ── 세포 · 계산용 메타데이터 · 원 edge ──
  nref <- atlas_nso_reference(ref)
  c92 <- readRDS(st$p92[["cells"]])$cells; c51 <- readRDS(st$p51[["cells"]])$cells
  chk <- nso_check_inputs(c92, c51, raw_ids = reader$obs("cell_id"), graph_ids = reader$ids(), reference = nref)
  for (r in seq_len(nrow(chk$checks))) ck$add(chk$checks$check[r], chk$checks$pass[r], chk$checks$detail[r])
  md <- nso_calculation_metadata(data.frame(chk$cells, donor = as.character(reader$obs("donor_id")), stringsAsFactors = FALSE))
  if (!identical(as.character(reader$obs("Study")), md$study)) stop("NRV 세포표 study 가 h5 obs/Study 와 다르다.")
  ck$add("metadata_matches_h5", TRUE, "block == R24 nrv_blocks(obs/donor_id) and study == obs/Study for every cell")
  n <- nrow(md); rm(c92, c51)
  if (!is.null(pins$n_cells)) { ok <- n == pins$n_cells; ck$add("denominator_263159", ok, sprintf("%d cells == pinned %d (all cells kept, original order)", n, pins$n_cells))
    if (!ok) stop("세포 분모가 고정값과 다르다: ", n) }
  he <- vapply(st, function(z) .nrv_sha(z[["edges"]]), ""); e <- readRDS(st$p92[["edges"]])
  if (he[["p92"]] != he[["p51"]] && !identical(e, readRDS(st$p51[["edges"]]))) stop("p92 와 p51 의 edges.rds 가 다르다 (같은 그래프여야 한다).")
  if (nrow(e) && any(e$i < 1L | e$j > n | e$i >= e$j)) stop("edge 가 1..n 상삼각 세포 번호가 아니다.")
  ck$add("edges_original_only", TRUE, sprintf("%d positive upper-triangle edges of the original graph (NRV edges.rds, identical for p92 and p51); no edge added, removed or reweighted", nrow(e)))
  # ── 기준선 (NSO attempt · connectivity; 원저자 열 없음) ──
  bfile <- .mbvq_manifest_file(a$base, "manifest.json", "lossless.rds"); bl <- readRDS(bfile)
  base_mk <- bl$reference$stage_markers; base_nd <- bl$programs$cells$n_detected; base_ids <- as.character(bl$metadata$cell_id); rm(bl)
  okb <- identical(base_ids, md$cell_id) && is.matrix(base_nd) && nrow(base_nd) == n && all(paste0("marker_", NSO_STAGES) %in% colnames(base_nd))
  ck$add("baseline_cells_identical_order", okb, sprintf("baseline NSO lossless: %d cell_id identical to the current cell order; detection matrix %d x %d", length(base_ids), NROW(base_nd), NCOL(base_nd)))
  if (!okb) stop("기준선 NSO 세포 순서 · 검출 행렬이 현재 세포표와 다르다.")
  rfile <- .mbvq_manifest_file(a$base_con, "manifest.json", "ordered_reachability.tsv")
  base_reach <- utils::read.delim(rfile, stringsAsFactors = FALSE, colClasses = "character", quote = "", na.strings = character(0))
  if (!all(c("criterion", "graph", "stage", "n_cells", "n_reachable_cells", "first_break") %in% names(base_reach))) stop("기준선 ordered_reachability.tsv 에 필요한 열이 없다.")
  # ── (3) 원 RNA 근거 (raw 1 회, ENSG 정확 대응) ──
  tid <- as.character(reader$var_ids()); tsym <- as.character(reader$symbols())
  if (length(tid) != length(tsym)) stop("raw/var 색인(ENSG)과 gene_symbols 길이가 다르다: ", length(tid), " vs ", length(tsym))
  if (any(grepl("^ENSG[0-9]+\\.[0-9]+$", tid))) stop("raw/var 색인에 ENSG 버전 접미사가 있다 (정확 대응 전에 규칙을 정해야 한다).")
  genes <- unique(rbind(ref$panels[, c("ensg", "symbol")], ref$transitions[, c("ensg", "symbol")]))
  if (anyDuplicated(genes$ensg) || anyDuplicated(genes$symbol)) stop("패널 · 변화 유전자의 ENSG ↔ symbol 이 일대일이 아니다.")
  lab <- atlas_target_labels(tid, genes); gm <- nrv_gene_map(lab, genes$symbol)
  ck$add("target_genes_exact_ensg", TRUE, sprintf("%d panel/transition genes; %d found in raw/var index by exact ENSG (%d with more than one column, summed as R25); %d missing (kept as missing_gene, never 0); %d of %d target ids start with ENSG",
                                                  nrow(genes), sum(lengths(gm) > 0), sum(lengths(gm) > 1), sum(lengths(gm) == 0), sum(startsWith(tid, "ENSG")), length(tid)))
  groups <- list(p92 = data.frame(cluster = md$cluster92, block = md$block, stringsAsFactors = FALSE), p51 = data.frame(cluster = md$cluster51, block = md$block, stringsAsFactors = FALSE))
  ge <- nso_gene_evidence(.nrv_checked_chunks(reader$raw_chunks(a$chunk), n), gm, groups)
  log("RNA streamed:", n, "cells;", sum(ge$genes$in_input), "/", nrow(genes), "atlas genes in input")
  pe <- nso_program_evidence(ge, nref)
  den <- tapply(pe$groups$n_cells, paste(pe$groups$partition, pe$groups$program_id), sum)
  ck$add("denominators_preserved", all(den == n), sprintf("every partition x program sums cluster x block cells to %d (%d sums); cell tables keep all %d cells in original order", n, length(den), n))
  # ── 원 edge 연결 (R62 규칙) ──
  have <- colnames(ge$cell$log); pl <- split(ref$panels$symbol, factor(ref$panels$stage, levels = NSO_STAGES))
  con <- atlas_connectivity(n, e, md$block, ge$cell$log, lapply(pl, function(g) intersect(g, have)), lengths(pl))
  det <- con$detection; det$missing_genes <- vapply(NSO_STAGES, function(s0) paste(setdiff(pl[[s0]], have), collapse = "|"), "")
  det$panel_genes <- vapply(NSO_STAGES, function(s0) paste(pl[[s0]], collapse = "|"), ""); det$panel_status <- ref$panel_status$status[match(NSO_STAGES, ref$panel_status$stage)]
  okr <- atlas_reach_consistent(con$reach); ck$add("reach_tables_consistent", okr, "reachable <= detected where both are evaluable; after the first unreached stage no later stage is reached")
  for (cr in names(NCA_CRITERIA)) for (gr in names(NCA_GRAPHS)) { z <- con$reach[con$reach$criterion == cr & con$reach$graph == gr, ]
    log(cr, gr, "first break:", if (nzchar(z$first_break[1])) paste0(z$first_break[1], " (", z$break_kind[1], ")") else "none", "| reachable:", paste(z$n_reachable_cells, collapse = ",")) }
  rm(e)
  # ── 기준선 비교 ──
  cmp <- atlas_baseline_comparison(base_mk, base_nd, base_reach, ref, con); rcmp <- atlas_reach_comparison(base_reach, con); ccmp <- atlas_cluster_comparison(md, base_nd, con)
  # ── 계산 표 → calc_manifest.json (여기서 결과를 고정한다) ──
  files <- character(0); rtp <- logical(0)
  wt <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); rtp[[f]] <<- mbvt_write_tsv(df, p0); files <<- c(files, p0); invisible(p0) }
  wt(atlas_availability(ref, ge, gm, tsym), "target_gene_availability.tsv"); wt(atlas_cell_table(md, ge, pe), "cell_atlas_evidence.tsv.gz")
  wt(pe$groups, "cluster_block_programs.tsv.gz"); wt(ge$groups, "cluster_block_genes.tsv.gz"); wt(pe$clusters, "cluster_programs.tsv"); wt(pe$programs, "program_definitions.tsv")
  wt(det, "detection_sets.tsv"); wt(con$stage_components, "stage_components.tsv"); wt(con$pair_edges, "adjacent_pair_edges.tsv"); wt(con$reach, "ordered_reachability.tsv")
  wt(cmp, "baseline_comparison.tsv"); wt(rcmp, "reach_comparison.tsv"); wt(ccmp, "cluster_comparison.tsv")
  p0 <- file.path(out, "lossless.rds")
  ref_small <- ref; ref_small$all_contrasts <- NULL
  saveRDS(list(contract = "neutrophil_atlas_rna_2026-10-19/codex_instruction.md", source_freeze_sha256 = hz,
               settings = list(min_cells = NSO_MIN_CELLS, chunk = a$chunk, criteria = NCA_CRITERIA, graphs = NCA_GRAPHS, params = ref$params),
               reference = ref_small, nso_reference = nref, metadata = md, gene_map = gm, gene_evidence = ge, programs = pe, connectivity = con, detection = det,
               baseline = list(stage_markers = base_mk, n_detected = base_nd, reach = base_reach), comparison = list(stage = cmp, reach = rcmp, cluster = ccmp)), p0, version = 3)
  files <- c(files, p0)
  cmf <- file.path(out, "calc_manifest.json")
  writeLines(jsonlite::toJSON(list(phase = "calculation tables fixed before the original UMAP and author labels are read", written_at = stamp(), source_freeze_sha256 = hz,
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  calc_files <- files; log("calculation fixed:", length(files), "files")
  # ── (4) 표시 (고정 뒤): 원 UMAP ──
  U <- reader$umap()
  if (!is.matrix(U) || nrow(U) != n || ncol(U) < 2L || any(!is.finite(U[, 1:2]))) stop("UMAP 세포 수 · 모양이 세포 표와 다르다 (", if (is.matrix(U)) paste(dim(U), collapse = " x ") else class(U)[1], " vs ", n, ").")
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  f1 <- file.path(fig, "fig01_atlas_panel_umap.png"); f2 <- file.path(fig, "fig02_baseline_vs_atlas.png")
  t1 <- atlas_draw_panel_umap(f1, U[, 1:2], pe, nref, det, ko, fam); t2 <- atlas_draw_comparison(f2, cmp, rcmp, ko, fam)
  files <- c(files, f1, f2); wt(rbind(t1, t2), "figure_text.tsv"); log("figures written"); rm(U)
  # ── 사후 (고정 뒤): 원저자 라벨 ──
  au <- reader$author(); if (length(au) != n) stop("원저자 라벨 수가 세포 수와 다르다: ", length(au), " vs ", n)
  wt(nso_posthoc_composition(md, au), "posthoc_author_cluster_composition.tsv"); wt(nso_posthoc_programs(au, pe), "posthoc_author_program_summary.tsv")
  calc_now <- vapply(calc_files, .nrv_sha, ""); cm0 <- jsonlite::fromJSON(cmf)$files
  ck$add("author_and_umap_read_after_fix", identical(unname(calc_now), unname(unlist(cm0[rel(calc_files)]))),
         "UMAP and author labels were opened only after calc_manifest.json was written; fixed calculation files unchanged afterwards")
  sz0 <- jsonlite::fromJSON(sfz)$files
  ck$add("source_frozen_before_target", identical(unname(vapply(sfiles, .nrv_sha, "")), unname(src_hash)) && identical(unname(unlist(sz0[rel(sfiles)])), unname(src_hash)),
         "source_freeze.json was written before any target reader call; source tables unchanged at the end")
  # ── provenance ──
  allinp <- c(stats::setNames(sf$path, paste0("source_", sf$role)), source_manifest = a$source_manifest, f,
              unlist(lapply(names(st), function(p) stats::setNames(st[[p]], paste0("nrv_", p, "_", names(st[[p]]))))),
              baseline_manifest = bm[["baseline"]], baseline_lossless = bfile, baseline_connectivity_manifest = bm[["baseline_connectivity"]], baseline_reach = rfile)
  sh <- vapply(seq_along(allinp), function(k) { x <- allinp[[k]]; if (x %in% f) unname(cur[match(x, f)]) else if (x %in% sf$path) sf$sha256[match(x, sf$path)] else .nrv_sha(x) }, "")
  wt(data.frame(input = names(allinp), path = unname(allinp), sha256 = sh, stringsAsFactors = FALSE), "input_sha256.tsv")
  wt(snap, "source_snapshot.tsv")
  pk <- c("Matrix", "hdf5r", "igraph", "digest", "jsonlite")
  wt(data.frame(item = c("R_executable", "R", "platform", "libPaths", "packages", "title_font", "chunk", "min_cells", "criteria", "source_params", "random_seed", "command",
                         "working_dir", "started_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pk, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), if (ko) MBVR_TITLE_FAMILY else "default",
                          as.character(a$chunk), as.character(NSO_MIN_CELLS), paste(names(NCA_CRITERIA), NCA_CRITERIA, sep = "=", collapse = " "),
                          paste(names(ref$params), unlist(ref$params), sep = "=", collapse = " "), "none (no random numbers are drawn)",
                          paste(c("64_run_atlas_rna.R", argv), collapse = " "), getwd(), format(t_start, "%Y-%m-%dT%H:%M:%S%z"), stamp()), stringsAsFactors = FALSE), "environment.tsv")
  ck$add("tables_roundtrip_ok", all(rtp), sprintf("%d TSVs re-read with the same rows, field counts and text columns (mbvt_write_tsv)", length(rtp)))
  cks <- ck$table(); wt(cks, "checks.tsv")
  a1 <- con$reach[con$reach$criterion == "ge1" & con$reach$graph == "all_edges", ]; b1 <- con$reach[con$reach$criterion == "ge1" & con$reach$graph == "same_block", ]
  fmt_fb <- function(z) if (nzchar(z$first_break[1])) sprintf("%s (%s)", z$first_break[1], z$break_kind[1]) else "끊김 없음"
  rd <- c("# 호중구 atlas RNA 패널 적용 — 실행 결과", "",
          "- 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md · design.md · plan Task 2. 설명서: W 루트 README_neutrophil_atlas_rna_ko.md.",
          "- 순서: source/ 표 → source_freeze.json → target 적용 · 연결 · 기준선 비교 → calc_manifest.json → 그림 · posthoc_ 표 → manifest.json.",
          sprintf("- 세포 %s (전체 분모, 원 순서). atlas 유전자 %d / %d 가 target 에 있음 (ENSG 정확 대응).", .nca_fmt_n(n), sum(ge$genes$in_input), nrow(genes)),
          paste0("- atlas 패널 크기 (target 안): ", paste(sprintf("%s %d (%d)", NSO_STAGES, det$n_panel_genes, det$n_observed_genes), collapse = " · "), "."),
          paste0("- ≥1 검출 세포: ", paste(sprintf("%s %s", NSO_STAGES, ifelse(is.na(det$n_cells_ge1), "판정불가", .nca_fmt_n(det$n_cells_ge1))), collapse = " · "), "."),
          sprintf("- ge1 순서 도달 첫 끊김: 모든 원 edge %s · 같은 block edge %s.", fmt_fb(a1), fmt_fb(b1)),
          "- 기준선 패널과 크기가 달라 검출 · 연결 증가를 단계 정확도 향상으로 읽지 않는다. 예측 stage_label 은 만들지 않았다. code 는 nominal collection group 이다.",
          "- 그림: figures/fig01_atlas_panel_umap.png (atlas 6 단계 패널, 같은 원 UMAP), figures/fig02_baseline_vs_atlas.png (기준선 vs atlas 패널 크기 · 검출 · 순서 도달).")
  p0 <- file.path(out, "README_ko.md"); writeLines(rd, p0); files <- c(files, p0)
  allf <- c(sfiles, sfz, files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), source_freeze_sha256 = hz, files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .atl_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.atl_bb, .f))
  quit(status = atlas_main(commandArgs(trailingOnly = TRUE)))
}
