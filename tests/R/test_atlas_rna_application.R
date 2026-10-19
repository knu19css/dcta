# =============================================================================
# test_atlas_rna_application.R — atlas RNA 고정 패널의 target 적용 · 원 edge 연결 · 기준선 비교 (R64–65) 테스트
#   계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md, design.md, plan Task 2
#   missing vs observed zero, 빈 패널 = unavailable (0 과 구분), reference · hash 불일치 중지, 전체 분모 · 순서 보존, 원 edge 만,
#   기존 output 덮어쓰기 거절, source 고정 전 target 읽기 금지 · 계산 고정 전 UMAP · 저자 라벨 읽기 금지
# =============================================================================
library(testthat)
atl_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- atl_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "25_nrv_expression.R",
            "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R",
            "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R", "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R",
            "64_run_atlas_rna.R", "65_atlas_rna_view.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

AE <- sprintf("ENSG%011d", 1:8)

test_that("target 대응: 패널 ENSG 와 정확히 같은 target 열에만 symbol 을 붙이고, 나머지 열은 ENSG 그대로 둔다", {
  lab <- atlas_target_labels(c(AE[1], "ENSG99", AE[3]), data.frame(ensg = AE[c(1, 3, 8)], symbol = c("S1", "S3", "S8"), stringsAsFactors = FALSE))
  expect_equal(lab, c("S1", "ENSG99", "S3"))
})

test_that("NSO 호환 reference: 패널 → 단계 marker (공유 표시), 선택 변화 → atlas 변화, 빈 패널 단계는 marker 0 개로 남는다", {
  ref <- list(panels = data.frame(stage = c("Myeloblast", "MY", "MY", "MM"), rank = c(1L, 1L, 2L, 1L), ensg = AE[c(1, 3, 8, 3)], symbol = c("S1", "S3", "S8", "S3"), stringsAsFactors = FALSE),
              transitions = data.frame(pair = c("PM>MY", "PM>MY"), direction = c("up", "down"), rank = 1L, ensg = AE[c(3, 2)], symbol = c("S3", "S2"), stringsAsFactors = FALSE))
  nr <- atlas_nso_reference(ref)
  expect_equal(nr$stages$stage_id, NSO_STAGES); expect_equal(sum(nr$stage_markers$stage_id == "PM"), 0L)
  expect_equal(nr$stage_markers$shared_with[nr$stage_markers$stage_id == "MY" & nr$stage_markers$gene == "S3"], "MM")
  expect_equal(nr$shared_markers$stages[nr$shared_markers$gene == "S3"], "MY|MM")
  expect_equal(nr$changes$gene[nr$changes$direction == "down"], "S2"); expect_true(all(nr$changes$evidence == "atlas")); expect_equal(unique(nr$changes$edge_id), "PM>MY")
  expect_equal(unique(nr$changes$from), "PM"); expect_equal(unique(nr$changes$to), "MY")
})

test_that("단계 가용성: 빈 패널 · target 에 유전자 없음 · 기준보다 적은 관측 유전자는 모두 판정불가이고 이유가 다르다", {
  expect_equal(atlas_stage_status(n_panel = c(0L, 2L, 2L, 3L, 1L), n_observed = c(0L, 0L, 1L, 3L, 1L), min_markers = 2L),
               c("unavailable_empty_panel", "unavailable_missing_genes", "unavailable_too_few_genes", "available", "unavailable_too_few_genes"))
  expect_equal(atlas_stage_status(1L, 1L, 1L), "available")
})

test_that("연결 표: 판정불가 단계는 0 이 아니라 결측이고, 첫 끊김이 판정불가 단계인지 실제 edge 부재인지 구분하며, 판정불가 뒤 도달은 평가하지 않는다", {
  S <- cbind(Myeloblast = c(TRUE, FALSE, FALSE), PM = c(FALSE, TRUE, FALSE), MY = c(FALSE, FALSE, FALSE)); e <- data.frame(i = c(1L, 2L), j = c(2L, 3L))
  r <- atlas_mark_unavailable(nca_ordered_reach(3L, e, S), unavailable = "MY", reason = "unavailable_empty_panel")
  expect_true(is.na(r$table$n_cells[3])); expect_true(is.na(r$table$n_reachable_cells[3])); expect_equal(r$table$stage_status[3], "unavailable_empty_panel")
  expect_equal(r$first_break, "PM>MY"); expect_equal(r$break_kind, "unavailable_stage")
  r2 <- atlas_mark_unavailable(nca_ordered_reach(3L, e[1, ], cbind(S[, 1:2], MY = c(FALSE, FALSE, TRUE))), unavailable = character(0))
  expect_equal(r2$table$n_cells[3], 1L); expect_equal(r2$first_break, "PM>MY"); expect_equal(r2$break_kind, "no_real_edge"); expect_equal(r2$table$stage_status, rep("available", 3))
  S4 <- cbind(S, MM = c(FALSE, FALSE, TRUE)); r3 <- atlas_mark_unavailable(nca_ordered_reach(3L, e, S4), unavailable = "MY", reason = "unavailable_missing_genes")
  expect_equal(r3$table$n_cells[4], 1L); expect_true(is.na(r3$table$n_reachable_cells[4])); expect_equal(r3$table$stage_status[4], "not_evaluable_after_unavailable")
  expect_equal(r3$table$stage_status[3], "unavailable_missing_genes")
  r4 <- atlas_mark_unavailable(nca_ordered_reach(3L, e, S4), unavailable = "Myeloblast")
  expect_equal(r4$first_break, "start>Myeloblast"); expect_true(all(is.na(r4$table$n_reachable_cells)))
  r5 <- atlas_mark_unavailable(nca_ordered_reach(3L, e[2, ], S4), unavailable = "MY")                 # 끊김이 먼저(실제) → 뒤의 0 은 그대로 둔다
  expect_equal(r5$first_break, "Myeloblast>PM"); expect_equal(r5$break_kind, "no_real_edge"); expect_equal(r5$table$n_reachable_cells[4], 0L)
})


# ── 소형 합성 세계 (실행기 전체) ──────────────────────────────────────────────
atl_world <- function() {
  w <- tempfile("atl_world"); dir.create(w); out <- file.path(w, "out", "attempt_01"); sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  # source: 6 단계 x code A · B · C + PMN(A). TPM 15 = 그 단계만 높음. E7 = HGNC 없음, E8 = MY 높음이나 target 에 없음.
  src <- file.path(w, "src"); dir.create(src)
  ct <- c(MB = "Myeloblast", PM = "Promyelocyte", MC = "Myelocyte", MM = "Metamyelocyte", B = "Band", S = "Segmented")
  cols <- c(as.vector(t(outer(names(ct), c("a", "b", "c"), paste, sep = "_"))), "PMN_a")
  stg <- c(rep(names(ct), each = 3), "PMN")
  peak <- c(E1 = "MB", E2 = "PM", E3 = "MC", E4 = "MM", E5 = "S", E6 = "", E7 = "MB", E8 = "MC")
  tpm <- t(vapply(names(peak), function(g) if (g == "E6") rep(3, length(cols)) else ifelse(stg == peak[[g]], 15, 0), numeric(length(cols))))
  tp <- data.frame(tpm, ENSG = AE, check.names = FALSE); names(tp)[seq_along(cols)] <- cols
  for (f in c("GSE294330_exp_mrn_tpm.csv.gz", "GSE294330_exp_mrn_counts.csv.gz")) { con <- gzfile(file.path(src, f), "w"); utils::write.csv(tp, con, row.names = FALSE); close(con) }
  tt <- function(k) paste0("mRNA.cellstage_", if (stg[k] == "PMN") "Neutrophil" else ct[[stg[k]]], "_donor_", toupper(sub("^.*_", "", cols[k])))
  soft <- c("^SERIES = GSE0", unlist(lapply(seq_along(cols), function(k) c(paste0("^SAMPLE = GSM", k), paste0("!Sample_title = ", tt(k)),
    paste0("!Sample_characteristics_ch1 = cell type: ", if (stg[k] == "PMN") "Neutrophil" else ct[[stg[k]]]), paste0("!Sample_description = Column name in raw_counts_mrn.csv: ", cols[k])))))
  writeLines(soft, file.path(src, "GSE294330_full.soft"))
  utils::write.table(data.frame(hgnc_id = paste0("HGNC:", c(1:6, 8)), symbol = paste0("S", c(1:6, 8)), locus_group = "protein-coding gene", ensembl_gene_id = AE[c(1:6, 8)], status = "Approved"),
                     file.path(src, "hgnc_complete_set.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  smf <- file.path(w, "source_downloads.json")                                               # 실제 docs 기록과 같은 형식: files = [{name, sha256}, ...]
  writeLines(jsonlite::toJSON(list(files = lapply(unname(ATLAS_FILES), function(f) list(name = f, sha256 = sha(file.path(src, f))))), auto_unbox = TRUE, pretty = TRUE), smf)
  # target 6 세포: 1–4 가 MB · PM · MY · MM 유전자 각 5 count, 5 = S6 만, 6 = 아무것도. S5(SN)는 모든 세포 0, E8 은 target 에 없음.
  ids <- paste0("cell", 1:6); donor <- c("d1", "d1", "d1", "d2", "d2", "d2"); tid <- c(AE[1:6], "ENSG99999999999"); tsym <- c(paste0("S", 1:6), "OTHER")
  X <- Matrix::Matrix(rbind(c(5, 0, 0, 0, 0, 0, 95), c(0, 5, 0, 0, 0, 0, 95), c(0, 0, 5, 0, 0, 0, 95), c(0, 0, 0, 5, 0, 0, 95), c(0, 0, 0, 0, 0, 5, 95), c(0, 0, 0, 0, 0, 0, 100)), sparse = TRUE)
  e <- data.frame(i = 1:5, j = 2:6, w = 1)                                                   # 사슬 1–2–3–4–5–6 (3–4 는 block 사이)
  M <- file.path(w, "M"); for (rel in c(NSO_GRAPH, NSO_H5AD)) dir.create(file.path(M, dirname(rel)), recursive = TRUE, showWarnings = FALSE)
  writeLines("graph-bytes", file.path(M, NSO_GRAPH)); writeLines(ids, file.path(M, NSO_GRAPH_IDS)); writeLines("h5-bytes", file.path(M, NSO_H5AD))
  O <- file.path(w, "O"); mk <- function(cl) data.frame(cell = 1:6, cell_id = ids, cluster = cl, block = donor, study = "s1", sorting = "Bulk", stringsAsFactors = FALSE)
  for (p in names(NSO_NRV)) { att <- file.path(O, NSO_NRV[[p]]); sd <- nrv_stage_dir(att, "structure"); dir.create(sd, recursive = TRUE)
    cl <- mk(if (p == "p92") c("A", "A", "B", "B", "C", "C") else rep("X", 6))
    saveRDS(list(cells = cl, bm = data.frame(block = c("d1", "d2"), study = "s1", sorting_class = "bulk", n_cells = c(3L, 3L)), folds = data.frame(), keys = data.frame()), file.path(sd, "cells_lossless.rds"))
    saveRDS(e, file.path(sd, "edges.rds")); saveRDS(list(full = nrv_pair_table(e, cl$cluster, cl$block, cl$study), folds = list()), file.path(sd, "pairs.rds"))
    mdu_write(data.frame(input = c("graph", "graph_ids", "h5ad"), path = file.path(M, c(NSO_GRAPH, NSO_GRAPH_IDS, NSO_H5AD)),
                         sha256 = c(sha(file.path(M, NSO_GRAPH)), sha(file.path(M, NSO_GRAPH_IDS)), sha(file.path(M, NSO_H5AD)))), file.path(sd, "input_sha256.tsv"))
    nrv_write_stage_manifest(att, "structure", files = file.path(sd, c("cells_lossless.rds", "edges.rds", "pairs.rds", "input_sha256.tsv")), upstream = "test") }
  # 기준선 (NSO attempt · connectivity): 기준 패널 크기 2, 기준 검출 = Myeloblast 세포 1 · 2
  bd <- file.path(w, "base"); dir.create(bd); nd <- matrix(0, 6, 6, dimnames = list(NULL, paste0("marker_", NSO_STAGES))); nd[1, 1] <- 2; nd[2, 1] <- 1
  saveRDS(list(metadata = data.frame(cell_id = ids, stringsAsFactors = FALSE), programs = list(cells = list(n_detected = nd)),
               reference = list(stage_markers = data.frame(stage_id = rep(NSO_STAGES, each = 2), gene = paste0("b", 1:12), stringsAsFactors = FALSE))), file.path(bd, "lossless.rds"))
  writeLines(jsonlite::toJSON(list(files = list(lossless.rds = sha(file.path(bd, "lossless.rds")))), auto_unbox = TRUE), file.path(bd, "manifest.json"))
  bc <- file.path(w, "basecon"); dir.create(bc)
  rt <- expand.grid(stage = NSO_STAGES, graph = c("all_edges", "same_block"), criterion = c("ge1", "ge2"), stringsAsFactors = FALSE)[, 3:1]
  rt$n_cells <- 7L; rt$n_reachable_cells <- 3L; rt$first_break <- ""; mdu_write(rt, file.path(bc, "ordered_reachability.tsv"))
  writeLines(jsonlite::toJSON(list(files = list(ordered_reachability.tsv = sha(file.path(bc, "ordered_reachability.tsv")))), auto_unbox = TRUE), file.path(bc, "manifest.json"))
  U <- cbind(c(0, 1, 2, 3, 4, 5), c(0, 0, 0, 0, 0, 0)); au <- c("GMP-Neut", "GMP-Neut", "Early GMP", NA, "CD14 Mono", "CD14 Mono")
  g_src <- function(what) if (!file.exists(file.path(out, "source_freeze.json"))) stop(what, " 를 source 고정 전에 읽었다")
  g_calc <- function(what) if (!file.exists(file.path(out, "calc_manifest.json"))) stop(what, " 를 계산 고정 전에 읽었다")
  reader <- list(ids = function() { g_src("graph ids"); ids }, obs = function(name) { g_src("obs"); switch(name, cell_id = ids, donor_id = donor, Study = rep("s1", 6), stop("obs: ", name)) },
                 symbols = function() { g_src("symbols"); tsym }, var_ids = function() { g_src("var ids"); tid }, raw_chunks = function(chunk) { g_src("raw"); nrv_chunks_matrix(X, chunk) },
                 umap = function() { g_calc("UMAP"); U }, author = function() { g_calc("저자 라벨"); au },
                 input_files = function() c(graph = file.path(M, NSO_GRAPH), graph_ids = file.path(M, NSO_GRAPH_IDS), h5ad = file.path(M, NSO_H5AD)))
  pins <- list(h5ad = sha(file.path(M, NSO_H5AD)), graph = sha(file.path(M, NSO_GRAPH)), baseline = sha(file.path(bd, "manifest.json")), baseline_connectivity = sha(file.path(bc, "manifest.json")))
  list(argv = c("--source-dir", src, "--source-manifest", smf, "--data-root", M, "--o-root", O, "--baseline-dir", bd, "--baseline-connectivity-dir", bc, "--out-dir", out, "--chunk", "4"),
       out = out, reader = reader, pins = pins, ids = ids, src = src)
}
atl_run <- function(w, pins = w$pins, argv = w$argv) { code <- NA_integer_; utils::capture.output(code <- atlas_main(argv, root = ROOT, reader = w$reader, pins = pins)); code }
atl_rd <- function(out, f) { p <- file.path(out, f); con <- if (grepl("\\.gz$", f)) gzfile(p) else p
  if (grepl("\\.csv(\\.gz)?$", f)) utils::read.csv(con, colClasses = "character", na.strings = character(0), check.names = FALSE)
  else utils::read.delim(con, colClasses = "character", quote = "", na.strings = character(0), check.names = FALSE) }

test_that("실행기 (소형 합성): source 고정 → target 적용 → 계산 고정 → 표시 · 사후 순서로 끝까지 돌고, 패널 · 가용성 · 연결 · 분모 · 비교가 손 계산과 같다", {
  w <- atl_world(); expect_equal(atl_run(w), 0L); out <- w$out
  src_files <- file.path("source", c("source_files.csv", "source_params.csv", "source_metadata.csv", "source_groups.csv", "source_contrasts_all.csv.gz", "source_transitions.csv",
                                     "source_panels.csv", "source_panel_status.csv", "source_exclusions.csv", "source_hgnc_map.csv.gz", "source_hgnc_summary.csv",
                                     "source_pb_pmn_groups.csv", "source_pb_pmn_panel_genes.csv", "source_reference.rds"))
  need <- c("source_freeze.json", src_files, "target_gene_availability.tsv", "cell_atlas_evidence.tsv.gz", "cluster_block_programs.tsv.gz", "cluster_block_genes.tsv.gz",
            "cluster_programs.tsv", "program_definitions.tsv", "detection_sets.tsv", "stage_components.tsv", "adjacent_pair_edges.tsv", "ordered_reachability.tsv",
            "baseline_comparison.tsv", "reach_comparison.tsv", "cluster_comparison.tsv", "calc_manifest.json", "figure_text.tsv", "posthoc_author_cluster_composition.tsv",
            "posthoc_author_program_summary.tsv", "checks.tsv", "input_sha256.tsv", "environment.tsv", "source_snapshot.tsv", "lossless.rds", "manifest.json", "README_ko.md",
            file.path("figures", c("fig01_atlas_panel_umap.png", "fig02_baseline_vs_atlas.png")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  pn <- atl_rd(out, "source/source_panels.csv"); ps <- atl_rd(out, "source/source_panel_status.csv")
  expect_equal(pn$symbol[pn$stage == "MY"], c("S3", "S8")); expect_equal(pn$symbol[pn$stage == "SN"], "S5"); expect_equal(ps$status[ps$stage == "BN"], "empty_no_candidate")
  expect_true(any(atl_rd(out, "source/source_exclusions.csv")$reason == "map_unmapped"))
  expect_equal(atl_rd(out, "source/source_pb_pmn_groups.csv")$code, "A"); expect_false(any(atl_rd(out, "source/source_groups.csv")$stage == "PB_PMN"))
  av <- atl_rd(out, "target_gene_availability.tsv")
  expect_equal(av$availability[av$symbol == "S8" & av$item_kind == "stage_panel"], "missing_gene"); expect_equal(av$availability[av$symbol == "S5" & av$item_kind == "stage_panel"], "observed_zero")
  expect_equal(av$availability[av$symbol == "S1" & av$item_kind == "stage_panel"], "observed"); expect_equal(av$n_cells_detected[av$symbol == "S8" & av$item_kind == "stage_panel"], "NA")
  ds <- atl_rd(out, "detection_sets.tsv"); my <- ds[ds$stage == "MY", ]
  expect_equal(c(my$n_panel_genes, my$n_observed_genes, my$missing_genes, my$status_ge1), c("2", "1", "S8", "available"))
  expect_equal(ds$status_ge1[ds$stage == "BN"], "unavailable_empty_panel"); expect_equal(ds$n_cells_ge1[ds$stage == "SN"], "0")
  rt <- atl_rd(out, "ordered_reachability.tsv"); a <- rt[rt$criterion == "ge1" & rt$graph == "all_edges", ]; s <- rt[rt$criterion == "ge1" & rt$graph == "same_block", ]
  expect_equal(a$n_reachable_cells, c("1", "1", "1", "1", "NA", "NA")); expect_equal(a$stage_status, c(rep("available", 4), "unavailable_empty_panel", "not_evaluable_after_unavailable"))
  expect_equal(a$n_cells[a$stage == "SN"], "0"); expect_equal(unique(a$first_break), "MM>BN"); expect_equal(unique(a$break_kind), "unavailable_stage")
  expect_equal(unique(s$first_break), "MY>MM"); expect_equal(unique(s$break_kind), "no_real_edge")       # 3–4 는 block 사이 edge 뿐
  g2 <- rt[rt$criterion == "ge2" & rt$graph == "all_edges", ]                                     # 관측 패널 유전자가 1 개뿐인 단계는 ≥2 기준에서 판정불가
  expect_equal(g2$stage_status, c(rep("unavailable_too_few_genes", 4), "unavailable_empty_panel", "unavailable_too_few_genes")); expect_equal(unique(g2$first_break), "start>Myeloblast")
  ce <- atl_rd(out, "cell_atlas_evidence.tsv.gz"); expect_equal(ce$cell_id, w$ids); expect_false(any(grepl("stage_label|predicted|author", names(ce))))
  expect_equal(ce$marker_BN_n_detected, rep("NA", 6)); expect_equal(ce$marker_SN_n_detected, rep("0", 6))
  cb <- atl_rd(out, "cluster_block_programs.tsv.gz"); expect_true(all(tapply(as.numeric(cb$n_cells), paste(cb$partition, cb$program_id), sum) == 6))
  cm <- atl_rd(out, "baseline_comparison.tsv"); r <- cm[cm$stage == "Myeloblast", ]
  expect_equal(as.numeric(c(r$baseline_panel_size, r$atlas_panel_size, r$baseline_n_ge1, r$atlas_n_ge1, r$overlap_ge1, r$baseline_reach_all_edges)), c(2, 1, 2, 1, 1, 3))
  expect_equal(cm$atlas_n_ge1[cm$stage == "BN"], "NA"); expect_equal(cm$atlas_panel_status[cm$stage == "BN"], "empty_no_candidate")
  expect_equal(cm$atlas_panel_in_target[cm$stage == "MY"], "1")
  cc <- atl_rd(out, "cluster_comparison.tsv"); z <- cc[cc$partition == "p92" & cc$cluster == "A" & cc$stage == "Myeloblast", ]
  expect_equal(as.numeric(c(z$n_cells, z$baseline_frac_ge1, z$atlas_frac_ge1)), c(2, 1, 0.5))
  ck <- atl_rd(out, "checks.tsv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true(all(c("source_frozen_before_target", "denominators_preserved", "baseline_cells_identical_order", "author_and_umap_read_after_fix") %in% ck$check))
  cmf <- names(jsonlite::fromJSON(file.path(out, "calc_manifest.json"))$files); expect_false(any(grepl("posthoc|figure|^source/", cmf)))
  sfz <- names(jsonlite::fromJSON(file.path(out, "source_freeze.json"))$files); expect_true(all(startsWith(sfz, "source/"))); expect_setequal(sfz, src_files)
  man <- jsonlite::fromJSON(file.path(out, "manifest.json"))$files
  expect_true(all(vapply(names(man), function(f) identical(digest::digest(file.path(out, f), algo = "sha256", file = TRUE), man[[f]]), TRUE)))
  expect_true(all(src_files %in% names(man)))
  ft <- atl_rd(out, "figure_text.tsv"); expect_false(any(ft$text %in% c("NA", "NaN", "Inf")))
})

test_that("실행기: 덮어쓰기 거절, 고정 hash · source 파일 불일치 중지, 인자 누락", {
  w <- atl_world(); expect_equal(atl_run(w), 0L); expect_equal(atl_run(w), 6L)
  w2 <- atl_world(); bad <- w2$pins; bad$baseline <- strrep("0", 64); expect_equal(atl_run(w2, pins = bad), 4L); expect_true(any(grepl("hash", readLines(file.path(w2$out, "run.log")))))
  expect_false(file.exists(file.path(w2$out, "calc_manifest.json")))
  w3 <- atl_world(); cat("x\n", file = file.path(w3$src, "GSE294330_full.soft"), append = TRUE); expect_equal(atl_run(w3), 4L)
  expect_true(any(grepl("source", readLines(file.path(w3$out, "run.log"))))); expect_false(file.exists(file.path(w3$out, "source_freeze.json")))
  w4 <- atl_world(); expect_equal(atl_run(w4, argv = w4$argv[1:4]), 2L)
  w5 <- atl_world(); bad <- w5$pins; bad$h5ad <- strrep("1", 64); expect_equal(atl_run(w5, pins = bad), 4L)
})
