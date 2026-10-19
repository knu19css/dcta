# =============================================================================
# test_atlas_rna_reference.R — atlas RNA source 전용 대비 · 패널 (R63) 테스트
#   계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md, design.md, plan Task 1
#   손 계산 예제: 같은 code x 단계 중복 열(중앙값), code 의 단계 결측(쌍 분모), TPM 0, 부호 섞임(일치율), 높은 쪽 TPM 부족,
#   HGNC 다중 대응 · 비코딩 제외(조용히 지우지 않고 제외 표), 빈 중간 단계 패널, 동률 ENSG 정렬, SOFT 열 대응 · 단계 · code · 병합 표기
# =============================================================================
library(testthat)
atl_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- atl_test_root()
for (f in c("58_nso_inputs.R", "63_atlas_rna_reference.R")) source(file.path(ROOT, "R", "branch_biology", f))

E <- sprintf("ENSG%011d", 1:9)
# 3 단계(Myeloblast · PM · MY), code A–D. A 의 Myeloblast 는 열 2 개(중복), D 는 PM 이 없다. log2(TPM + 1): 0 · 1 · 3 · 7 · 15 · 31 → 0 · 1 · 2 · 3 · 4 · 5
atl_toy <- function() {
  cols <- c("A_MB_1", "A_MB_2", "B_MB", "C_MB", "D_MB", "A_PM", "B_PM", "C_PM", "A_MY", "B_MY", "C_MY", "D_MY")
  v <- rbind(c(7, 23, 15, 15, 15, 1, 1, 1, 1, 1, 1, 1),          # E1 Myeloblast 높음 (A 는 중앙값 15)
             c(0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7),               # E2 MY 높음 (TPM 0 포함)
             c(3, 3, 3, 3, 3, 15, 15, 0, 3, 3, 3, 3),             # E3 부호 섞임 (일치율 2/3)
             c(1, 1, 1, 1, 1, 15, 15, 15, 12, 12, 12, 12),        # E4 PM 이 Myeloblast 보다만 높음 (PM→MY 차이 -0.30)
             c(0, 0, 0, 0, 0, 31, 31, 31, 31, 31, 31, 31),        # E5 HGNC 다중 대응
             c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),               # E6 HGNC 미대응 · 모두 0
             c(15, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0),          # E7 비코딩 (Myeloblast 높음)
             c(0, 0, 0, 0, 0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),  # E8 높은 쪽 TPM < 1
             c(0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7))               # E9 E2 와 같은 값 (동률 → ENSG 순서)
  tpm <- matrix(v, nrow = 9, dimnames = list(E, cols))
  md <- data.frame(column = cols, gsm = paste0("GSM", seq_along(cols)), title = cols, cell_type = "", stage = rep(c("Myeloblast", "PM", "MY"), c(5, 3, 4)),
                   code = c("A", "A", "B", "C", "D", "A", "B", "C", "A", "B", "C", "D"), subfraction = "", merge_note = "", in_bm = TRUE, stringsAsFactors = FALSE)
  map <- data.frame(ensg = E, symbol = c("G1", "G2", "G3", "G4", "", "", "G7", "G8", "G9"), hgnc_id = "", n_hgnc_rows = c(1L, 1L, 1L, 1L, 2L, 0L, 1L, 1L, 1L),
                    locus_group = c(rep("protein-coding gene", 4), "", "", "non-coding RNA", "protein-coding gene", "protein-coding gene"),
                    map_status = c("one_to_one", "one_to_one", "one_to_one", "one_to_one", "multi_hgnc", "unmapped", "one_to_one", "one_to_one", "one_to_one"), stringsAsFactors = FALSE)
  map$use_for_panels <- map$map_status == "one_to_one" & map$locus_group == "protein-coding gene"
  list(tpm = tpm, metadata = md, map = map)
}
ST3 <- c("Myeloblast", "PM", "MY")

test_that("code x 단계 요약: 같은 code x 단계의 여러 열은 TPM 중앙값으로 묶고, 없는 단계는 만들지 않는다", {
  g <- atlas_code_stage_tpm(atl_toy()$tpm, atl_toy()$metadata)
  expect_equal(unname(g$tpm_med[E[1], "A|Myeloblast"]), 15); expect_equal(g$groups$n_columns[g$groups$group == "A|Myeloblast"], 2L)
  expect_false("D|PM" %in% g$groups$group); expect_equal(nrow(g$groups), 11L)
})

test_that("인접 대비: 쌍이 있는 code 만 쓰고(분모), 중앙값 · 부호 일치율 · 높은 쪽 TPM 코드 수 · 후보 판정이 손 계산과 같다", {
  ref <- atlas_build_reference(atl_toy(), stages = ST3); a <- ref$all_contrasts; r <- function(p, e) a[a$pair == p & a$ensg == e, ]
  expect_equal(unique(a$n_pairs[a$pair == "Myeloblast>PM"]), 3L); expect_equal(r("Myeloblast>PM", E[1])$codes, "A|B|C")   # D 는 PM 이 없다
  expect_equal(r("Myeloblast>PM", E[1])$median_diff, -3); expect_equal(r("Myeloblast>PM", E[1])$direction, "down"); expect_true(r("Myeloblast>PM", E[1])$candidate)
  expect_equal(r("Myeloblast>PM", E[3])$median_diff, 2); expect_equal(r("Myeloblast>PM", E[3])$sign_agree, 2 / 3); expect_false(r("Myeloblast>PM", E[3])$candidate)
  expect_equal(r("Myeloblast>PM", E[8])$n_codes_high_tpm, 0L); expect_false(r("Myeloblast>PM", E[8])$candidate)              # 차이 0.85 이지만 PM TPM 0.8 < 1
  expect_equal(r("PM>MY", E[4])$median_diff, log2(13) - 4); expect_false(r("PM>MY", E[4])$candidate)                        # |차이| < 0.5
  expect_false(r("Myeloblast>PM", E[6])$candidate); expect_equal(r("Myeloblast>PM", E[6])$median_diff, 0)                    # TPM 0 끼리
  expect_equal(r("Myeloblast>PM", E[1])$diff_A, -3); expect_true(is.na(r("Myeloblast>PM", E[1])$diff_D)); expect_equal(r("Myeloblast>PM", E[1])$tpm_from_A, 15)
  expect_equal(nrow(a), 9L * 2L)                                                              # 모든 유전자 x 모든 인접 쌍을 보존
})

test_that("선택 변화와 단계 패널: 일대일 protein-coding 만 쓰고 제외는 이유와 함께 남기며, 빈 중간 단계 패널은 억지로 채우지 않는다", {
  ref <- atlas_build_reference(atl_toy(), stages = ST3); tr <- ref$transitions
  expect_equal(tr$ensg[tr$pair == "Myeloblast>PM" & tr$direction == "up"], E[4]); expect_equal(tr$ensg[tr$pair == "Myeloblast>PM" & tr$direction == "down"], E[1])
  expect_equal(tr$ensg[tr$pair == "PM>MY" & tr$direction == "up"], c(E[2], E[9])); expect_equal(sum(tr$pair == "PM>MY" & tr$direction == "down"), 0L)
  ex <- ref$exclusions
  expect_equal(ex$reason[ex$ensg == E[5] & ex$kind == "transition"], "map_multi_hgnc"); expect_true(all(ex$reason[ex$ensg == E[7]] == "not_protein_coding"))
  expect_true(any(ex$kind == "panel" & ex$ensg == E[7] & ex$stage_or_pair == "Myeloblast"))
  p <- ref$panels; ps <- ref$panel_status
  expect_equal(p$ensg[p$stage == "Myeloblast"], E[1]); expect_equal(p$min_adjacent_effect[p$stage == "Myeloblast"], 3)
  expect_equal(p$ensg[p$stage == "MY"], c(E[2], E[9]))                                       # 동률 효과 → ENSG 순서
  expect_equal(sum(p$stage == "PM"), 0L); expect_equal(ps$status[ps$stage == "PM"], "empty_no_candidate")   # E4 는 Myeloblast 쪽만 만족
  expect_equal(ps$neighbors[ps$stage == "PM"], "Myeloblast|MY"); expect_equal(ps$neighbors[ps$stage == "Myeloblast"], "PM")
  r1 <- atlas_build_reference(atl_toy(), params = utils::modifyList(ATLAS_PARAMS, list(top_panel = 1L)), stages = ST3)
  expect_equal(r1$panels$ensg[r1$panels$stage == "MY"], E[2])                                 # 상한은 순위대로 자른다
})

test_that("쌍의 code 가 3 개보다 적으면 그 대비에서는 후보가 없다 (분모 기준)", {
  s <- atl_toy(); keep <- !(s$metadata$stage == "MY" & s$metadata$code %in% c("C", "D")); s$tpm <- s$tpm[, keep]; s$metadata <- s$metadata[keep, ]
  ref <- atlas_build_reference(s, stages = ST3)
  expect_equal(unique(ref$all_contrasts$n_pairs[ref$all_contrasts$pair == "PM>MY"]), 2L); expect_false(any(ref$all_contrasts$candidate[ref$all_contrasts$pair == "PM>MY"]))
  expect_equal(ref$panel_status$status[ref$panel_status$stage == "MY"], "empty_no_candidate")
})

test_that("HGNC 대응: 정확한 ENSG 일대일만 one_to_one, 여러 행이면 multi_hgnc, 없으면 unmapped", {
  hg <- data.frame(symbol = c("G1", "G5a", "G5b", "G7"), ensembl_gene_id = c(E[1], E[5], E[5], E[7]), locus_group = c("protein-coding gene", "protein-coding gene", "protein-coding gene", "non-coding RNA"),
                   hgnc_id = paste0("HGNC:", 1:4), stringsAsFactors = FALSE)
  m <- atlas_map_hgnc(E[c(1, 5, 6, 7)], hg)
  expect_equal(m$map_status, c("one_to_one", "multi_hgnc", "unmapped", "one_to_one")); expect_equal(m$symbol, c("G1", "", "", "G7"))
  expect_equal(m$use_for_panels, c(TRUE, FALSE, FALSE, FALSE)); expect_equal(m$n_hgnc_rows, c(1L, 2L, 0L, 1L))
})

atl_write_soft <- function(path, samples) {
  lines <- c("^SERIES = GSE0", unlist(lapply(samples, function(s) c(paste0("^SAMPLE = ", s$gsm), paste0("!Sample_title = ", s$title), paste0("!Sample_characteristics_ch1 = cell type: ", s$ct),
    paste0("!Sample_description = Library name: ", s$col), if (!is.null(s$col)) paste0("!Sample_description = Column name in raw_counts_mrn.csv: ", s$col), "!Sample_library_strategy = RNA-Seq"))))
  writeLines(lines, path)
}
test_that("source 읽기: SOFT 의 열 이름으로 행렬 열을 대응하고 단계 · code · 세부 분획 · 병합 표기를 읽으며, 대응이 어긋나면 멈춘다", {
  d <- tempfile("atl_src"); dir.create(d)
  tp <- data.frame(MC_a_2 = c(1, 2), R1_eB_C = c(3, 4), PMN_a = c(5, 6), ENSG = E[1:2], check.names = FALSE)
  con <- gzfile(file.path(d, "GSE294330_exp_mrn_tpm.csv.gz"), "w"); utils::write.csv(tp, con, row.names = FALSE); close(con)
  con <- gzfile(file.path(d, "GSE294330_exp_mrn_counts.csv.gz"), "w"); utils::write.csv(tp, con, row.names = FALSE); close(con)
  atl_write_soft(file.path(d, "GSE294330_full.soft"), list(list(gsm = "GSM1", title = "mRNA.cellstage_Myelocyte.early_donor_A", ct = "Myelocyte", col = "MC_a_2"),
    list(gsm = "GSM2", title = "mRNA.cellstage_Band_donor_C_mergedInFinalWith_eB_c", ct = "Band", col = "R1_eB_C"),
    list(gsm = "GSM3", title = "mRNA.cellstage_Neutrophil_donor_A", ct = "Neutrophil", col = "PMN_a"),
    list(gsm = "GSM4", title = "miRNA.cellstage.Band.unspecified_donor.A", ct = "Band", col = NULL)))
  utils::write.table(data.frame(symbol = c("G1", "G2"), ensembl_gene_id = E[1:2], locus_group = "protein-coding gene", hgnc_id = c("HGNC:1", "HGNC:2"), status = "Approved"),
                     file.path(d, "hgnc_complete_set.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  s <- atlas_read_source(d); m <- s$metadata
  expect_equal(m$column, c("MC_a_2", "R1_eB_C", "PMN_a")); expect_equal(m$stage, c("MY", "BN", "PB_PMN")); expect_equal(m$code, c("A", "C", "A"))
  expect_equal(m$subfraction, c("early", "", "")); expect_equal(m$merge_note, c("", "eB_c", "")); expect_equal(m$in_bm, c(TRUE, TRUE, FALSE)); expect_equal(m$gsm, c("GSM1", "GSM2", "GSM3"))
  expect_equal(dim(s$tpm), c(2L, 3L)); expect_equal(rownames(s$tpm), E[1:2]); expect_true(all(s$map$map_status == "one_to_one"))
  atl_write_soft(file.path(d, "GSE294330_full.soft"), list(list(gsm = "GSM1", title = "mRNA.cellstage_Myelocyte_donor_A", ct = "Myelocyte", col = "MC_a_2")))
  expect_error(atlas_read_source(d), "R1_eB_C")                                               # 행렬 열이 SOFT 에 없으면 멈춘다
})

test_that("실제 GSE294330 (있으면): mRNA 40 열, 단계별 열 수 · code, PMN 분리, 인접 대비 code 수 4 · 4 · 4 · 4 · 3", {
  d <- file.path(ROOT, "data", "external", "granulopoiesis_atlas_2026-10-19"); skip_if_not(file.exists(file.path(d, "GSE294330_full.soft")))
  s <- atlas_read_source(d); m <- s$metadata
  expect_equal(nrow(m), 40L); expect_equal(as.integer(table(m$stage)[c("Myeloblast", "PM", "MY", "MM", "BN", "SN", "PB_PMN")]), c(4L, 4L, 6L, 8L, 9L, 4L, 5L))
  cs <- function(st) paste(sort(unique(m$code[m$stage == st])), collapse = "")
  expect_equal(vapply(c("Myeloblast", "PM", "MY", "MM", "BN", "SN"), cs, ""), c(Myeloblast = "ABCE", PM = "ABCE", MY = "ABCDE", MM = "ABDE", BN = "ABCDE", SN = "ACE"))
  expect_equal(sum(m$in_bm), 35L); expect_equal(nrow(s$tpm), 61552L)
  g <- atlas_code_stage_tpm(s$tpm[1:5, m$in_bm], m[m$in_bm, ]); n <- function(a, b) length(intersect(g$groups$code[g$groups$stage == a], g$groups$code[g$groups$stage == b]))
  expect_equal(c(n("Myeloblast", "PM"), n("PM", "MY"), n("MY", "MM"), n("MM", "BN"), n("BN", "SN")), c(4L, 4L, 4L, 4L, 3L))
})
