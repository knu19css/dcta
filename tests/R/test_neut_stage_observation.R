# =============================================================================
# test_neut_stage_observation.R — 호중구 성숙 단계 관측 (R58–61) 테스트
#   계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md
#   Task 1: reference 6 단계 · 순서 · 공유 marker · PB 분리 · 0 행 오류, 입력 ID 검사, 계산용 메타데이터 정보 분리
#   Task 2: 정규화 · 평균 logRNA · 검출률 · 경계 평균 손 계산 대조, missing_gene / 관측 0, 공유 marker, 단일 유전자 프로그램,
#           빈 경계, 한 공여자, 관측 쌍 없음, 방향 = 부호 반전, B 내부 단절 (R24 nrv_internal_coordinate 와 같은 상태), spot-check 원 count
#   Task 3: 소형 합성 자료로 실행기 전체 — 결과 폴더 중복 중지, 전체 cell ID · 분모 보존, 원 UMAP 순서 · 세포 수 대조,
#           미검출 프로그램 · 빈 edge 그림 경로 (글자 'NA' 없음), 결과 고정 전 저자 라벨 · UMAP 입력 금지, 고정 hash 불일치 중지
# =============================================================================
library(testthat)
nso_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- nso_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
            "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
            "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

# ── 고정 예제: 실제 02_nodes.csv / 03_edges.csv 와 같은 열 (값은 손으로 줄였다) ──
nso_fixture_nodes <- function() {
  ids <- c("HSC", "MPP", "CMP", "GMP", "Myeloblast", "PM", "MY", "MM", "BN", "SN", "PB_PMN")
  mk <- c("HLF;HOPX", "CDK6;MYC", "FLT3;CD164", "SPI1;CEBPA", "MPO;ELANE;PRTN3", "DEFA4;BPI", "LTF;CAMP;LCN2", "LTF;LCN2;MMP9", "MMP9;FCGR3B",
          "FCGR3B;FPR1;CXCR2", "FCGR3B;FPR1;CXCR2")
  data.frame(node_id = ids, stage_order = 1:11, display_name = ids, full_name = ids, description_ko = "설명",
             compartment = c(rep("bone_marrow", 10), "peripheral_blood"), name = ids, marker_genes = mk, program_ko = "프로그램", state_note_ko = "메모",
             reference_doi = "10.0/x", adult_evidence_scope = c(rep("adult_human_bone_marrow_confirmed", 4), rep("human_bone_marrow_or_blood; donor_age_not_confirmed", 7)),
             x = 0, y = 0, fill_color = "#FFFFFF", node_label = ids, stringsAsFactors = FALSE, check.names = FALSE)
}
nso_fixture_edges <- function() {
  fr <- c("HSC", "MPP", "CMP", "GMP", "Myeloblast", "PM", "MY", "MM", "BN", "SN"); to <- c("MPP", "CMP", "GMP", "Myeloblast", "PM", "MY", "MM", "BN", "SN", "PB_PMN")
  data.frame(edge_id = sprintf("E%02d", 1:10), from = fr, to = to,
             edge_type = c(rep("classical_model_projection", 3), "classical_to_morphology_projection", rep("maturation", 5), "marrow_to_blood"), process_ko = "과정",
             up_genes = c("", "", "", "", "DEFA4;BPI", "LTF;CAMP;LCN2", "MMP9", "MMP9", "", ""),
             down_genes = c("", "", "", "", "", "MPO;ELANE;PRTN3", "MKI67;TOP2A", "LTF;CAMP;LCN2", "CEBPE;MMP9", "CEBPA;CEBPE;LTF;MMP9"),
             conditional_up_genes = c("CDK6", "FLT3", "SPI1", "MPO;ELANE", "", "", "", "FCGR3B;FPR1", "FPR1;CXCR2", ""),
             conditional_down_genes = c("HLF", "HLF", "GATA2", "", "MPO;ELANE", "", "", "", "", ""), note_ko = "메모",
             direction_basis = c(rep("adult_human_RNA_programs_projected_to_classical_edge", 3), "candidate_mapping_between_progenitor_and_morphology",
                                 "granule_wave_with_conditional_decline", "granule_wave_literature", "maturation_program_literature",
                                 "granule_wave_and_candidate_receptor_trend", "reported_decline_and_candidate_receptor_trend", "marrow_blood_transcriptome_comparison"),
             adult_evidence_scope = "범위", reference_doi = "10.0/x", source = fr, target = to, interaction = c(rep("differentiates_to", 9), "marrow_to_blood"),
             name = paste(fr, to), edge_label = "", line_type = "SOLID", label_position = "C", stringsAsFactors = FALSE, check.names = FALSE)
}
nso_write_fixture <- function(nodes = nso_fixture_nodes(), edges = nso_fixture_edges()) {
  d <- tempfile("nso_ref"); dir.create(d)
  utils::write.csv(nodes, file.path(d, "02_nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(edges, file.path(d, "03_edges.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  c(nodes = file.path(d, "02_nodes.csv"), edges = file.path(d, "03_edges.csv"))
}
nso_fixture_ref <- function(...) { p <- nso_write_fixture(...); nso_reference(p[["nodes"]], p[["edges"]]) }


# ── Task 1: reference ────────────────────────────────────────────────────────
test_that("reference: 형태 6 단계를 Myeloblast→PM→MY→MM→BN→SN 순서로 고정하고, 단계 사이 변화는 E05–E09 만, PB 이동(E10)과 PB_PMN 은 따로 둔다", {
  ref <- nso_fixture_ref()
  expect_equal(ref$stages$stage_id, c("Myeloblast", "PM", "MY", "MM", "BN", "SN")); expect_equal(ref$stages$stage_rank, 1:6)
  expect_equal(ref$stages$source_stage_order, 5:10)
  expect_equal(sort(unique(ref$changes$edge_id)), c("E05", "E06", "E07", "E08", "E09"))
  expect_false("PB_PMN" %in% c(ref$stages$stage_id, ref$stage_markers$stage_id, ref$changes$to))
  expect_equal(ref$pb$node$node_id, "PB_PMN"); expect_equal(unique(ref$pb$changes$edge_id), "E10")
  expect_equal(sort(ref$pb$changes$gene), c("CEBPA", "CEBPE", "LTF", "MMP9")); expect_true(all(ref$pb$changes$direction == "down"))
  e06 <- ref$changes[ref$changes$edge_id == "E06", ]
  expect_equal(e06$gene[e06$direction == "up"], c("LTF", "CAMP", "LCN2")); expect_equal(e06$gene[e06$direction == "down"], c("MPO", "ELANE", "PRTN3"))
  expect_true(all(e06$evidence == "core"))
  e05 <- ref$changes[ref$changes$edge_id == "E05", ]
  expect_equal(e05$gene[e05$evidence == "conditional"], c("MPO", "ELANE")); expect_true(all(e05$direction[e05$evidence == "conditional"] == "down"))
  tr <- ref$transitions
  expect_equal(paste(tr$from, tr$to), c("Myeloblast PM", "PM MY", "MY MM", "MM BN", "BN SN")); expect_equal(tr$edge_id, sprintf("E%02d", 5:9))
  expect_equal(tr$n_core_up, c(2L, 3L, 1L, 1L, 0L)); expect_equal(tr$n_conditional_up, c(0L, 0L, 0L, 2L, 2L))
  e <- nso_fixture_edges(); tr2 <- nso_fixture_ref(edges = e[e$edge_id != "E07", ])$transitions          # 빠진 단계 사이 연결은 지우지 않고 표시
  expect_equal(tr2$status[tr2$from == "MY"], "no_reference_edge"); expect_equal(tr2$edge_id[tr2$from == "MY"], "")
  sc <- ref$scope
  expect_equal(nrow(sc), 21L)                                                 # 11 node + 10 edge 모두 행으로 남는다
  expect_equal(sc$role[sc$id == "E04"], "entry_projection"); expect_equal(sc$role[sc$id == "E10"], "marrow_to_blood")
  expect_equal(sc$role[sc$id == "E02"], "upstream_projection"); expect_equal(sc$role[sc$id == "E07"], "stage_transition")
  expect_equal(sc$role[sc$id == "GMP"], "upstream_context"); expect_equal(sc$role[sc$id == "MM"], "morphology_stage")
  expect_equal(sc$role[sc$id == "PB_PMN"], "marrow_to_blood_target")
})

test_that("reference: 공유 marker 를 지우지 않는다 (LTF·LCN2 = MY·MM, MMP9 = MM·BN, FCGR3B = BN·SN), PB 와의 공유는 따로 표시한다", {
  ref <- nso_fixture_ref(); sm <- ref$shared_markers; mk <- ref$stage_markers
  expect_equal(sm$stages[sm$gene == "LTF"], "MY|MM"); expect_equal(sm$stages[sm$gene == "LCN2"], "MY|MM")
  expect_equal(sm$stages[sm$gene == "MMP9"], "MM|BN"); expect_equal(sm$stages[sm$gene == "FCGR3B"], "BN|SN")
  expect_false(any(c("MPO", "CAMP", "DEFA4") %in% sm$gene))                  # 한 단계에만 있고 PB 와도 공유하지 않는다
  expect_equal(sm$n_stages[sm$gene == "FPR1"], 1L); expect_true(sm$pb_shared[sm$gene == "FPR1"]); expect_false(sm$pb_shared[sm$gene == "LTF"])
  expect_equal(sum(mk$gene == "LTF"), 2L)                                     # 두 단계 행이 모두 남는다
  expect_equal(mk$shared_with[mk$stage_id == "MY" & mk$gene == "LTF"], "MM"); expect_equal(mk$shared_with[mk$stage_id == "MY" & mk$gene == "CAMP"], "")
  expect_equal(ref$stages$n_markers, c(3L, 2L, 3L, 3L, 2L, 3L))
})

test_that("reference: 단계 순서가 틀리거나 단계가 빠지거나, 단계 사이 edge 가 다음 단계로 가지 않으면 멈춘다", {
  n <- nso_fixture_nodes(); n$stage_order[n$node_id == "MY"] <- 9L; n$stage_order[n$node_id == "BN"] <- 7L
  expect_error(nso_fixture_ref(nodes = n), "순서")
  n <- nso_fixture_nodes(); expect_error(nso_fixture_ref(nodes = n[n$node_id != "MM", ]), "MM")
  e <- nso_fixture_edges(); e$to[e$edge_id == "E06"] <- "MM"; expect_error(nso_fixture_ref(edges = e), "E06")      # PM → MM: 한 단계를 건너뜀
  e <- nso_fixture_edges(); e$from[e$edge_id == "E07"] <- "MM"; e$to[e$edge_id == "E07"] <- "MY"; expect_error(nso_fixture_ref(edges = e), "E07")   # 역방향
  e <- nso_fixture_edges(); e$edge_id[e$edge_id == "E08"] <- "E07"; expect_error(nso_fixture_ref(edges = e), "중복")
})

test_that("reference: 0 행 표와 필수 열 누락은 이유와 함께 멈추고, marker 가 빈 단계는 지우지 않고 표시한다", {
  expect_error(nso_fixture_ref(nodes = nso_fixture_nodes()[0, ]), "0 행")
  expect_error(nso_fixture_ref(edges = nso_fixture_edges()[0, ]), "0 행")
  n <- nso_fixture_nodes(); n$marker_genes <- NULL; expect_error(nso_fixture_ref(nodes = n), "marker_genes")
  e <- nso_fixture_edges(); e$conditional_up_genes <- NULL; expect_error(nso_fixture_ref(edges = e), "conditional_up_genes")
  n <- nso_fixture_nodes(); n$marker_genes[n$node_id == "BN"] <- ""; ref <- nso_fixture_ref(nodes = n)
  expect_equal(nrow(ref$stages), 6L); expect_equal(ref$stages$n_markers[ref$stages$stage_id == "BN"], 0L)
  expect_equal(ref$stages$marker_status[ref$stages$stage_id == "BN"], "no_reference_markers"); expect_equal(ref$stages$marker_status[ref$stages$stage_id == "MM"], "markers_listed")
  expect_false("BN" %in% ref$stage_markers$stage_id); expect_false("MMP9" %in% ref$shared_markers$gene)   # BN 이 비면 MMP9 는 MM 에만 남는다
})

test_that("reference (실제 참고표가 있으면): 6 단계, PM 유전자 순서, E05–E09 26 행, 공유 4 유전자, PB 분리, 기준 hash", {
  d <- NRV_REF_DIR; skip_if_not(file.exists(file.path(d, "02_nodes.csv")))
  ref <- nso_reference(file.path(d, "02_nodes.csv"), file.path(d, "03_edges.csv"))
  expect_equal(ref$stages$stage_id, c("Myeloblast", "PM", "MY", "MM", "BN", "SN")); expect_equal(ref$stages$source_stage_order, 5:10)
  expect_equal(ref$stage_markers$gene[ref$stage_markers$stage_id == "PM"], c("DEFA4", "BPI", "DEFA1", "DEFA3"))
  expect_equal(nrow(ref$changes), 26L); expect_equal(as.integer(table(ref$changes$edge_id)[c("E05", "E06", "E07", "E08", "E09")]), c(6L, 6L, 4L, 6L, 4L))
  expect_equal(sort(ref$shared_markers$gene[ref$shared_markers$n_stages >= 2]), c("FCGR3B", "LCN2", "LTF", "MMP9"))
  expect_equal(sort(ref$shared_markers$gene[ref$shared_markers$pb_shared]), c("CXCR2", "FCGR3B", "FPR1"))
  expect_equal(ref$files$sha256[ref$files$file == "02_nodes.csv"], "8059caa90f7bf67561f2fd89603b493032435a295a34c2e327209adec83368f8")
  expect_equal(ref$files$sha256[ref$files$file == "03_edges.csv"], "73f9f61fe35cdac0c9162bca930c6294737712c5e2973685f76c8b44eeadb58b")
})


# ── Task 1: 입력 검사 · 계산용 메타데이터 ─────────────────────────────────────
nso_mk_cells <- function(ids, cl, block = c("b1", "b1", "b2", "b2")) data.frame(cell = seq_along(ids), cell_id = ids, cluster = cl, block = block[seq_along(ids)],
                                                                             study = "s1", sorting = "Bulk", stringsAsFactors = FALSE)

test_that("입력 검사: 네 ID 목록(p92, p51, raw, graph)이 중복 없이 같은 순서일 때만 통과하고, 재정렬하지 않고 멈춘다", {
  ref <- nso_fixture_ref(); ids <- paste0("c", 1:4)
  c92 <- nso_mk_cells(ids, c("0", "0", "1", "2")); c51 <- nso_mk_cells(ids, c("0", "0", "1", "1"))
  ck <- nso_check_inputs(c92, c51, raw_ids = ids, graph_ids = ids, reference = ref)
  expect_equal(ck$cell_id, ids); expect_true(all(ck$checks$pass)); expect_true(nrow(ck$checks) >= 5L)
  expect_equal(ck$cells$cluster92, c("0", "0", "1", "2")); expect_equal(ck$cells$cluster51, c("0", "0", "1", "1")); expect_equal(ck$cells$block, c("b1", "b1", "b2", "b2"))
  expect_error(nso_check_inputs(c92, c51, raw_ids = rev(ids), graph_ids = ids, reference = ref), "순서")
  expect_error(nso_check_inputs(c92, c51, raw_ids = ids, graph_ids = ids[1:3], reference = ref), "순서")
  d <- c92; d$cell_id[2] <- "c1"; expect_error(nso_check_inputs(d, c51, ids, ids, ref), "중복")
  d <- c51; d$block[1] <- "b9"; expect_error(nso_check_inputs(c92, d, ids, ids, ref), "block")
  d <- c92; d$cluster[3] <- NA; expect_error(nso_check_inputs(d, c51, ids, ids, ref), "군집")
  d <- c92; d$cell <- c(1L, 3L, 2L, 4L); expect_error(nso_check_inputs(d, c51, ids, ids, ref), "cell")
  d <- c92; d$author_label <- "x"; expect_error(nso_check_inputs(d, c51, ids, ids, ref), "author_label", fixed = TRUE)
  r2 <- ref; r2$stages <- r2$stages[c(2, 1, 3:6), ]; expect_error(nso_check_inputs(c92, c51, ids, ids, r2), "순서")
})

test_that("계산용 메타데이터: cell_id · cluster92 · cluster51 · donor · study · block 만 돌려주고, 저자 라벨 · UMAP 열이 들어오면 멈춘다", {
  ids <- paste0("c", 1:4)
  cells <- data.frame(cell_id = ids, cluster92 = c("0", "0", "1", "2"), cluster51 = c("0", "0", "1", "1"), donor = c("d1", "d1", "d2", "d2"), study = "s1",
                      block = c("d1", "d1", "d2", "d2"), sorting = "Bulk", cell = 1:4, stringsAsFactors = FALSE)
  md <- nso_calculation_metadata(cells)
  expect_equal(names(md), c("cell_id", "cluster92", "cluster51", "donor", "study", "block")); expect_equal(md$cell_id, ids); expect_equal(md$cluster92, cells$cluster92)
  for (bad in c("author_label", "AuthorCellType", "umap_1", "X_umap", "cell_type")) { d <- cells; d[[bad]] <- "x"; expect_error(nso_calculation_metadata(d), bad, fixed = TRUE) }
  d <- cells; d$donor <- c("Oetjen2018_C", "Oetjen2018_Ck", "d2", "d2"); d$block <- c("Oetjen2018_C_block", "Oetjen2018_C_block", "d2", "d2")
  expect_equal(nso_calculation_metadata(d)$block, d$block)                    # R24 반복 표본 block 은 허용
  d$block[2] <- "Oetjen2018_Ck"; expect_error(nso_calculation_metadata(d), "block")
  d <- cells; d$donor <- NULL; expect_error(nso_calculation_metadata(d), "donor")
  d <- cells; d$cell_id[4] <- "c1"; expect_error(nso_calculation_metadata(d), "중복")
})


# ── Task 2 고정 예제: 6 세포 x 7 유전자 정수 count (library = 100, 50, 100, 100, 50, 100; 손으로 계산할 수 있게) ──
nso_t2_symbols <- c("MPO", "ELANE", "LTF", "CAMP", "MMP9", "FCGR3B", "OTHER")      # CAMP 는 입력에 있으나 전부 0, PRTN3 · LCN2 등은 입력에 없다
nso_t2_counts <- function() Matrix::Matrix(rbind(c(9, 1, 0, 0, 0, 0, 90), c(4, 0, 1, 0, 0, 0, 45), c(0, 0, 3, 0, 1, 0, 96),
                                                 c(0, 0, 9, 0, 4, 0, 87), c(0, 0, 0, 0, 2, 2, 46), c(1, 0, 0, 0, 0, 9, 90)), sparse = TRUE)
nso_t2_cells <- function() data.frame(cluster = c("A", "A", "A", "B", "B", "B"), block = c("b1", "b1", "b2", "b1", "b2", "b2"), study = "s1", stringsAsFactors = FALSE)
nso_t2_genes <- c("MPO", "ELANE", "PRTN3", "DEFA4", "BPI", "LTF", "CAMP", "LCN2", "MMP9", "FCGR3B", "FPR1", "CXCR2", "MKI67", "TOP2A", "CEBPE")
nso_t2_ge <- function(cells = nso_t2_cells(), X = nso_t2_counts())
  nso_gene_evidence(nrv_chunks_matrix(X, chunk = 4L), nrv_gene_map(nso_t2_symbols, nso_t2_genes), groups = list(pX = cells[, c("cluster", "block")]))
# 세포 edge (상삼각): A 내부 1–2, 2–3 · B 내부 4–5, 5–6 · A–B 2–4, 3–5 (w = 2), 3–6, 1–4
nso_t2_edges <- function() data.frame(i = c(1L, 2L, 4L, 5L, 2L, 3L, 3L, 1L), j = c(2L, 3L, 5L, 6L, 4L, 5L, 6L, 4L), w = c(1, 1, 1, 1, 1, 2, 1, 1))

test_that("gene evidence: 정규화 log1p(10^4 count / 전체 library) · 군집 x block 평균 · 검출률이 손 계산과 같고, 입력에 없는 유전자는 0 으로 채우지 않는다", {
  ge <- nso_t2_ge()
  expect_equal(ge$cell$library, c(100, 50, 100, 100, 50, 100))
  expect_equal(unname(ge$cell$log[1, "MPO"]), log1p(900)); expect_equal(unname(ge$cell$log[2, "MPO"]), log1p(800)); expect_equal(unname(ge$cell$log[6, "FCGR3B"]), log1p(900))
  expect_equal(unname(ge$cell$det[1, "ELANE"]), 1); expect_equal(unname(ge$cell$det[2, "ELANE"]), 0)
  expect_equal(ge$genes$status[ge$genes$gene == "PRTN3"], "missing_gene"); expect_equal(ge$genes$status[ge$genes$gene == "CAMP"], "in_input")
  expect_equal(ge$n_cells, 6L)
  g <- ge$groups; r <- function(cl, b, gn) g[g$partition == "pX" & g$cluster == cl & g$block == b & g$gene == gn, ]
  expect_equal(r("A", "b1", "MPO")$mean_log, (log1p(900) + log1p(800)) / 2); expect_equal(r("A", "b1", "MPO")$det_rate, 1)
  expect_equal(r("A", "b1", "ELANE")$det_rate, 0.5); expect_equal(r("A", "b1", "ELANE")$n_cells, 2L)
  expect_equal(r("B", "b2", "FCGR3B")$mean_log, (log1p(400) + log1p(900)) / 2)
  expect_equal(r("A", "b2", "CAMP")$n_detected, 0); expect_equal(r("A", "b2", "CAMP")$mean_log, 0)
  expect_false(any(g$gene %in% c("PRTN3", "LCN2")))
  X2 <- nso_t2_counts(); X2[1, 1] <- 9.5; expect_error(nso_t2_ge(X = X2), "raw count")
  expect_error(nso_gene_evidence(nrv_chunks_matrix(nso_t2_counts()), nrv_gene_map(nso_t2_symbols, nso_t2_genes), groups = list(pX = nso_t2_cells()[1:5, ])), "세포 수")
})

test_that("program evidence: missing_gene 과 관측 0 을 구별하고, 입력 유전자가 없는 프로그램은 0 점이 아니라 unavailable 이다", {
  ge <- nso_t2_ge(); pe <- nso_program_evidence(ge, nso_fixture_ref()); d <- pe$programs
  expect_equal(d$n_expected[d$program_id == "marker_Myeloblast"], 3L); expect_equal(d$n_observed[d$program_id == "marker_Myeloblast"], 2L)
  expect_equal(d$missing_genes[d$program_id == "marker_Myeloblast"], "PRTN3"); expect_equal(d$status[d$program_id == "marker_PM"], "no_genes_in_input")
  expect_equal(d$program_kind[d$program_id == "E06_core_down"], "transition_core_down"); expect_equal(d$from[d$program_id == "E06_core_down"], "PM")
  gp <- pe$groups; q <- function(cl, b, p) gp[gp$cluster == cl & gp$block == b & gp$program_id == p, ]
  expect_true(is.na(q("A", "b1", "marker_PM")$mean_log)); expect_equal(q("A", "b1", "marker_PM")$status, "unavailable_missing_genes")
  expect_equal(q("A", "b1", "marker_Myeloblast")$mean_log, ((log1p(900) + log1p(800)) / 2 + log1p(100) / 2) / 2)
  expect_equal(q("A", "b1", "marker_Myeloblast")$n_detected_genes, 2L); expect_equal(q("A", "b1", "marker_Myeloblast")$detected_genes, "MPO|ELANE")
  expect_equal(q("B", "b2", "marker_MY")$status, "observed_zero"); expect_equal(q("B", "b2", "marker_MY")$mean_log, 0)   # LTF · CAMP 모두 0 count
  expect_equal(q("B", "b2", "marker_MY")$detected_genes, "")
  expect_equal(q("A", "b2", "marker_MY")$mean_log, log1p(300) / 2)                    # LTF log1p(300), CAMP 관측 0, LCN2 입력 없음
  cm <- pe$cells$mean_log; cn <- pe$cells$n_detected
  expect_equal(unname(cm[1, "marker_Myeloblast"]), (log1p(900) + log1p(100)) / 2); expect_equal(unname(cn[1, "marker_Myeloblast"]), 2)
  expect_true(all(is.na(cm[, "marker_PM"]))); expect_equal(pe$cells$unavailable, rep("marker_PM:no_genes_in_input", 6))
  expect_equal(colnames(cm), paste0("marker_", NSO_STAGES))
})

test_that("program evidence: 공유 marker 는 두 단계 프로그램에 같은 값으로 남고, 유전자 1 개 프로그램은 따로 표시한다", {
  ge <- nso_t2_ge(); pe <- nso_program_evidence(ge, nso_fixture_ref(), min_cells = 1); d <- pe$programs; gp <- pe$groups
  q <- function(cl, b, p) gp[gp$cluster == cl & gp$block == b & gp$program_id == p, ]
  expect_equal(d$shared_genes[d$program_id == "marker_MY"], "LTF|LCN2"); expect_equal(d$shared_genes[d$program_id == "marker_MM"], "LTF|LCN2|MMP9")
  expect_true(d$single_gene_program[d$program_id == "E07_core_up"]); expect_false(d$single_gene_program[d$program_id == "E06_core_up"])
  expect_equal(q("B", "b1", "E07_core_up")$mean_log, log1p(400))                     # MMP9 하나
  expect_equal(q("B", "b1", "marker_MY")$mean_log, log1p(900) / 2)                   # LTF log1p(900), CAMP 0
  expect_equal(q("B", "b1", "marker_MM")$mean_log, (log1p(900) + log1p(400)) / 2)    # LTF, MMP9 (LCN2 입력 없음)
  cp <- pe$clusters; z <- cp[cp$partition == "pX" & cp$cluster == "A" & cp$program_id == "marker_Myeloblast", ]
  expect_equal(z$n_cells, 3L); expect_equal(z$n_blocks, 2L); expect_equal(z$n_blocks_min_cells, 2L)
  expect_equal(z$mean_log_pooled, ((log1p(900) + log1p(800)) / 3 + log1p(100) / 3) / 2)
  expect_equal(z$mean_log_equal_block, (((log1p(900) + log1p(800)) / 2 + log1p(100) / 2) / 2 + 0) / 2)
  z5 <- nso_program_evidence(ge, nso_fixture_ref())$clusters                          # 기본 min_cells = 5: 적격 block 없음 → equal-block 은 NA
  expect_true(is.na(z5$mean_log_equal_block[z5$cluster == "A" & z5$program_id == "marker_Myeloblast"]))
})

test_that("관측 연결: 군집 쌍 · 경계 세포 수 · 경계 평균 · 같은 block 대비가 손 계산과 같고, 상태(빈 경계 · 적격 block 없음 · 한 block)를 남긴다", {
  ge <- nso_t2_ge(); pe <- nso_program_evidence(ge, nso_fixture_ref())
  oc <- nso_observed_connections(nso_t2_edges(), nso_t2_cells(), ge, programs = pe$programs, min_cells = 1)
  expect_equal(nrow(oc$pairs), 1L); expect_equal(oc$pairs$cluster_a, "A"); expect_equal(oc$pairs$cluster_b, "B")
  expect_equal(oc$pairs$n_cell_edges, 4L); expect_equal(oc$pairs$weight_sum, 5); expect_equal(oc$pairs$n_boundary_a, 3L); expect_equal(oc$pairs$n_boundary_b, 3L)
  f <- oc$features; L <- f[f$feature == "LTF", ]
  expect_equal(L$mean_boundary_a, (log1p(200) + log1p(300)) / 3); expect_equal(L$mean_boundary_b, log1p(900) / 3)
  expect_equal(L$delta_pooled, log1p(900) / 3 - (log1p(200) + log1p(300)) / 3)
  expect_equal(L$delta_equal_block, ((log1p(900) - log1p(200) / 2) + (0 - log1p(300))) / 2); expect_equal(L$n_eligible_blocks, 2L); expect_equal(L$status, "ok")
  expect_equal(f$status[f$feature == "PRTN3"], "missing_gene"); expect_true(is.na(f$delta_pooled[f$feature == "PRTN3"]))
  bn <- f[f$feature == "marker_BN", ]
  d_mmp9 <- (2 * log1p(400)) / 3 - log1p(100) / 3; d_fcgr <- (log1p(400) + log1p(900)) / 3
  expect_equal(bn$feature_kind, "program"); expect_equal(bn$delta_pooled, (d_mmp9 + d_fcgr) / 2)
  expect_equal(f$status[f$feature == "marker_PM"], "unavailable_missing_genes"); expect_true(is.na(f$delta_pooled[f$feature == "marker_PM"]))
  o5 <- nso_observed_connections(nso_t2_edges(), nso_t2_cells(), ge, programs = pe$programs)          # 기본 min_cells = 5
  expect_equal(o5$features$status[o5$features$feature == "LTF"], "no_eligible_block"); expect_true(is.na(o5$features$delta_equal_block[o5$features$feature == "LTF"]))
  c1 <- nso_t2_cells(); c1$block <- "b1"                                               # 한 공여자: block 하나
  o1 <- nso_observed_connections(nso_t2_edges(), c1, nso_t2_ge(cells = c1), programs = pe$programs, min_cells = 1); L1 <- o1$features[o1$features$feature == "LTF", ]
  expect_equal(L1$status, "single_block"); expect_equal(L1$n_eligible_blocks, 1L); expect_equal(L1$delta_equal_block, L1$delta_pooled)
  c2 <- nso_t2_cells(); c2$block <- c("b1", "b1", "b1", "b2", "b2", "b2")              # A–B edge 가 모두 다른 block 사이: 같은 block 경계가 비었다
  o2 <- nso_observed_connections(nso_t2_edges(), c2, nso_t2_ge(cells = c2), programs = pe$programs, min_cells = 1); L2 <- o2$features[o2$features$feature == "LTF", ]
  expect_equal(L2$status, "no_same_block_boundary"); expect_equal(L2$n_blocks_with_boundary, 0L); expect_true(is.na(L2$delta_equal_block))
  expect_equal(L2$delta_pooled, L$delta_pooled)                                        # 전체 경계 대비는 block 과 무관
})

test_that("관측 연결: 쌍이 없으면 0 행 표를 오류 없이 돌려주고, 대비 방향은 저장 순서의 부호 반전일 뿐이다", {
  ge <- nso_t2_ge(); pe <- nso_program_evidence(ge, nso_fixture_ref()); e <- nso_t2_edges(); cl <- nso_t2_cells()
  e_in <- e[cl$cluster[e$i] == cl$cluster[e$j], ]
  o0 <- nso_observed_connections(e_in, cl, ge, programs = pe$programs, min_cells = 1)
  expect_equal(nrow(o0$pairs), 0L); expect_equal(nrow(o0$features), 0L); expect_true(all(c("cluster_a", "feature", "delta_pooled", "status") %in% names(o0$features)))
  o <- nso_observed_connections(e, cl, ge, programs = pe$programs, min_cells = 1)
  cz <- cl; cz$cluster[cz$cluster == "A"] <- "Z"                                       # 이름만 바꿔 저장 순서를 뒤집는다 (B < Z)
  oz <- nso_observed_connections(e, cz, nso_t2_ge(cells = cz), programs = pe$programs, min_cells = 1)
  expect_equal(oz$pairs$cluster_a, "B"); expect_equal(oz$features$delta_pooled[oz$features$feature == "LTF"], -o$features$delta_pooled[o$features$feature == "LTF"])
  expect_equal(oz$features$delta_equal_block[oz$features$feature == "LTF"], -o$features$delta_equal_block[o$features$feature == "LTF"])
  expect_true(nzchar(o$meta$orientation)); expect_equal(o$meta$min_cells, 1)
})

test_that("B 내부 연속성: 모든 관측 삼중쌍 A–B–C 의 상태가 R24 nrv_internal_coordinate 와 같다 (그래프 연결 ≠ 내부 연속)", {
  cl <- c("A", "B", "B", "B", "B", "C", "D")                                         # 1 = a1, 2–5 = B, 6 = c1, 7 = d1
  base <- data.frame(i = c(1L, 2L, 4L, 5L, 4L), j = c(2L, 3L, 5L, 6L, 7L), w = 1)     # B = {2,3} | {4,5} 두 성분; a1–b1, b4–c1, b3–d1
  variants <- list(disconnected = base, bridged = rbind(base, data.frame(i = 3L, j = 4L, w = 1)),
                   nested = data.frame(i = c(1L, 2L, 2L, 4L, 4L), j = c(2L, 3L, 6L, 5L, 7L), w = 1))
  want <- list(disconnected = c("A|C" = "internal_disconnected", "A|D" = "internal_disconnected", "C|D" = "span_ok"),
               bridged = c("A|C" = "span_ok", "A|D" = "span_ok", "C|D" = "span_ok"),
               nested = c("A|C" = "no_internal_span", "A|D" = "internal_disconnected", "C|D" = "internal_disconnected"))
  for (v in names(variants)) { e <- variants[[v]]; ic <- nso_internal_continuity(e, cl); tb <- ic$triples[ic$triples$cluster_B == "B", ]
    got <- stats::setNames(tb$status, paste(tb$cluster_A, tb$cluster_C, sep = "|"))
    expect_equal(got[names(want[[v]])], want[[v]], info = v)
    for (k in seq_len(nrow(tb))) expect_equal(tb$status[k], nrv_internal_coordinate(e, cl, "B", tb$cluster_A[k], tb$cluster_C[k])$status, info = paste(v, k)) }
  ic <- nso_internal_continuity(variants$disconnected, cl); s <- ic$clusters
  expect_equal(s$n_components[s$cluster == "B"], 2L); expect_equal(s$n_internal_disconnected[s$cluster == "B"], 2L); expect_equal(s$n_neighbors[s$cluster == "B"], 3L)
  expect_equal(nrow(ic$triples[ic$triples$cluster_B == "A", ]), 0L)                    # 이웃이 하나뿐인 군집에는 삼중쌍이 없다
})

test_that("gene evidence: keep_cells 로 고른 세포의 원 count 와 library 를 그대로 남긴다 (실제 RNA spot-check 용)", {
  ge <- nso_gene_evidence(nrv_chunks_matrix(nso_t2_counts(), chunk = 4L), nrv_gene_map(nso_t2_symbols, nso_t2_genes),
                          groups = list(pX = nso_t2_cells()[, c("cluster", "block")]), keep_cells = c(5L, 1L))
  k <- ge$kept
  expect_equal(k$cell, c(1L, 5L)); expect_equal(k$library, c(100, 50))
  expect_equal(unname(k$counts[, "MPO"]), c(9, 0)); expect_equal(unname(k$counts[, "MMP9"]), c(0, 2)); expect_equal(unname(k$counts[, "FCGR3B"]), c(0, 2))
  expect_false("PRTN3" %in% colnames(k$counts))
  expect_null(nso_t2_ge()$kept$counts)                                                  # 고르지 않으면 남기지 않는다
})


# ── Task 3: 실행기 · 그림 (소형 합성 자료로 처음부터 끝까지; 실제 자료는 같은 실행기를 attempt_01 에) ──
#' 소형 세계: reference 고정 예제, O 의 NRV structure 단계(manifest 포함), R57 · R54 비교 결과(posthoc 저자 열 포함 — 읽히면 안 된다),
#' M 의 graph / graph_ids / h5ad 자리 파일, 가짜 reader (author() · umap() 은 calc_manifest.json 이 없으면 멈춘다)
nso_t3_world <- function(edges = nso_t2_edges(), umap_rows = 6L) {
  w <- tempfile("nso_world"); dir.create(w); out <- file.path(w, "out", "attempt_01"); sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
  ref <- nso_write_fixture(); ids <- paste0("cell", 1:6); donor <- c("d1", "d1", "d2", "d1", "d2", "d2")
  mk <- function(cl) data.frame(cell = 1:6, cell_id = ids, cluster = cl, block = donor, study = "s1", sorting = "Bulk", stringsAsFactors = FALSE)
  cells <- list(p92 = mk(c("A", "A", "A", "B", "B", "B")), p51 = mk(rep("X", 6)))
  M <- file.path(w, "M"); for (rel in c(NSO_GRAPH, NSO_H5AD)) dir.create(file.path(M, dirname(rel)), recursive = TRUE, showWarnings = FALSE)
  writeLines("graph-bytes", file.path(M, NSO_GRAPH)); writeLines(ids, file.path(M, NSO_GRAPH_IDS)); writeLines("h5-bytes", file.path(M, NSO_H5AD))
  O <- file.path(w, "O")
  for (p in names(cells)) { att <- file.path(O, NSO_NRV[[p]]); sd <- nrv_stage_dir(att, "structure"); dir.create(sd, recursive = TRUE); cl <- cells[[p]]
    saveRDS(list(cells = cl, bm = data.frame(block = c("d1", "d2"), study = "s1", sorting_class = "bulk", n_cells = c(3L, 3L)), folds = data.frame(), keys = data.frame()),
            file.path(sd, "cells_lossless.rds"))
    saveRDS(edges, file.path(sd, "edges.rds")); saveRDS(list(full = nrv_pair_table(edges, cl$cluster, cl$block, cl$study), folds = list()), file.path(sd, "pairs.rds"))
    mdu_write(data.frame(input = c("graph", "graph_ids", "h5ad"), path = file.path(M, c(NSO_GRAPH, NSO_GRAPH_IDS, NSO_H5AD)),
                         sha256 = c(sha(file.path(M, NSO_GRAPH)), sha(file.path(M, NSO_GRAPH_IDS)), sha(file.path(M, NSO_H5AD)))), file.path(sd, "input_sha256.tsv"))
    nrv_write_stage_manifest(att, "structure", files = file.path(sd, c("cells_lossless.rds", "edges.rds", "pairs.rds", "input_sha256.tsv")), upstream = "test") }
  cma <- file.path(w, "cma"); fdt <- file.path(w, "fdt"); dir.create(file.path(cma, "p92"), recursive = TRUE); dir.create(file.path(cma, "p51")); dir.create(fdt)
  ct <- function(p, cl) data.frame(partition = p, cluster = cl, n_cells = 3, status = "annotated", predicted_type = paste0("T_", cl), predicted_M = 0.5, warnings = "",
                                   annotation_kind = "marker_only_provisional", posthoc_author_top_label = "SHOULD_NOT_BE_READ", posthoc_author_top_frac = 1, n_author_labels = 1L,
                                   stringsAsFactors = FALSE)
  mdu_write(ct("p92", c("A", "B")), file.path(cma, "p92", "cluster_marker_annotation.tsv")); mdu_write(ct("p51", "X"), file.path(cma, "p51", "cluster_marker_annotation.tsv"))
  writeLines(jsonlite::toJSON(list(files = list("p92/cluster_marker_annotation.tsv" = sha(file.path(cma, "p92", "cluster_marker_annotation.tsv")),
                                                "p51/cluster_marker_annotation.tsv" = sha(file.path(cma, "p51", "cluster_marker_annotation.tsv")))), auto_unbox = TRUE),
             file.path(cma, "manifest.json"))
  Mx <- function(cl) matrix(seq_len(length(cl) * length(NSO_COMPETING)) / 10, length(cl), dimnames = list(cl, NSO_COMPETING))
  saveRDS(list(scores = list(p92 = list(n = c(A = 3, B = 3), M = Mx(c("A", "B"))), p51 = list(n = c(X = 6), M = Mx("X"))),
               author_composition = data.frame(author_label = "SHOULD_NOT_BE_READ")), file.path(fdt, "lossless.rds"))
  writeLines(jsonlite::toJSON(list(files = list(lossless.rds = sha(file.path(fdt, "lossless.rds")))), auto_unbox = TRUE), file.path(fdt, "manifest.json"))
  U <- cbind(c(0, 1, 2, 5, 6, 7), c(0, 1, 0, 5, 6, 5)); au <- c("GMP-Neut", "GMP-Neut", "Early GMP", "CD14 Mono", NA, "CD14 Mono")
  guard <- function(what) if (!file.exists(file.path(out, "calc_manifest.json"))) stop(what, " 를 결과 고정 전에 읽었다")
  reader <- list(ids = function() ids, obs = function(name) switch(name, cell_id = ids, donor_id = donor, Study = rep("s1", 6), stop("모르는 obs: ", name)),
                 symbols = function() nso_t2_symbols, raw_chunks = function(chunk) nrv_chunks_matrix(nso_t2_counts(), chunk),
                 umap = function() { guard("UMAP"); U[seq_len(umap_rows), , drop = FALSE] }, author = function() { guard("저자 라벨"); au },
                 input_files = function() c(graph = file.path(M, NSO_GRAPH), graph_ids = file.path(M, NSO_GRAPH_IDS), h5ad = file.path(M, NSO_H5AD)))
  pins <- list(reference = c(nodes = sha(ref[["nodes"]]), edges = sha(ref[["edges"]])), graph = sha(file.path(M, NSO_GRAPH)), h5ad = sha(file.path(M, NSO_H5AD)))
  list(argv = c("--data-root", M, "--o-root", O, "--reference-dir", dirname(ref[["nodes"]]), "--out-dir", out, "--cma-dir", cma, "--fdt-dir", fdt, "--chunk", "4"),
       out = out, reader = reader, pins = pins, ids = ids)
}
nso_t3_run <- function(w, pins = w$pins, argv = w$argv) { code <- NA_integer_; utils::capture.output(code <- nso_main(argv, root = ROOT, reader = w$reader, pins = pins)); code }
nso_read <- function(out, f) utils::read.delim(if (grepl("\\.gz$", f)) gzfile(file.path(out, f)) else file.path(out, f), colClasses = "character", quote = "", check.names = FALSE)
nso_lines <- function(out, f) readLines(if (grepl("\\.gz$", f)) gzfile(file.path(out, f)) else file.path(out, f))

test_that("실행기 (소형 합성 자료): 끝까지 돌고, 전체 cell ID · 분모 · 관측 쌍을 보존하며, 예측 단계 라벨 · 저자 열을 계산 표에 만들지 않는다", {
  w <- nso_t3_world(); expect_equal(nso_t3_run(w), 0L); out <- w$out
  need <- c("stage_reference_availability.tsv", "cell_stage_evidence.tsv.gz", "cluster_block_programs.tsv.gz", "cluster_block_genes.tsv.gz", "cluster_programs.tsv",
            "observed_pairs.tsv", "observed_edge_evidence.tsv.gz", "internal_continuity.tsv.gz", "internal_cluster_summary.tsv", "stage_mapping_limitations.tsv",
            "cluster_reference_comparison.tsv", "spot_check_cells.tsv", "figure_text.tsv", "posthoc_author_cluster_composition.tsv", "posthoc_author_program_summary.tsv",
            "input_sha256.tsv", "source_snapshot.tsv", "environment.tsv", "checks.tsv", "lossless.rds", "calc_manifest.json", "manifest.json", "README_ko.md",
            file.path("figures", c("fig01_neut_stage_evidence_umap.png", "fig02_observed_connections_umap.png")))
  expect_true(all(file.exists(file.path(out, need))), info = paste(need[!file.exists(file.path(out, need))], collapse = ", "))
  ce <- nso_read(out, "cell_stage_evidence.tsv.gz")
  expect_equal(ce$cell_id, w$ids); expect_equal(ce$cluster92, c("A", "A", "A", "B", "B", "B")); expect_equal(ce$r57_type_p92, c("T_A", "T_A", "T_A", "T_B", "T_B", "T_B"))
  expect_false(any(grepl("stage_label|predicted_stage|author", names(ce)))); expect_equal(as.numeric(ce$log_MPO[1]), log1p(900))
  expect_equal(as.numeric(ce$marker_Myeloblast_mean_log[1]), (log1p(900) + log1p(100)) / 2); expect_equal(unique(ce$unavailable_reason), "marker_PM:no_genes_in_input")
  for (f in c("cell_stage_evidence.tsv.gz", "cluster_reference_comparison.tsv", "cluster_programs.tsv", "lossless.rds")) if (!grepl("rds$", f)) expect_false(any(grepl("SHOULD_NOT_BE_READ", nso_lines(out, f))), info = f)
  expect_false(any(grepl("SHOULD_NOT_BE_READ", unlist(readRDS(file.path(out, "lossless.rds"))$comparison))))
  cb <- nso_read(out, "cluster_block_programs.tsv.gz"); s <- tapply(as.numeric(cb$n_cells), paste(cb$partition, cb$program_id), sum)
  expect_true(all(s == 6)); expect_equal(length(s), 2L * length(unique(cb$program_id)))   # 분할 x 프로그램마다 전체 분모 6
  oe <- nso_read(out, "observed_edge_evidence.tsv.gz")
  expect_equal(unique(oe$partition), "p92"); expect_equal(unique(paste(oe$cluster_a, oe$cluster_b)), "A B")   # p51 은 군집 하나: 관측 쌍 0
  expect_equal(nso_read(out, "observed_pairs.tsv")$partition, "p92")
  ck <- nso_read(out, "checks.tsv"); expect_true(all(ck$pass == "TRUE"), info = paste(ck$check[ck$pass != "TRUE"], collapse = ", "))
  expect_true(all(c("inputs_pinned", "pairs_match_stored_nrv_p92", "pairs_match_stored_nrv_p51", "tables_roundtrip_ok", "author_and_umap_read_after_fix") %in% ck$check))
  man <- jsonlite::fromJSON(file.path(out, "manifest.json"))$files
  expect_true(all(vapply(names(man), function(f) identical(digest::digest(file.path(out, f), algo = "sha256", file = TRUE), man[[f]]), TRUE)))
  sp <- nso_read(out, "spot_check_cells.tsv"); expect_equal(sp$cell_id, c("cell1", "cell4")); expect_equal(as.numeric(sp$count_MPO), c(9, 0))
  expect_equal(as.numeric(sp$library_size), c(100, 100)); expect_equal(as.numeric(sp$log_LTF[2]), log1p(900))
  pa <- nso_read(out, "posthoc_author_cluster_composition.tsv"); a <- pa[pa$partition == "p92" & pa$cluster == "A", ]
  expect_equal(as.numeric(a$n_cells[a$author_label == "GMP-Neut"]), 2); expect_true("(no author label)" %in% pa$author_label)
})

test_that("표시: 원 UMAP 은 세포 수를 대조하고, 미검출 프로그램과 관측 쌍 없음도 그림에 이유를 적는다 (글자 'NA' 없음)", {
  w <- nso_t3_world(edges = nso_t2_edges()[1:4, ]); expect_equal(nso_t3_run(w), 0L)       # 군집 안 edge 만: 관측 쌍 없음
  ft <- nso_read(w$out, "figure_text.tsv")
  expect_true(any(ft$element == "unavailable" & ft$panel == "PM"))                         # DEFA4 · BPI 가 입력에 없다
  expect_true(any(ft$element == "no_observed_pairs")); expect_false(any(ft$text %in% c("NA", "NaN", "Inf")))
  ttl <- ft$text[ft$element == "title" & ft$panel %in% NSO_STAGES]                           # 패널 제목에 단계 이름을 두 번 쓰지 않는다 (full_name 이 단계 이름으로 시작할 때)
  expect_equal(length(ttl), 6L); expect_false(any(grepl("^([^ ]+) · \\1", ttl)))
  expect_true(all(file.info(file.path(w$out, "figures", c("fig01_neut_stage_evidence_umap.png", "fig02_observed_connections_umap.png")))$size > 0))
  w2 <- nso_t3_world(umap_rows = 5L); expect_equal(nso_t3_run(w2), 4L)                      # UMAP 세포 수가 다르면 멈춘다
  expect_true(any(grepl("UMAP", readLines(file.path(w2$out, "run.log")))))
})

test_that("결과 고정 전 저자 라벨 · UMAP 입력 금지, 결과 폴더 덮어쓰기 금지, 고정 hash 불일치 · 인자 누락 중지", {
  w <- nso_t3_world(); expect_equal(nso_t3_run(w), 0L)         # reader 의 author() · umap() 은 calc_manifest.json 이 없으면 멈춘다 → 통과 = 고정 뒤에 읽음
  cm <- jsonlite::fromJSON(file.path(w$out, "calc_manifest.json"))$files
  expect_true(all(c("cell_stage_evidence.tsv.gz", "observed_edge_evidence.tsv.gz", "stage_reference_availability.tsv") %in% names(cm)))
  expect_false(any(grepl("posthoc|figure", names(cm))))
  expect_equal(nso_t3_run(w), 6L)                                                          # 같은 폴더에 다시 쓰지 않는다
  w3 <- nso_t3_world(); bad <- w3$pins; bad$h5ad <- strrep("0", 64); expect_equal(nso_t3_run(w3, pins = bad), 4L)
  expect_true(any(grepl("hash", readLines(file.path(w3$out, "run.log")))))
  w4 <- nso_t3_world(); expect_equal(nso_t3_run(w4, argv = w4$argv[1:6]), 2L)               # --out-dir 없음
})
