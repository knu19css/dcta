# =============================================================================
# test_nso_connectivity_audit.R — 단계 marker 검출 세포의 원 그래프 연결 감사 (R62) 테스트
#   계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_connectivity_instruction.md
#   검출 집합(1 개 이상 · 2 개 이상) · 유도 부분그래프 성분 · 인접 쌍 edge(겹침 · A-only↔B-only · 같은 block) ·
#   순서 도달: 쌍마다 edge 는 있어도 중간 단계 성분이 끊기면 불가 / 전체 가능 / 겹침만 있고 edge 가 없으면 zero-hop 불가 /
#   cross-block edge 로만 가능 / 무관한 세포를 거친 우회 불가 / 실제 edge 가 있으면 같은 세포가 다음 층에 다시 나올 수 있음
# =============================================================================
library(testthat)
nca_test_root <- function() {
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}
ROOT <- nca_test_root()
for (f in c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R", "27_run_nrv.R",
            "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R", "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R",
            "58_nso_inputs.R", "62_nso_connectivity_audit.R"))
  source(file.path(ROOT, "R", "branch_biology", f))

test_that("검출 집합: marker log > 0 인 유전자 수가 기준 이상인 세포 (1 개 이상 · 2 개 이상), 단계 겹침은 그대로 두고 빠진 marker 는 멈춘다", {
  L <- Matrix::sparseMatrix(i = c(1, 1, 2, 3, 3, 3), j = c(1, 2, 2, 1, 3, 4), x = c(1.2, 0.5, 2.0, 0.7, 0.1, 0.9), dims = c(4, 4),
                            dimnames = list(NULL, c("G1", "G2", "G3", "G4")))
  mk <- list(S1 = c("G1", "G2"), S2 = c("G3", "G4"))
  s1 <- nca_detection_sets(L, mk, 1L); s2 <- nca_detection_sets(L, mk, 2L)
  expect_equal(unname(s1[, "S1"]), c(TRUE, TRUE, TRUE, FALSE)); expect_equal(unname(s1[, "S2"]), c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(unname(s2[, "S1"]), c(TRUE, FALSE, FALSE, FALSE)); expect_equal(unname(s2[, "S2"]), c(FALSE, FALSE, TRUE, FALSE))
  expect_true(s1[3, "S1"] && s1[3, "S2"])                                                   # 한 세포가 두 단계에 동시에 속해도 그대로
  expect_error(nca_detection_sets(L, list(S1 = c("G1", "G9")), 1L), "G9")
})

# 4 세포 예제: 1 = A only, 2 = A·B, 3 = B only, 4 = B only (block b2). edge 1–2, 1–3, 2–3, 3–4, 1–4
nca_ex_edges <- function() data.frame(i = c(1L, 1L, 2L, 3L, 1L), j = c(2L, 3L, 3L, 4L, 4L), w = 1)
nca_ex_block <- c("b1", "b1", "b1", "b2")

test_that("인접 쌍 edge: A–B edge 수 · 양끝 고유 세포 · 같은 block · 겹침 세포 · A-only↔B-only edge 가 손 계산과 같다", {
  A <- c(TRUE, TRUE, FALSE, FALSE); B <- c(FALSE, TRUE, TRUE, TRUE)
  p <- nca_pair_edges(nca_ex_edges(), A, B, nca_ex_block)
  expect_equal(p$n_overlap_cells, 1L); expect_equal(p$n_edges, 4L)                          # (1,2) (1,3) (2,3) (1,4); (3,4) 는 B–B 라 아님
  expect_equal(p$n_endpoints_A, 2L); expect_equal(p$n_endpoints_B, 3L); expect_equal(p$n_same_block_edges, 3L)
  expect_equal(p$n_edges_Aonly_Bonly, 2L); expect_equal(p$n_edges_with_overlap_endpoint, 2L)
  ps <- nca_pair_edges(nca_ex_edges(), A, B, nca_ex_block, same_block = TRUE)               # 같은 block edge 만: (1,4) 가 빠진다
  expect_equal(ps$n_edges, 3L); expect_equal(ps$n_edges_Aonly_Bonly, 1L)
})

test_that("연결 성분: 전체 그래프와 단계 유도 부분그래프를 구별하고, 같은 block 제한은 단계 안 edge 에도 적용한다", {
  w <- nca_components(4L, nca_ex_edges()); expect_equal(w$n_components, 1L); expect_equal(w$largest, 4L)
  B <- c(FALSE, TRUE, TRUE, TRUE)
  cb <- nca_components(4L, nca_ex_edges(), keep = B); expect_equal(cb$n_components, 1L); expect_equal(cb$largest, 3L); expect_equal(cb$n_isolated, 0L)
  expect_true(is.na(cb$membership[1]))                                                      # 집합 밖 세포는 성분이 없다
  cs <- nca_components(4L, nca_ex_edges(), keep = B, block = nca_ex_block, same_block = TRUE)
  expect_equal(cs$n_components, 2L); expect_equal(cs$largest, 2L); expect_equal(cs$n_isolated, 1L)   # (3,4) 가 block 사이라 빠짐 → 세포 4 고립
})

test_that("순서 도달 (지시서 6 번 기대 예제): 쌍마다 edge 가 있어도 중간 단계 성분이 끊기면 전체 순서는 불가, B 안 edge 가 있으면 가능", {
  S <- cbind(A = c(TRUE, FALSE, FALSE, FALSE), B = c(FALSE, TRUE, TRUE, FALSE), C = c(FALSE, FALSE, FALSE, TRUE))
  e1 <- data.frame(i = c(1L, 3L), j = c(2L, 4L))                                            # A–B: 1–2, B–C: 3–4, B 안 edge 없음
  expect_equal(nca_pair_edges(e1, S[, "A"], S[, "B"], rep("b", 4))$n_edges, 1L); expect_equal(nca_pair_edges(e1, S[, "B"], S[, "C"], rep("b", 4))$n_edges, 1L)
  r <- nca_ordered_reach(4L, e1, S)
  expect_equal(r$table$n_components, c(1L, 2L, 1L)); expect_equal(r$table$n_reachable_cells, c(1L, 1L, 0L)); expect_equal(r$first_break, "B>C")
  r2 <- nca_ordered_reach(4L, rbind(e1, data.frame(i = 2L, j = 3L)), S)
  expect_equal(r2$table$n_reachable_cells, c(1L, 2L, 1L)); expect_equal(r2$first_break, ""); expect_true(all(r2$table$reached))
})

test_that("순서 도달: 겹침만 있고 원 edge 가 없으면 zero-hop 전진 불가, 무관한 세포를 거친 우회 불가, 실제 edge 가 있으면 같은 세포가 다음 층에 다시 나올 수 있다", {
  S3 <- cbind(A = c(TRUE, FALSE), B = c(TRUE, FALSE), C = c(FALSE, TRUE))                    # 세포 1 = A · B, 세포 2 = C, edge 1–2 만
  r3 <- nca_ordered_reach(2L, data.frame(i = 1L, j = 2L), S3)
  expect_equal(r3$table$n_reachable_cells, c(1L, 0L, 0L)); expect_equal(r3$first_break, "A>B")
  S5 <- cbind(A = c(TRUE, FALSE, FALSE), B = c(FALSE, FALSE, TRUE))                         # 세포 2 는 어느 단계도 아님: 1–2–3 우회
  r5 <- nca_ordered_reach(3L, data.frame(i = c(1L, 2L), j = c(2L, 3L)), S5)
  expect_equal(r5$table$n_reachable_cells, c(1L, 0L)); expect_equal(r5$first_break, "A>B")
  S6 <- cbind(A = c(TRUE, FALSE), B = c(TRUE, TRUE))                                        # 세포 1 = A · B, 세포 2 = B, edge 1–2
  r6 <- nca_ordered_reach(2L, data.frame(i = 1L, j = 2L), S6)
  expect_equal(r6$table$n_reachable_cells, c(1L, 2L)); expect_equal(r6$table$n_reappearing_cells, c(0L, 1L))   # 세포 1 은 실제 edge 1–2 로 B 층에 다시 나온다
})

test_that("순서 도달: cross-block edge 로만 이어지는 순서는 모든 edge 에서는 가능하고 같은 block edge 만으로는 끊긴다", {
  S4 <- cbind(A = c(TRUE, FALSE, FALSE), B = c(FALSE, TRUE, FALSE), C = c(FALSE, FALSE, TRUE)); e4 <- data.frame(i = c(1L, 2L), j = c(2L, 3L)); b4 <- c("b1", "b2", "b2")
  ra <- nca_ordered_reach(3L, e4, S4, block = b4); expect_equal(ra$table$n_reachable_cells, c(1L, 1L, 1L)); expect_equal(ra$first_break, "")
  rs <- nca_ordered_reach(3L, e4, S4, block = b4, same_block = TRUE); expect_equal(rs$table$n_reachable_cells, c(1L, 0L, 0L)); expect_equal(rs$first_break, "A>B")
})

test_that("그림 축: 막대 좌표가 log10(n + 1) 이므로 눈금은 원 세포 수 0, 1, 10, 100 … 을 실제 좌표 log10(값 + 1) 에 둔다 (0 이 1 자리에 보이지 않게)", {
  tk <- .nca_log_ticks(log10(60178) * 1.18)
  expect_equal(tk$value, c(0, 1, 10, 100, 1000, 10000, 100000)); expect_equal(tk$at, log10(c(1, 2, 11, 101, 1001, 10001, 100001)))
  expect_equal(tk$label, c("0", "1", "10", "100", "1,000", "10,000", "100,000"))
  expect_equal(.nca_log_ticks(log10(5) * 1.18)$value, c(0, 1))                                # 작은 축 (최대 4 세포, 위끝 0.825): 0 과 1 만, 10 (좌표 1.04) 은 축 밖
})

test_that("실행기: 인자 확인, 결과 폴더가 비어 있지 않으면 덮어쓰지 않는다", {
  expect_equal(nca_main(c("--out-dir")), 2L)
  d <- tempfile("nca_out"); dir.create(d); writeLines("x", file.path(d, "keep.txt"))
  expect_equal(nca_main(c("--in-dir", "x", "--o-root", "y", "--out-dir", d)), 6L)
})
