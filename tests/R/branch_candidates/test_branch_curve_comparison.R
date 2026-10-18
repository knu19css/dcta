# =============================================================================
# test_branch_curve_comparison.R — 단일 곡선 대 분기 곡선 비교 (작은 합성자료만 쓴다)
# 실행: Rscript tests/R/branch_candidates/run_tests.R   (저장소 루트 기준 경로만 쓴다)
# =============================================================================

root <- Sys.getenv("RBC_REPO_ROOT", ".")
source(file.path(root, "R", "branch_candidates", "branch_curve_comparison.R"))
source(file.path(root, "R", "branch_candidates", "simulate_branch_data.R"))

small <- function(...) rbc_simulate(cells_per_donor = 90L, n_genes = 6L, ...)
run <- function(d, ...) compare_branch_curves(d$counts, d$cell_data, d$truth$split_point, ...)

# ── 입력 검증 ────────────────────────────────────────────────────────────────
test_that("count columns and cell_data rows must name the same cells in the same order", {
  d <- small(scenario = "two_arms")
  bad <- d; colnames(bad$counts)[1] <- "not_a_cell"
  expect_error(run(bad), "cell_id")
  bad <- d; bad$counts <- bad$counts[, c(2, 1, 3:ncol(bad$counts))]
  expect_error(run(bad), "cell_id")
  bad <- d; rownames(bad$counts)[2] <- rownames(bad$counts)[1]
  expect_error(run(bad), "gene names")
})

test_that("counts must be non-negative integers", {
  d <- small(scenario = "two_arms")
  bad <- d; bad$counts[1, 1] <- -1
  expect_error(run(bad), "non-negative integer")
  bad <- d; storage.mode(bad$counts) <- "double"; bad$counts[1, 1] <- 1.5
  expect_error(run(bad), "non-negative integer")
  bad <- d; bad$counts[2, 3] <- NA
  expect_error(run(bad), "non-negative integer")
})

test_that("donor, arm, position and library size are required and valid", {
  d <- small(scenario = "two_arms")
  bad <- d; bad$cell_data$donor_id[3] <- NA; expect_error(run(bad), "donor_id")
  bad <- d; bad$cell_data$arm[4] <- ""; expect_error(run(bad), "arm")
  bad <- d; bad$cell_data$library_size[5] <- 0; expect_error(run(bad), "library_size")
  bad <- d; bad$cell_data$position[6] <- 1.2; expect_error(run(bad), "position")
  bad <- d; bad$cell_data$library_size <- NULL; expect_error(run(bad), "library_size")
})

test_that("split_point must lie strictly inside the observed positions", {
  d <- small(scenario = "two_arms")
  for (sp in list(NA_real_, 0, 1, -0.1, 1.3, c(0.4, 0.5), max(d$cell_data$position)))
    expect_error(compare_branch_curves(d$counts, d$cell_data, sp), "split_point")
})

test_that("trunk must precede the split and at least two arms must follow it", {
  d <- small(scenario = "two_arms")
  bad <- d; i <- which(bad$cell_data$arm != "trunk")[1]; bad$cell_data$arm[i] <- "trunk"
  expect_error(run(bad), "trunk cells after split_point")
  bad <- d; j <- which(bad$cell_data$arm == "trunk")[1]; bad$cell_data$arm[j] <- "A"
  expect_error(run(bad), "arm cells before split_point")
  bad <- d; bad$cell_data$arm[bad$cell_data$arm != "trunk"] <- "A"
  expect_error(run(bad), "at least two")
})

test_that("unsupported options are rejected instead of silently reinterpreted", {
  d <- small(scenario = "two_arms")
  expect_error(run(d, geometry_source = "per_fold"), "geometry_source")
  expect_error(run(d, gene_provenance = "independent"), "gene_provenance")
  expect_error(run(d, spline_df = 0L), "spline_df")
  expect_error(run(d, folds = rep(1L, 3)), "folds")
})

# ── 모형 공간 ────────────────────────────────────────────────────────────────
test_that("branch terms are symmetric sum-to-zero deviations that vanish on the trunk", {
  C <- rbc_arm_contrasts(c("B", "C", "A"))
  expect_identical(dim(C), c(3L, 2L))
  expect_equal(unname(colSums(C)), c(0, 0))
  expect_identical(rownames(C), c("A", "B", "C"))
  d <- small(scenario = "two_arms", n_arms = 3L)
  cd <- d$cell_data; sp <- d$truth$split_point
  basis <- rbc_spline_basis(cd$position, 3L)
  Cm <- rbc_arm_contrasts(setdiff(unique(cd$arm), "trunk"))
  Xs <- rbc_design(cd$position, cd$arm, sp, basis, Cm, "single")
  Xb <- rbc_design(cd$position, cd$arm, sp, basis, Cm, "branch")
  expect_identical(ncol(Xs), 1L + 3L)
  expect_identical(ncol(Xb) - ncol(Xs), 2L * (3L - 1L))
  extra <- Xb[, -seq_len(ncol(Xs)), drop = FALSE]
  expect_true(all(extra[cd$arm == "trunk", ] == 0))
  expect_true(all(extra[cd$position <= sp, ] == 0))
  expect_equal(unname(Xb[, seq_len(ncol(Xs))]), unname(Xs))
})

test_that("renaming or reordering arms and cells leaves predictions and losses unchanged", {
  d <- small(scenario = "two_arms", n_arms = 3L)
  r1 <- run(d)
  d2 <- d; d2$cell_data$arm <- c(trunk = "trunk", A = "zeta", B = "alpha", C = "mid")[d$cell_data$arm]
  r2 <- run(d2)
  o <- rev(seq_len(ncol(d$counts))); d3 <- d; d3$counts <- d$counts[, o]; d3$cell_data <- d$cell_data[o, ]; rownames(d3$cell_data) <- NULL
  r3 <- run(d3)
  for (r in list(r2, r3)) {
    expect_identical(r$gene_fold$status, r1$gene_fold$status)
    expect_equal(r$gene_fold$single_loss, r1$gene_fold$single_loss, tolerance = 1e-6)
    expect_equal(r$gene_fold$branch_loss, r1$gene_fold$branch_loss, tolerance = 1e-6)
    expect_equal(r$summary$delta, r1$summary$delta, tolerance = 1e-6)
  }
  k <- match(r1$cells$cell_id, r3$cells$cell_id)
  expect_equal(r3$cell_nll$branch[, k], r1$cell_nll$branch, tolerance = 1e-6)
})

# ── 평가 방식 ────────────────────────────────────────────────────────────────
test_that("several donors give leave-one-donor-out with disjoint donors, not a within-donor split", {
  d <- small(scenario = "two_arms", n_donors = 3L)
  r <- run(d)
  expect_identical(r$summary$evaluation_scheme, "leave_one_donor_out")
  expect_identical(r$summary$evaluation_tier, "fixed_graph_lodo_stability")
  expect_identical(r$summary$n_donors_total, 3L)
  expect_true(all(r$folds$n_eval_donors == 1L))
  expect_true(all(r$folds$n_overlap_donors == 0L))
  expect_false(any(grepl("within_donor", c(r$summary$evaluation_scheme, r$summary$evaluation_tier, r$folds$evaluation_scheme))))
  expect_false(r$summary$generalization_verified)
  expect_true(r$summary$exploratory)
})

test_that("one donor with replicates is a replicate holdout, not a donor generalisation", {
  d <- rbc_simulate(scenario = "one_donor_replicates", cells_per_donor = 240L, n_genes = 6L)
  r <- run(d)
  expect_identical(r$summary$evaluation_scheme, "replicate_holdout")
  expect_identical(r$summary$evaluation_tier, "fixed_graph_replicate_stability")
  expect_identical(r$summary$n_donors_total, 1L)
  expect_false(r$summary$generalization_verified)
  expect_true(r$summary$exploratory)
})

test_that("one donor without replicates is a seeded within-donor cell split", {
  d <- rbc_simulate(scenario = "one_donor", cells_per_donor = 240L, n_genes = 6L)
  r <- run(d)
  expect_identical(r$summary$evaluation_scheme, "within_donor_cell_split")
  expect_identical(r$summary$evaluation_tier, "fixed_graph_cell_split_stability")
  expect_false(r$summary$generalization_verified)
  expect_true(r$summary$exploratory)
  expect_identical(r$params$k_used, 5L)
  expect_true(all(r$cells$evaluated))
  expect_identical(nrow(r$folds), 5L)
  for (a in unique(d$cell_data$arm)) expect_identical(length(unique(r$cells$fold_id[d$cell_data$arm == a])), 5L)
})

test_that("an arm confined to one donor leaves the folds unavailable instead of predicting an unseen arm", {
  d <- small(scenario = "arm_donor_overlap")
  r <- run(d)
  expect_true(all(r$folds$status == "unavailable"))
  expect_true(all(grepl("arm_missing_in_training", r$folds$reason)))
  expect_identical(r$summary$label, "unresolved")
  expect_identical(r$summary$n_genes_compared, 0L)
  expect_false(any(r$cells$evaluated))
})

test_that("user folds report donor overlap between training and evaluation", {
  d <- small(scenario = "two_arms", n_donors = 3L)
  f <- rep(c(1L, 2L, 3L), length.out = nrow(d$cell_data))
  r <- run(d, folds = f)
  expect_identical(r$summary$evaluation_scheme, "user_folds")
  expect_true(all(r$folds$n_overlap_donors == 3L))
})

# ── 기대 동작 ────────────────────────────────────────────────────────────────
test_that("the single path is a curve, not a line", {
  d <- small(scenario = "single_curve")
  expect_gt(run(d, spline_df = 1L)$summary$single_loss, run(d, spline_df = 3L)$summary$single_loss)
})

test_that("a smooth single curved path gives no meaningful branch advantage", {
  d <- small(scenario = "single_curve")
  r <- run(d)
  expect_lt(r$summary$delta, 0.01)
})

test_that("diverging arms give the branch model a held-out prediction advantage", {
  d <- small(scenario = "two_arms")
  r <- run(d)
  expect_gt(r$summary$delta, 0)
  expect_identical(r$summary$label, "branch_preferred_exploratory")
  g <- r$genes[r$genes$status == "compared", ]
  expect_gt(mean(g$delta[g$gene %in% d$truth$branch_genes]), mean(g$delta[!g$gene %in% d$truth$branch_genes]))
})

test_that("all-zero genes and failed fits are excluded from both models, never counted as evidence", {
  d <- small(scenario = "two_arms", zero_genes = 1L)
  z <- d$truth$zero_genes
  only_d1 <- rownames(d$counts)[1]; d$counts[only_d1, d$cell_data$donor_id != "D1"] <- 0L
  r <- run(d)
  expect_identical(r$genes$status[r$genes$gene == z], "excluded")
  expect_identical(r$genes$reason[r$genes$gene == z], "all_zero_gene")
  expect_identical(r$genes$status[r$genes$gene == only_d1], "excluded")
  expect_match(r$genes$reason[r$genes$gene == only_d1], "all_zero_in_training")
  expect_true(all(is.na(r$genes$delta[r$genes$status == "excluded"])))
  expect_identical(r$summary$n_genes_compared, sum(r$genes$status == "compared"))
  expect_identical(r$summary$n_genes_excluded, 2L)
  ok <- r$gene_fold$status == "ok"
  expect_true(all(r$gene_fold$n_eval_cells[ok] > 0))
  expect_true(all(is.na(r$gene_fold$single_loss[!ok]) & is.na(r$gene_fold$branch_loss[!ok])))
})

test_that("both models are scored on the same evaluation cells", {
  d <- small(scenario = "two_arms")
  r <- run(d)
  expect_identical(dim(r$cell_nll$single), dim(r$cell_nll$branch))
  expect_identical(is.na(r$cell_nll$single), is.na(r$cell_nll$branch))
  ev <- r$cells$evaluated
  expect_false(anyNA(r$cell_nll$single[r$genes$status == "compared", ev]))
})

test_that("the same seed reproduces results and the global random stream is untouched", {
  d <- rbc_simulate(scenario = "one_donor", cells_per_donor = 240L, n_genes = 4L)
  set.seed(99); a <- stats::runif(1)
  set.seed(99); r1 <- run(d, seed = 7L); b <- stats::runif(1)
  expect_identical(a, b)
  r2 <- run(d, seed = 7L)
  expect_identical(r1$cells$fold_id, r2$cells$fold_id)
  expect_identical(r1$summary, r2$summary)
})

test_that("cell order is preserved and out-of-range evaluation cells are flagged", {
  d <- small(scenario = "two_arms", n_donors = 3L)
  r <- run(d)
  expect_identical(r$cells$cell_id, d$cell_data$cell_id)
  cd <- d$cell_data; manual <- logical(nrow(cd))
  for (dn in unique(cd$donor_id)) { tr <- cd$position[cd$donor_id != dn]; i <- cd$donor_id == dn; manual[i] <- cd$position[i] < min(tr) | cd$position[i] > max(tr) }
  expect_identical(r$cells$out_of_training_range, manual)
  expect_identical(sum(r$folds$n_eval_out_of_range), sum(manual))
})

test_that("gene provenance defaults to unknown and is carried into the summary", {
  d <- small(scenario = "two_arms")
  expect_identical(run(d)$summary$gene_provenance, "unknown")
  expect_identical(run(d, gene_provenance = "same_as_route_genes")$summary$gene_provenance, "same_as_route_genes")
})

test_that("labels never claim probabilities or verified branching", {
  d <- small(scenario = "two_arms")
  r <- run(d)
  expect_true(all(c(r$summary$label, r$genes$label) %in% c("branch_preferred_exploratory", "single_preferred_exploratory", "unresolved")))
  expect_false(any(grepl("p_value|pvalue|prob|accuracy|verified_branch", names(r$summary))))
})
