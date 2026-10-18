# 실행: Rscript tests/R/branch_candidates/run_tests.R
# 저장소 루트는 이 파일 위치에서 세 단계 위로 정한다 (절대경로를 코드에 적지 않는다).
args <- commandArgs(trailingOnly = FALSE)
self <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(normalizePath(self)), "..", "..", ".."))
Sys.setenv(RBC_REPO_ROOT = root)
suppressPackageStartupMessages(library(testthat))
res <- testthat::test_file(file.path(root, "tests", "R", "branch_candidates", "test_branch_curve_comparison.R"), reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
bad <- sum(df$failed) + sum(df$error)
cat(sprintf("\nbranch_candidates tests: %d blocks, %d expectations, %d expectations failed, %d errors, %d skipped, %d warnings\n",
            nrow(df), sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning)))
quit(save = "no", status = if (bad > 0) 1L else 0L)
