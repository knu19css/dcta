# =============================================================================
# simulate_branch_data.R — 테스트 · 공개 예제용 작은 합성 count 자료 (실제 세포 자료가 아니다)
# =============================================================================
# scenario
#   single_curve         : 모든 세포가 하나의 휜 곡선(log 평균 = a + b·sin(pi·p))을 따른다. 분기 뒤 arm 이름만 붙는다
#   two_arms             : 분기 뒤 앞쪽 절반 유전자(branch_genes)가 arm 마다 반대 방향으로 갈라진다
#   arm_donor_overlap    : donor 2 명, 분기 뒤 arm 이 donor 와 완전히 겹친다 (D1 = A, D2 = B)
#   one_donor            : donor 1 명, 반복 없음
#   one_donor_replicates : donor 1 명, replicate_id 3 개 (R1 · R2 · R3)
# 반환: list(counts = 유전자 × 세포 정수 행렬, cell_data, truth = list(split_point, branch_genes, zero_genes, sparse_genes))
# 전역 난수 상태는 바꾸지 않는다.
# =============================================================================

rbc_simulate <- function(scenario = c("single_curve", "two_arms", "arm_donor_overlap", "one_donor", "one_donor_replicates"),
                         n_donors = 3L, cells_per_donor = 120L, n_genes = 12L, n_arms = 2L, split_point = 0.5,
                         theta = 8, zero_genes = 0L, sparse_genes = 0L, seed = 1L) {
  scenario <- match.arg(scenario)
  if (scenario %in% c("one_donor", "one_donor_replicates")) n_donors <- 1L
  if (scenario == "arm_donor_overlap") { n_donors <- 2L; n_arms <- 2L }
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE); old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  arms <- LETTERS[seq_len(n_arms)]; n <- n_donors * cells_per_donor
  donor <- rep(sprintf("D%d", seq_len(n_donors)), each = cells_per_donor); cell <- sprintf("cell%04d", seq_len(n))
  pos <- stats::runif(n); post <- pos > split_point
  arm <- rep("trunk", n); arm[post] <- arms[sample.int(n_arms, sum(post), replace = TRUE)]
  if (scenario == "arm_donor_overlap") arm[post] <- ifelse(donor[post] == "D1", "A", "B")
  rep_id <- if (scenario == "one_donor_replicates") paste0("R", (seq_len(n) - 1L) %% 3L + 1L) else donor
  lib <- round(stats::rlnorm(n, log(3000), 0.3)); sf <- lib / stats::median(lib); h <- pmax(pos - split_point, 0)
  genes <- sprintf("G%02d", seq_len(n_genes)); nb <- if (scenario == "single_curve") 0L else as.integer(ceiling(n_genes / 2))
  base <- stats::runif(n_genes, log(3), log(10)); amp <- stats::runif(n_genes, 0.8, 1.5); dev <- stats::runif(n_genes, 1.0, 1.6)
  donor_eff <- matrix(stats::rnorm(n_genes * n_donors, 0, 0.1), n_genes, n_donors, dimnames = list(genes, unique(donor)))
  eta <- outer(base, rep(1, n)) + outer(amp, sin(pi * pos)) + donor_eff[, donor, drop = FALSE]
  if (nb > 0L) {
    sgn <- stats::setNames(seq(1, -1, length.out = n_arms), arms); s <- ifelse(arm == "trunk", 0, sgn[arm])
    eta[seq_len(nb), ] <- eta[seq_len(nb), , drop = FALSE] + outer(dev[seq_len(nb)], s * h / (1 - split_point))
  }
  mu <- exp(eta) * rep(sf, each = n_genes)
  y <- matrix(stats::rnbinom(length(mu), size = theta, mu = mu), n_genes, n, dimnames = list(genes, cell))
  zg <- if (zero_genes > 0L) utils::tail(genes, zero_genes) else character(0)
  sg <- if (sparse_genes > 0L) utils::head(utils::tail(genes, zero_genes + sparse_genes), sparse_genes) else character(0)
  for (g in sg) y[g, ] <- stats::rnbinom(n, size = theta, mu = 0.05 * sf)
  y[zg, ] <- 0
  storage.mode(y) <- "integer"
  list(counts = y,
       cell_data = data.frame(cell_id = cell, donor_id = donor, position = pos, arm = arm, library_size = lib, replicate_id = rep_id, stringsAsFactors = FALSE),
       truth = list(split_point = split_point, branch_genes = genes[seq_len(nb)], zero_genes = zg, sparse_genes = sg))
}
