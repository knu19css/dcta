# =============================================================================
# 17_marker_donor_uncertainty.R — 기존 marker-only 근거의 donor 기반 불확실성 진단
# =============================================================================
#
# 무엇을 하는가
#   기존 51 개 coarse 군집(CTX unit_res4.0_seed0)의 marker 점수(A, B, M)를 원 RNA 에서 다시 만든다.
#   그다음 (군집, donor) 별로 나눠, 같은 marker 근거가 donor 마다 얼마나 재현되는지를 표로 남긴다.
#   기존 판정·패널·대조군·theta·군집은 바꾸지 않는다. 새 cutoff 를 만들지 않고 라벨을 바꾸지 않는다.
#
# 식 (legacy module_v1 / context-v2 와 같다)
#   x_cg     = log1p(count_cg / library_c * 1e4)    library = 전체 유전자 합, 중복 symbol 은 log 전에 합침
#   pos_c(s) = mean_{g in pos(s)} x_cg - mean_{j in ctrl(s+)} x_cj   (ctrl 은 저장된 0-based 원 열 그대로)
#   neg_c(s) = 같은 방식의 음성 module, 음성 패널이 비어 있으면 0, 쓸 수 없으면 NA
#   A_G(s), B_G(s) = 집단 G 안 세포 평균,  M_G(s) = A_G(s) - max(0, B_G(s) - theta(s))  (theta 는 저장값)
#   세포별 clipped M 의 평균이 아니다 (max 가 비선형이라 값이 달라진다).
#
# donor 층 (새로 더하는 진단, 해석은 README)
#   donor d 의 A_d, B_d, M_d  → equal-donor 평균(donor 같은 무게)과 세포 가중 pooled 를 이름을 달리해 둘 다 둔다.
#   donor 재표본(군집 안 donor 전체 프로필을 복원 추출) → 평균의 백분위 95% 구간, 승자 빈도(0 점 unassigned 포함, 동률 질량 보존)
#   정확한 donor 부호 검정: A 대 0, 모든 상태 쌍의 M_s - M_t. Holm(군집 안)과 BY(가족별 전체) 보정 p 를 함께 둔다.
#   donor / study 하나 빼기 영향 분석, 유전자별 기여(유전자 평균 - module 대조군 평균)
#   사후 저자 대조는 위 표를 고정한 뒤에만 한다(누설 차단).
#
# 이것은 "맞을 확률" 이 아니다. 같은 패널·대조군·군집 조건에서 donor 간 재현성을 기술할 뿐이다.
# =============================================================================

suppressPackageStartupMessages(library(Matrix))

MDU_TARGET_SUM <- 1e4
MDU_TIE_TOL <- 1e-12                     # 승자 동률 (legacy TIE_ATOL 과 같은 값)
MDU_UNASSIGNED <- "unassigned"
# 수치적 0 (Codex 14, 20): 부호 검정·유전자 기여 부호에서 부동소수 계산 잡음을 부호로 세지 않기 위한 공학적 guard.
#   값 = 2^14 x 기계 정밀도 x 2 x log1p(1e4) ≈ 6.7e-11. log 발현 척도(세포당 0 ~ log1p(1e4))와 기계 정밀도에 맞춘 보수적 설정이다.
#   모든 입력 크기에 대해 증명한 엄밀한 오차 상한이 아니다(clipped M 차이의 크기를 일반적으로 제한하지도 않는다).
#   생물학적 점수 cutoff 나 유의성 기준이 아니며, 관측된 작은 차이 분포에 맞춰 고른 값도 아니다.
MDU_NUM_TOL <- 2^14 * .Machine$double.eps * 2 * log1p(MDU_TARGET_SUM)
.mdu_sign3 <- function(x, tol = MDU_NUM_TOL) ifelse(is.na(x), NA_integer_, ifelse(x > tol, 1L, ifelse(x < -tol, -1L, 0L)))


# ── 0. 공통 검사 ────────────────────────────────────────────────────────────
#' 세포 ID 두 벡터가 순서까지 정확히 같아야 한다. 재정렬하지 않는다.
mdu_check_ids <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  if (anyDuplicated(a) || anyDuplicated(b)) stop("세포 ID 에 중복이 있다.")
  if (length(a) != length(b) || !identical(a, b)) stop("세포 ID 순서/길이가 다르다 (", length(a), " vs ", length(b), "). 재정렬하지 않는다.")
  TRUE
}

#' 적합 단계 세포 표에 저자 라벨·DAG 정보가 들어오지 않게 한다.
mdu_check_fit_columns <- function(cols) {
  bad <- cols[grepl("author|celltype|cell_type|label|prior|dag|ontology", cols, ignore.case = TRUE)]
  if (length(bad)) stop("저자/DAG 정보 열은 적합 단계에 넣지 않는다: ", paste(bad, collapse = ", "))
  TRUE
}

.mdu_split <- function(x) { x <- if (is.na(x)) "" else as.character(x); v <- strsplit(x, "|", fixed = TRUE)[[1]]; v[nzchar(v)] }


# ── 1. module 정의 (저장된 패널·대조군을 그대로 묶는다) ──────────────────────
#' 저장 state_table + controls 를 원 유전자 축(symbols)에 묶는다.
#' controls$col 은 0-based 원 열 번호다. 범위를 벗어나거나 패널 유전자를 가리키면 멈춘다.
mdu_module_definition <- function(state_table, controls, symbols) {
  st <- state_table; symbols <- as.character(symbols); nG <- length(symbols); states <- as.character(st$state)
  pos <- setNames(lapply(st$pos_genes, .mdu_split), states); neg <- setNames(lapply(st$neg_genes, .mdu_split), states)
  pos_ok <- setNames(st$pos_status == "ok", states)
  neg_mode <- setNames(ifelse(st$neg_status == "empty", "empty", ifelse(st$neg_status == "ok", "ok", "unavailable")), states)
  genes <- sort(unique(c(unlist(pos[pos_ok]), unlist(neg[neg_mode == "ok"]))))
  gcols <- lapply(genes, function(g) which(symbols == g))
  miss <- genes[lengths(gcols) == 0]
  if (length(miss)) stop("패널 유전자가 원 유전자 축에 없다: ", paste(miss, collapse = ", "))
  P <- Matrix::sparseMatrix(i = unlist(gcols), j = rep(seq_along(genes), lengths(gcols)), x = 1, dims = c(nG, length(genes)),
                            dimnames = list(NULL, genes))                                             # 중복 symbol 열을 합친다
  ctrl <- data.frame(module = as.character(controls$module), col = as.integer(controls$col), stringsAsFactors = FALSE)
  if (anyNA(ctrl$col) || any(ctrl$col < 0 | ctrl$col >= nG)) stop("대조군 열이 범위를 벗어났다 (0-based, 유전자 ", nG, " 개).")
  if (any(symbols[ctrl$col + 1L] %in% genes)) stop("대조군 열이 패널 유전자를 가리킨다.")
  ccols <- sort(unique(ctrl$col + 1L))
  qmat <- function(sign, use) {
    Q <- matrix(0, length(ccols), length(states), dimnames = list(NULL, states))
    for (s in states[use]) {
      c1 <- ctrl$col[ctrl$module == paste0(s, sign)] + 1L
      if (!length(c1)) stop("module ", s, sign, " 의 대조군이 없다.")
      if (anyDuplicated(c1)) stop("module ", s, sign, " 의 대조군 열이 중복이다.")
      Q[match(c1, ccols), s] <- 1 / length(c1)
    }
    Q
  }
  W <- function(sets, use) {
    M <- matrix(0, length(genes), length(states), dimnames = list(genes, states))
    for (s in states[use]) M[sets[[s]], s] <- 1 / length(sets[[s]])
    M
  }
  list(states = states, genes = genes, symbols = symbols, P = P, ctrl_cols = ccols,
       Qpos = qmat("+", pos_ok), Qneg = qmat("-", neg_mode == "ok"), Wpos = W(pos, pos_ok), Wneg = W(neg, neg_mode == "ok"),
       pos = pos, neg = neg, pos_ok = pos_ok, neg_mode = neg_mode, measured = pos_ok)
}


#' 패널 가용성 기록: 저장 state_table 의 listed/available 수와 상태를 그대로 남기고,
#' 원 패널(panel_adapted: name, markers, neg_markers)과 listed 수가 정확히 맞을 때만 빠진 유전자 이름을 적는다.
#' 맞지 않으면 추정하지 않고 멈춘다. alias 교정이나 패널 변경은 하지 않는다.
mdu_panel_availability <- function(state_table, panel) {
  ix <- match(state_table$state, panel$name)
  if (anyNA(ix)) stop("원 패널에 없는 상태가 있다: ", paste(state_table$state[is.na(ix)], collapse = ", "))
  do.call(rbind, lapply(seq_len(nrow(state_table)), function(i) {
    r <- state_table[i, ]; p <- panel[ix[i], ]
    do.call(rbind, lapply(c("pos", "neg"), function(side) {
      listed <- .mdu_split(if (side == "pos") p$markers else p$neg_markers)
      avail <- .mdu_split(if (side == "pos") r$pos_genes else r$neg_genes)
      nl <- as.integer(if (side == "pos") r$n_pos_listed else r$n_neg_listed); na <- as.integer(if (side == "pos") r$n_pos_available else r$n_neg_available)
      if (length(listed) != nl || length(avail) != na || !all(avail %in% listed))
        stop("원 패널과 state_table 이 맞지 않는다: ", r$state, " ", side)
      data.frame(state = r$state, side = side, n_listed = nl, n_available = na, status = if (side == "pos") r$pos_status else r$neg_status,
                 available_genes = paste(avail, collapse = "|"), missing_genes = paste(setdiff(listed, avail), collapse = "|"), stringsAsFactors = FALSE)
    }))
  }))
}


# ── 2. 세포별 module (chunk 하나; 전체 dense RNA 를 만들지 않는다) ─────────────
#' X: 세포 x 원 유전자 count (sparse). 반환: library, panel_log, det, pos/neg 대조군 평균, pos/neg module
mdu_chunk_modules <- function(X, def) {
  X <- methods::as(X, "CsparseMatrix")
  if (ncol(X) != length(def$symbols)) stop("count 행렬 열 수가 유전자 축과 다르다.")
  if (length(X@x) && (any(!is.finite(X@x)) || any(X@x < 0))) stop("count 에 비유한 값이나 음수가 있다.")
  lib <- Matrix::rowSums(X)
  if (any(lib <= 0)) stop("library 합이 0 인 세포가 있다.")
  cnt <- as.matrix(X %*% def$P)                                  # 중복 symbol 합 (log 전에)
  panel_log <- log1p(cnt / lib * MDU_TARGET_SUM); dimnames(panel_log) <- list(NULL, def$genes)
  det <- (cnt > 0) * 1; dimnames(det) <- dimnames(panel_log)
  Xc <- X[, def$ctrl_cols, drop = FALSE]
  Xc <- Matrix::Diagonal(x = MDU_TARGET_SUM / lib) %*% Xc; Xc <- methods::as(Xc, "CsparseMatrix"); Xc@x <- log1p(Xc@x)
  pc <- as.matrix(Xc %*% def$Qpos); nc <- as.matrix(Xc %*% def$Qneg)
  k <- def$states
  pos <- panel_log %*% def$Wpos - pc; neg <- panel_log %*% def$Wneg - nc
  pos[, !def$pos_ok] <- NA; pc[, !def$pos_ok] <- NA
  neg[, def$neg_mode == "empty"] <- 0; nc[, def$neg_mode == "empty"] <- NA
  neg[, def$neg_mode == "unavailable"] <- NA; nc[, def$neg_mode == "unavailable"] <- NA
  dimnames(pos) <- dimnames(neg) <- dimnames(pc) <- dimnames(nc) <- list(NULL, k)
  list(library = unname(lib), panel_log = panel_log, det = det, pos_ctrl = pc, neg_ctrl = nc, pos = pos, neg = neg)
}

#' 원 h5ad 의 raw/X(CSR)를 chunk 로 읽어 집단별 합을 모은다. group: 세포마다 1..n_groups 정수.
mdu_stream_h5 <- function(h5_path, def, group, n_groups, chunk = 5000L, log = function(...) NULL) {
  f <- hdf5r::H5File$new(h5_path, "r"); on.exit(f$close_all())
  rx <- f[["raw/X"]]
  if (!identical(as.character(hdf5r::h5attr(rx, "encoding-type")), "csr_matrix")) stop("raw/X 가 csr_matrix 가 아니다.")
  shp <- as.integer(hdf5r::h5attr(rx, "shape")); n <- shp[1]
  if (shp[2] != length(def$symbols) || length(group) != n) stop("raw/X 모양이 유전자 축/세포 수와 다르다.")
  indptr <- as.numeric(rx[["indptr"]]$read())
  if (length(indptr) != n + 1 || indptr[1] != 0 || any(diff(indptr) < 0)) stop("raw/X indptr 가 올바르지 않다.")
  K <- length(def$states); G <- length(def$genes)
  acc <- list(n = numeric(n_groups),
              pos = matrix(0, n_groups, K, dimnames = list(NULL, def$states)), neg = matrix(0, n_groups, K, dimnames = list(NULL, def$states)),
              pos_ctrl = matrix(0, n_groups, K, dimnames = list(NULL, def$states)), neg_ctrl = matrix(0, n_groups, K, dimnames = list(NULL, def$states)),
              panel_log = matrix(0, n_groups, G, dimnames = list(NULL, def$genes)), det = matrix(0, n_groups, G, dimnames = list(NULL, def$genes)))
  ds <- rx[["data"]]; di <- rx[["indices"]]; nonint <- 0
  for (s0 in seq(0, n - 1, by = chunk)) {
    s1 <- min(s0 + chunk, n); a <- indptr[s0 + 1]; b <- indptr[s1 + 1]; nr <- s1 - s0
    x <- if (b > a) as.numeric(ds[(a + 1):b]) else numeric(0); j <- if (b > a) as.integer(di[(a + 1):b]) + 1L else integer(0)
    if (any(x != round(x))) nonint <- nonint + sum(x != round(x))
    X <- Matrix::sparseMatrix(i = rep.int(seq_len(nr), diff(indptr[(s0 + 1):(s1 + 1)])), j = j, x = x, dims = c(nr, length(def$symbols)))
    cm <- mdu_chunk_modules(X, def); g <- group[(s0 + 1):s1]
    ug <- sort(unique(g))
    acc$n[ug] <- acc$n[ug] + tabulate(g, n_groups)[ug]
    for (k in c("pos", "neg", "pos_ctrl", "neg_ctrl", "panel_log", "det")) {
      rs <- rowsum(cm[[k]], g, reorder = TRUE); acc[[k]][as.integer(rownames(rs)), ] <- acc[[k]][as.integer(rownames(rs)), ] + rs
    }
    if ((s0 / chunk) %% 10 == 0) log(sprintf("stream rows %d-%d / %d", s0 + 1, s1, n))
  }
  if (nonint > 0) stop("raw/X 에 정수가 아닌 값이 ", nonint, " 개 있다: 원 count 가 아니다.")
  acc
}


# ── 3. 집계와 기준 대조 ─────────────────────────────────────────────────────
#' 집단 평균 A, B (세포별 module 의 평균)
mdu_group_means <- function(A, B, group) {
  n <- table(group); g <- names(n)
  list(A = rowsum(A, group)[g, , drop = FALSE] / as.numeric(n), B = rowsum(B, group)[g, , drop = FALSE] / as.numeric(n),
       n = setNames(as.numeric(n), g))
}

#' M = A - max(0, B - theta) (열마다 theta), NA 는 그대로
mdu_marker_M <- function(A, B, theta) {
  th <- theta[colnames(A)]
  if (anyNA(th)) stop("theta 가 없는 상태가 있다: ", paste(colnames(A)[is.na(th)], collapse = ", "))
  A - pmax(0, sweep(B, 2, th, "-"))
}

#' 다시 만든 pooled A/B/M 을 저장 기준(unit_state_evidence)과 비교한다. 어긋나면 멈춘다.
mdu_validate_reference <- function(A, B, ref, tol = 1e-9) {
  ref <- ref[ref$state %in% colnames(A), ]
  theta <- tapply(ref$theta, ref$state, function(x) { u <- unique(x[!is.na(x)]); if (length(u) > 1) stop("theta 가 군집마다 다르다."); if (length(u)) u else NA_real_ })
  theta <- setNames(as.numeric(theta[colnames(A)]), colnames(A))
  mat <- function(v) { m <- matrix(NA_real_, nrow(A), ncol(A), dimnames = dimnames(A)); m[cbind(match(as.character(ref$unit), rownames(A)), match(ref$state, colnames(A)))] <- v; m }
  RA <- mat(ref$A); RB <- mat(ref$B); RM <- mat(ref$M); M <- mdu_marker_M(A, B, theta)
  d <- c(A = NA, B = NA, M = NA)
  for (k in names(d)) {
    x <- list(A = A, B = B, M = M)[[k]]; r <- list(A = RA, B = RB, M = RM)[[k]]
    if (!identical(is.na(x), is.na(r))) stop("저장 기준값과 결측 위치가 다르다: ", k)
    d[[k]] <- if (all(is.na(x))) 0 else max(abs(x - r), na.rm = TRUE)
    if (d[[k]] > tol) stop(sprintf("저장 기준값과 다르다: %s 최대 오차 %.3g > %.1g", k, d[[k]], tol))
  }
  list(max_abs_diff = d, theta = theta, M = M, tol = tol)
}


# ── 4. donor 요약 ────────────────────────────────────────────────────────────
#' 입력: (cluster, donor, study, n_cells, state, A, B, theta) 긴 표. donor 가 없는 군집은 행이 없다(반대표가 아니다).
mdu_donor_summary <- function(dl) {
  dl$M <- dl$A - pmax(0, dl$B - dl$theta)
  key <- paste(dl$cluster, dl$state, sep = "\r")
  do.call(rbind, lapply(split(dl, factor(key, levels = unique(key))), function(x) {
    x <- x[order(x$donor), ]; ok <- !is.na(x$A); w <- x$n_cells
    Ap <- sum(w * x$A) / sum(w); Bp <- sum(w * x$B) / sum(w)
    data.frame(cluster = x$cluster[1], state = x$state[1], n_donors = sum(ok), n_studies = length(unique(x$study[ok])),
               n_cells = sum(w), donor_cells = paste0(x$donor, ":", x$n_cells, collapse = ";"),
               A_donor_mean = mean(x$A[ok]), A_donor_sd = if (sum(ok) > 1) stats::sd(x$A[ok]) else NA_real_,
               B_donor_mean = mean(x$B[ok]), M_donor_mean = mean(x$M[ok]), M_donor_sd = if (sum(ok) > 1) stats::sd(x$M[ok]) else NA_real_,
               A_pooled_cell_weighted = Ap, B_pooled_cell_weighted = Bp,
               M_pooled_cell_weighted = Ap - max(0, Bp - x$theta[1]), theta = x$theta[1], stringsAsFactors = FALSE)
  }))
}


# ── 5. 승자 (0 점 unassigned 선택지, 동률 질량 보존) ───────────────────────────
#' 측정 점수가 하나도 없으면(모두 NA) unassigned 가 이긴 것으로 보지 않고 모두 NA 를 돌려준다 (Codex 17).
mdu_winner_mass <- function(score, tol = MDU_TIE_TOL) {
  s <- c(score, setNames(0, MDU_UNASSIGNED)); out <- setNames(numeric(length(s)), names(s))
  if (all(is.na(score))) return(out * NA)
  ok <- !is.na(s); top <- max(s[ok]); w <- ok & s >= top - tol
  out[w] <- 1 / sum(w); out
}
.mdu_winner_label <- function(mass) paste(names(mass)[mass > 0], collapse = "|")


# ── 6. donor 재표본 ──────────────────────────────────────────────────────────
#' A, M: donor x 상태 행렬 (한 군집). 같은 donor 추출을 모든 상태에 쓴다. R 회(>= 2), seed 고정.
#'  - 상태별 평균: 그 상태 값이 있는 donor 만 쓴다. nonmissing donor < 2 이면 CI 를 내지 않는다(insufficient_donors).
#'  - 쌍 대비: donor 마다 M_s - M_t 를 먼저 만든다(둘 다 있는 공통 donor). 그다음 같은 추출로 평균한다 (Codex 사전검토 3).
#'  - 승자: 0 점 unassigned 포함, 동률은 질량을 나눈다. MC SE = 반복별 질량의 SD / sqrt(R) (Codex 사전검토 4).
#'    동률 반복 비중을 따로 남긴다. MC SE 는 모의실험 오차이며 생물학적 오차가 아니다.
mdu_bootstrap_cluster <- function(A, M, R = 2000L, seed = 20261011L) {
  R <- as.integer(R)
  if (length(R) != 1 || is.na(R) || R < 2L) stop("재표본 횟수 R 은 2 이상의 정수여야 한다.")
  nd <- nrow(M); k <- colnames(M); pairs <- if (length(k) > 1) t(utils::combn(k, 2)) else matrix(character(0), 0, 2)
  q <- function(x) unname(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, type = 7))
  stat_of <- function(n_nonmiss, x = NULL) {
    if (n_nonmiss == 0) return("all_missing")
    if (n_nonmiss < 2) return("insufficient_donors")
    if (length(unique(x[!is.na(x)])) == 1) "degenerate" else "ok"
  }
  nonmiss <- colSums(!is.na(M)); nonmissA <- colSums(!is.na(A))
  D <- if (nrow(pairs)) M[, pairs[, 1], drop = FALSE] - M[, pairs[, 2], drop = FALSE] else matrix(numeric(0), nd, 0)
  npair <- colSums(!is.na(D))
  pair_base <- data.frame(state_s = pairs[, 1], state_t = pairs[, 2], n_donors_paired = as.integer(npair), stringsAsFactors = FALSE)
  if (nd < 2) {
    type <- data.frame(state = k, n_donors_nonmissing = as.integer(nonmiss), A_ci_low = NA_real_, A_ci_high = NA_real_,
                       M_ci_low = NA_real_, M_ci_high = NA_real_, M_boot_sd = NA_real_,
                       M_boot_status = vapply(nonmiss, stat_of, character(1)), row.names = NULL, stringsAsFactors = FALSE)
    return(list(status = "insufficient_donors", n_replicates = R, n_donors = nd, seed = seed, type = type,
                pairs = cbind(pair_base, diff_ci_low = rep(NA_real_, nrow(pairs)), diff_ci_high = rep(NA_real_, nrow(pairs)),
                              diff_boot_status = vapply(npair, stat_of, character(1))),
                winner = data.frame(candidate = c(k, MDU_UNASSIGNED), freq = NA_real_, mc_se = NA_real_, tie_draw_frac = NA_real_),
                tie_draw_frac = NA_real_, n_draws_available = NA_integer_))
  }
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  idx <- matrix(sample.int(nd, nd * R, replace = TRUE), R, nd)
  W <- t(apply(idx, 1, tabulate, nbins = nd))
  bmean <- function(X) { X0 <- X; X0[is.na(X0)] <- 0; den <- W %*% (!is.na(X)); m <- (W %*% X0) / den; m[den == 0] <- NA; colnames(m) <- colnames(X); m }
  Ab <- bmean(A); Mb <- bmean(M)
  ci <- function(B, n_nm) {
    st <- vapply(seq_len(ncol(B)), function(j) stat_of(n_nm[j], B[, j]), character(1))
    lo <- hi <- rep(NA_real_, ncol(B)); ok <- st %in% c("ok", "degenerate")
    for (j in which(ok)) { v <- q(B[, j]); lo[j] <- v[1]; hi[j] <- v[2] }
    list(lo = lo, hi = hi, st = st, sd = ifelse(ok, apply(B, 2, stats::sd, na.rm = TRUE), NA_real_))
  }
  cA <- ci(Ab, nonmissA); cM <- ci(Mb, nonmiss)
  type <- data.frame(state = k, n_donors_nonmissing = as.integer(nonmiss), A_ci_low = cA$lo, A_ci_high = cA$hi,
                     M_ci_low = cM$lo, M_ci_high = cM$hi, M_boot_sd = cM$sd, M_boot_status = cM$st, row.names = NULL, stringsAsFactors = FALSE)
  pr <- if (nrow(pairs)) { cD <- ci(bmean(D), npair); cbind(pair_base, diff_ci_low = cD$lo, diff_ci_high = cD$hi, diff_boot_status = cD$st) } else
    cbind(pair_base, diff_ci_low = numeric(0), diff_ci_high = numeric(0), diff_boot_status = character(0))
  mass <- apply(Mb, 1, mdu_winner_mass)                       # (상태 + unassigned) x R, 측정 점수가 없는 반복은 NA 열
  avail <- !is.na(colSums(mass)); na_ <- sum(avail)          # 빈도의 분모 = 쓸 수 있는 반복 수 (Codex 17)
  if (na_ == 0) {
    winner <- data.frame(candidate = rownames(mass), freq = NA_real_, mc_se = NA_real_, tie_draw_frac = NA_real_, row.names = NULL, stringsAsFactors = FALSE)
    return(list(status = "all_missing", n_replicates = R, n_donors = nd, seed = seed, type = type, pairs = pr, winner = winner,
                tie_draw_frac = NA_real_, n_draws_available = 0L))
  }
  ms <- mass[, avail, drop = FALSE]; tie_draw <- colSums(ms > 0) > 1
  winner <- data.frame(candidate = rownames(ms), freq = rowMeans(ms), mc_se = if (na_ > 1) apply(ms, 1, stats::sd) / sqrt(na_) else NA_real_,
                       tie_draw_frac = rowMeans(ms > 0 & ms < 1), row.names = NULL, stringsAsFactors = FALSE)
  list(status = "ok", n_replicates = R, n_donors = nd, seed = seed, type = type, pairs = pr, winner = winner, tie_draw_frac = mean(tie_draw),
       n_draws_available = as.integer(na_))
}


# ── 7. 정확한 donor 부호 검정 ────────────────────────────────────────────────
#' |x| <= MDU_NUM_TOL 는 부호가 아니라 0 으로 센다. 정확한 0 과 수치적 0 을 따로 남긴다 (Codex 14).
mdu_sign_test <- function(x, tol = MDU_NUM_TOL) {
  s <- .mdu_sign3(x, tol)
  np <- sum(s == 1L, na.rm = TRUE); nn <- sum(s == -1L, na.rm = TRUE); nz <- sum(s == 0L, na.rm = TRUE); nm <- sum(is.na(x)); nn_ <- np + nn + nz
  ze <- sum(x == 0, na.rm = TRUE)
  st <- if (nn_ == 0) "no_donors" else if (nn_ == 1) "single_donor" else if (np + nn == 0) "all_zero" else "ok"
  p <- switch(st, no_donors = NA_real_, single_donor = NA_real_, all_zero = 1, ok = stats::binom.test(np, np + nn, 0.5)$p.value)
  data.frame(n_pos = as.integer(np), n_neg = as.integer(nn), n_zero = as.integer(nz), n_zero_exact = as.integer(ze),
             n_zero_numerical = as.integer(nz - ze), n_missing = as.integer(nm), n_nonmissing = as.integer(nn_), p_exact = p, status = st,
             stringsAsFactors = FALSE)
}

#' 한 군집: A 대 0 (상태마다), M_s - M_t (모든 비순서 쌍)
mdu_sign_tests_cluster <- function(A, M, cluster) {
  k <- colnames(M)
  a <- do.call(rbind, lapply(k, function(s) cbind(data.frame(cluster = cluster, family = "A_support", state = s), mdu_sign_test(A[, s]))))
  pw <- if (length(k) > 1) {
    pr <- t(utils::combn(k, 2))
    do.call(rbind, lapply(seq_len(nrow(pr)), function(i) cbind(data.frame(cluster = cluster, family = "pairwise", state_s = pr[i, 1], state_t = pr[i, 2]),
                                                                mdu_sign_test(M[, pr[i, 1]] - M[, pr[i, 2]]))))
  } else NULL
  list(a_support = a, pairwise = pw)
}

#' Holm (군집 x 가족 안) 와 BY (가족 전체). NA p 는 보정에서 빠진다. 아무 행도 지우지 않는다.
mdu_adjust <- function(df) {
  df$p_holm_cluster <- NA_real_; df$p_by_global <- NA_real_
  for (g in split(seq_len(nrow(df)), paste(df$cluster, df$family, sep = "\r"))) {
    ok <- g[!is.na(df$p_exact[g])]; if (length(ok)) df$p_holm_cluster[ok] <- stats::p.adjust(df$p_exact[ok], "holm")
  }
  for (g in split(seq_len(nrow(df)), df$family)) {
    ok <- g[!is.na(df$p_exact[g])]; if (length(ok)) df$p_by_global[ok] <- stats::p.adjust(df$p_exact[ok], "BY")
  }
  df
}


# ── 8. donor / study 하나 빼기 ───────────────────────────────────────────────
#' 전체 승자가 동률이면 대비를 임의로 고르지 않는다(NA, contrast_status = full_tie). 남는 donor 가 없으면 insufficient_donors.
mdu_leave_one_out <- function(M, groups, what = c("donor", "study")) {
  what <- match.arg(what); ug <- sort(unique(groups))
  if (what == "study" && length(ug) < 2) return(data.frame(what = character(0), left_out = character(0), n_left = integer(0), winner = character(0),
                                                            n_tied = integer(0), full_winner = character(0), full_winner_score_after = numeric(0),
                                                            best_other_score_after = numeric(0), contrast_after = numeric(0), contrast_status = character(0)))
  mean_na <- function(X) { m <- colMeans(X, na.rm = TRUE); m[is.nan(m)] <- NA; m }
  full <- mdu_winner_mass(mean_na(M)); full_na <- all(is.na(full))
  full_tie <- !full_na && sum(full > 0) > 1; fw <- if (full_tie || full_na) NA_character_ else names(full)[full > 0]
  empty <- function(base, status) cbind(base, winner = NA_character_, n_tied = NA_integer_, full_winner_score_after = NA_real_,
                                        best_other_score_after = NA_real_, contrast_after = NA_real_, contrast_status = status)
  do.call(rbind, lapply(ug, function(g) {
    keep <- groups != g; nl <- as.integer(sum(keep))
    base <- data.frame(what = what, left_out = g, n_left = nl, full_winner = if (full_na) NA_character_ else .mdu_winner_label(full), stringsAsFactors = FALSE)
    if (nl == 0) return(empty(base, "insufficient_donors"))
    s <- mean_na(M[keep, , drop = FALSE])
    if (full_na || all(is.na(s))) return(empty(base, "all_missing"))                 # 측정 점수가 없으면 unassigned 확정처럼 보이지 않게 (Codex 17)
    ms <- mdu_winner_mass(s); sc <- c(s, setNames(0, MDU_UNASSIGNED))
    if (full_tie) return(cbind(base, winner = .mdu_winner_label(ms), n_tied = sum(ms > 0), full_winner_score_after = NA_real_,
                               best_other_score_after = NA_real_, contrast_after = NA_real_, contrast_status = "full_tie"))
    oth <- sc[names(sc) != fw]; bo <- if (all(is.na(oth))) NA_real_ else max(oth, na.rm = TRUE)
    cbind(base, winner = .mdu_winner_label(ms), n_tied = sum(ms > 0), full_winner_score_after = unname(sc[fw]),
          best_other_score_after = bo, contrast_after = unname(sc[fw]) - bo, contrast_status = "ok")
  }))
}


# ── 9. 유전자 근거 ──────────────────────────────────────────────────────────
#' 입력: (cluster, donor, state, side, gene, gene_mean, ctrl_mean, detect_frac). 기여 = gene_mean - ctrl_mean
mdu_gene_summary <- function(ge) {
  ge$contribution <- ge$gene_mean - ge$ctrl_mean
  key <- paste(ge$cluster, ge$state, ge$side, ge$gene, sep = "\r")
  do.call(rbind, lapply(split(ge, factor(key, levels = unique(key))), function(x) {
    v <- x$contribution; s <- .mdu_sign3(v)
    data.frame(cluster = x$cluster[1], state = x$state[1], side = x$side[1], gene = x$gene[1], n_donors = sum(!is.na(v)),
               n_donors_pos = sum(s == 1L, na.rm = TRUE), n_donors_neg = sum(s == -1L, na.rm = TRUE), n_donors_zero = sum(s == 0L, na.rm = TRUE),
               n_donors_zero_exact = sum(v == 0, na.rm = TRUE), n_donors_zero_numerical = sum(s == 0L & v != 0, na.rm = TRUE),
               frac_donors_pos = mean(s == 1L, na.rm = TRUE), contribution_median = stats::median(v, na.rm = TRUE),
               contribution_mean = mean(v, na.rm = TRUE), detect_frac_median = stats::median(x$detect_frac, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
}


# ── 9b. 안정성 대 사후 저자 일치 (군집 표) ──────────────────────────────────
#' calls(unit, call, status), per_cluster(저자 대조 군집별), meta(재표본 군집 상태), win(재표본 승자 빈도)
#' boot_freq_existing_call_equal_donor: equal-donor 평균 재표본에서 기존 판정(no-call 이면 unassigned)이 이긴 빈도.
#' 기존 판정은 세포 가중 pooled 점수에서 나왔으므로, 이 빈도는 donor 안정성, donor 가중 방식,
#' 음성 penalty max(0, B - theta) 의 집계 순서(pooled B 에 한 번 vs donor 마다)를 함께 반영한다.
mdu_stability_table <- function(calls, per_cluster, meta, win) {
  s <- data.frame(cluster = as.character(calls$unit), existing_call = as.character(calls$call), call_status = as.character(calls$status), stringsAsFactors = FALSE)
  s$existing_call[is.na(s$existing_call)] <- ""
  pc <- per_cluster; names(pc)[names(pc) != "cluster"] <- paste0("author_", names(pc)[names(pc) != "cluster"])
  s <- merge(s, pc, by = "cluster", all.x = TRUE, sort = FALSE)
  mt <- data.frame(cluster = as.character(meta$cluster), bootstrap_status = meta$status, n_donors = meta$n_donors, n_replicates = meta$n_replicates,
                   seed = meta$seed, tie_draw_frac = meta$tie_draw_frac, stringsAsFactors = FALSE)
  s <- merge(s, mt, by = "cluster", all.x = TRUE, sort = FALSE)
  cand <- ifelse(s$existing_call == "", MDU_UNASSIGNED, s$existing_call); key <- paste(win$cluster, win$candidate)
  s$boot_freq_existing_call_equal_donor <- win$freq[match(paste(s$cluster, cand), key)]
  s$boot_mc_se_existing_call_equal_donor <- win$mc_se[match(paste(s$cluster, cand), key)]
  s[order(suppressWarnings(as.integer(s$cluster)), s$cluster), ]
}


# ── 10. 사후 저자 대조 (적합 표를 고정한 뒤에만 부른다) ────────────────────────
#' 저자 fine 라벨 → 패널 상태. 한 라벨이 두 상태에 대응하면 멈춘다.
mdu_expected_state <- function(map) {
  lab <- lapply(map$exact_fine, .mdu_split); st <- rep(map$state, lengths(lab)); lab <- unlist(lab)
  d <- lab[duplicated(lab)]
  if (length(d)) stop("저자 라벨이 둘 이상의 상태에 대응한다: ", paste(unique(d), collapse = ", "))
  setNames(st, lab)
}

mdu_author_audit <- function(cells, calls, map) {
  ex <- mdu_expected_state(map)
  ix <- match(as.character(cells$cluster), as.character(calls$unit))
  if (anyNA(ix)) stop("판정 표에 없는 군집이 있다.")
  call <- as.character(calls$call[ix]); call[is.na(call)] <- ""; status <- as.character(calls$status[ix])
  expected <- unname(ex[as.character(cells$author)])
  outcome <- ifelse(is.na(expected), "unmapped_reference",
                    ifelse(call == "", ifelse(startsWith(status, "unavailable"), "unavailable", "held"), ifelse(call == expected, "correct", "wrong")))
  cells <- data.frame(cells, expected_state = expected, call = call, call_status = status, outcome = outcome, stringsAsFactors = FALSE)
  cnt <- function(x) {
    o <- x$outcome; ev <- sum(o != "unmapped_reference"); C <- sum(o == "correct"); W <- sum(o == "wrong")
    data.frame(n_total = nrow(x), n_evaluable = ev, n_correct = C, n_wrong = W, n_held = sum(o == "held"), n_unavailable = sum(o == "unavailable"),
               n_unmapped_reference = sum(o == "unmapped_reference"), recovery = if (ev) C / ev else NA_real_,
               coverage = if (ev) (C + W) / ev else NA_real_, precision = if (C + W) C / (C + W) else NA_real_)
  }
  # 자료와 묶음 변수를 함께 받는다 (Codex 사전검토 7: 타입별 분모는 평가 가능 세포만)
  by <- function(d, v) { s <- split(d, factor(v, levels = sort(unique(v)))); do.call(rbind, lapply(names(s), function(k) cbind(key = k, cnt(s[[k]])))) }
  pd <- by(cells, cells$donor); names(pd)[1] <- "donor"; pc <- by(cells, as.character(cells$cluster)); names(pc)[1] <- "cluster"
  ev <- cells[!is.na(cells$expected_state), ]; pt <- by(ev, ev$expected_state); names(pt)[1] <- "expected_state"
  # 정답 타입별 층에서 C/(C+W) 는 PPV 가 아니라 그 타입 안의 조건부 회수율이다 (Codex 15)
  names(pt)[names(pt) == "precision"] <- "conditional_recovery"
  # 예측 타입별 PPV: 평가 가능 세포 중 그 상태로 판정된 세포가 분모. 판정 0 이면 NA (Codex 15)
  types <- sort(unique(c(ev$expected_state, ev$call[ev$call != ""])))
  pp <- do.call(rbind, lapply(types, function(t) {
    k <- ev$call == t; n <- sum(k); C <- sum(k & ev$outcome == "correct")
    data.frame(predicted_state = t, n_called_evaluable = as.integer(n), n_correct = as.integer(C), n_wrong = as.integer(n - C),
               ppv = if (n) C / n else NA_real_, stringsAsFactors = FALSE)
  }))
  al <- ifelse(is.na(cells$author), "<NA>", as.character(cells$author))
  cda <- stats::aggregate(list(n_cells = rep(1L, nrow(cells))), list(cluster = as.character(cells$cluster), donor = cells$donor, author_label = al), sum)
  list(cells = cells, overall = cnt(cells), per_donor = pd, per_cluster = pc, per_type = pt, per_predicted = pp, cluster_donor_author = cda)
}
