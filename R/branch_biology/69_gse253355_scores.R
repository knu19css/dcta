# =============================================================================
# 69_gse253355_scores.R — GSE253355: 고정 atlas 6x10 패널의 유전자 대응 · 세포 점수 · 군집 점수 (raw count, 모든 유전자 library)
# =============================================================================
#
# 계약: codex_gse253355_instruction.md 2–3 항.
#   - 패널은 attempt_02 source 고정본(atlas_stage_markers.csv) 그대로다. 유전자 이름은 정확한 symbol 만 대응한다(이름 변경 · 대체 · 보충 없음).
#     없는 유전자는 missing 으로 남긴다. HGNC 이전 · 별칭 symbol 이 자료에 있으면 진단 열로만 적고 점수에는 쓰지 않는다. GENE.1 같은 접미사 변형도 진단으로만 적는다.
#   - 세포 값 = log1p(1e4 · raw count / 모든 유전자 library). 단계 점수(주) = 목록 유전자 수(10)를 분모로 한 평균(측정 안 된 유전자는 0 으로 더한다, 보수적).
#     진단 = 측정된 유전자 수를 분모로 한 평균. 측정된 유전자가 하나도 없으면 두 값 모두 결측(판정불가).
#   - 군집 점수 = 0 인 세포를 포함한 모든 세포 평균. 보조 = donor 동일 가중(donor x 군집 >= min_cells 세포인 donor 만, 없으면 결측 — pooled 로 대체하지 않는다).
# =============================================================================

#' 패널 유전자 → 자료 유전자 대응. hgnc = data.frame(symbol, prev_symbol, alias_symbol) (없으면 NULL)
gsx_panel_map <- function(features, panels, hgnc = NULL) {
  features <- as.character(features)
  d <- data.frame(stage = as.character(panels$stage), rank = as.integer(panels$rank), symbol = as.character(panels$symbol), stringsAsFactors = FALSE)
  if (!is.null(panels$ensg)) d$ensg <- as.character(panels$ensg)
  ne <- vapply(d$symbol, function(g) sum(features == g), 0L, USE.NAMES = FALSE)
  if (any(ne > 1L)) stop("자료 유전자 이름에 같은 symbol 이 두 번 이상 있다: ", paste(d$symbol[ne > 1L], collapse = ", "))
  d$status <- ifelse(ne == 1L, "exact", "missing"); d$feature_index <- ifelse(ne == 1L, match(d$symbol, features), NA_integer_)
  d$suffix_variants_in_data <- vapply(d$symbol, function(g) { pre <- paste0(g, "."); hit <- startsWith(features, pre) & grepl("^[0-9]+$", substring(features, nchar(pre) + 1L))
    paste(features[hit], collapse = "|") }, "", USE.NAMES = FALSE)                                # GENE.1 같은 접미사 변형 (정규식 이스케이프 없이)
  d$alias_candidates_in_data <- if (is.null(hgnc)) rep("", nrow(d)) else vapply(d$symbol, function(g) { z <- hgnc[hgnc$symbol == g, , drop = FALSE]
    if (!nrow(z) || g %in% features) return("")                                                     # 정확히 측정된 유전자는 별칭 진단을 하지 않는다 (다른 유전자 이름과 혼동 방지)
    cand <- trimws(unlist(strsplit(c(as.character(z$prev_symbol), as.character(z$alias_symbol)), "|", fixed = TRUE))); cand <- unique(cand[nzchar(cand) & !is.na(cand)])
    paste(intersect(cand, features), collapse = "|") }, "", USE.NAMES = FALSE)
  d$used_in_score <- d$status == "exact"
  d$note <- ifelse(d$status == "exact", "exact symbol in raw counts", "not measured under this symbol: counted as 0 in the fixed-denominator score, excluded from the available-gene diagnostic; alias / suffix columns are diagnostics only")
  d
}

#' 세포 점수. counts = 유전자 x 세포 raw dgCMatrix (모든 유전자), pm = gsx_panel_map(). 반환: library · fixed · available · n_detected (세포 x 단계) · stage_genes · gene_log · gene_det
gsx_cell_scores <- function(counts, pm, stages = NSO_STAGES) {
  lib <- Matrix::colSums(counts); if (any(lib <= 0)) stop("library 합이 0 인 세포가 있다: ", sum(lib <= 0), " 개")
  idx <- unique(pm$feature_index[pm$status == "exact"]); n <- ncol(counts)
  sub <- methods::as(counts[idx, , drop = FALSE], "CsparseMatrix"); rownames(sub) <- rownames(counts)[idx]
  Lg <- methods::as(sub %*% Matrix::Diagonal(x = 1e4 / lib), "CsparseMatrix"); Lg@x <- log1p(Lg@x); rownames(Lg) <- rownames(sub); colnames(Lg) <- colnames(counts)
  Dg <- methods::as(sub > 0, "CsparseMatrix"); Dg <- methods::as(Dg, "dMatrix"); rownames(Dg) <- rownames(sub)
  mk <- function() matrix(NA_real_, n, length(stages), dimnames = list(colnames(counts), stages)); fixed <- mk(); avail <- mk(); nd <- mk()
  sg <- do.call(rbind, lapply(stages, function(s) { z <- pm[pm$stage == s, , drop = FALSE]; g <- z$symbol[z$status == "exact"]
    data.frame(stage = s, n_listed = nrow(z), n_measured = length(g), genes_measured = paste(g, collapse = "|"), genes_missing = paste(z$symbol[z$status != "exact"], collapse = "|"),
               stringsAsFactors = FALSE) }))
  for (s in stages) { z <- pm[pm$stage == s, , drop = FALSE]; g <- z$symbol[z$status == "exact"]; if (!length(g)) next
    S <- Matrix::colSums(Lg[g, , drop = FALSE]); fixed[, s] <- S / nrow(z); avail[, s] <- S / length(g); nd[, s] <- Matrix::colSums(Dg[g, , drop = FALSE]) }
  list(library = unname(as.numeric(lib)), fixed = fixed, available = avail, n_detected = nd, stage_genes = sg, gene_log = Lg, gene_det = Dg)
}

#' 군집 점수 (긴 표): 모든 세포 평균(0 포함), donor 동일 가중(donor x 군집 >= min_cells), 검출 세포 비율. M = 세포 x 단계 주 점수, N = 검출 유전자 수, A = 관측 유전자 분모 점수(선택)
gsx_cluster_scores <- function(M, N, cluster, donor, min_cells = 5L, A = NULL) {
  cluster <- as.character(cluster); donor <- as.character(donor); u <- fdt_cluster_order_safe(cluster); f <- factor(cluster, levels = u)
  nc <- as.integer(table(f)); key <- factor(paste(cluster, donor, sep = "\r")); kc <- sub("\r.*$", "", levels(key)); kn <- as.integer(table(key))
  nd_all <- as.integer(tapply(donor, f, function(z) length(unique(z)))); nd_ok <- as.integer(vapply(u, function(k) sum(kn[kc == k] >= min_cells), 0L))
  rows <- lapply(colnames(M), function(s) { v <- M[, s]; av <- !all(is.na(v)); dv <- N[, s]
    eb <- if (av) { km <- tapply(v, key, mean); vapply(u, function(k) { w <- kc == k & kn >= min_cells; if (any(w)) mean(km[w]) else NA_real_ }, 0) } else rep(NA_real_, length(u))
    d <- data.frame(cluster = u, stage = s, n_cells = nc, n_donors = nd_all, n_donors_min_cells = nd_ok,
                    mean_all_cells = if (av) as.numeric(tapply(v, f, mean)) else NA_real_, mean_donor_equal = unname(eb),
                    frac_cells_ge1 = if (av) as.numeric(tapply(dv >= 1, f, mean)) else NA_real_, frac_cells_ge2 = if (av) as.numeric(tapply(dv >= 2, f, mean)) else NA_real_,
                    stringsAsFactors = FALSE)
    if (!is.null(A)) d$mean_available_genes <- if (av) as.numeric(tapply(A[, s], f, mean)) else NA_real_
    d })
  d <- do.call(rbind, rows); d <- d[order(match(d$cluster, u), match(d$stage, colnames(M))), , drop = FALSE]; rownames(d) <- NULL; d
}

#' 패널 유전자별 측정 · 검출 (전체 세포)
gsx_gene_detection <- function(pm, cs) {
  Lg <- cs$gene_log; Dg <- cs$gene_det; n <- length(cs$library)
  d <- pm[, c("stage", "rank", "symbol", "status", "alias_candidates_in_data", "suffix_variants_in_data")]
  m <- d$status == "exact"; g <- d$symbol
  d$n_cells_detected <- ifelse(m, as.numeric(Matrix::rowSums(Dg)[match(g, rownames(Dg))]), NA_real_); d$frac_cells_detected <- d$n_cells_detected / n
  d$mean_log_all_cells <- ifelse(m, as.numeric(Matrix::rowMeans(Lg)[match(g, rownames(Lg))]), NA_real_)
  d$availability <- ifelse(!m, "missing_gene", ifelse(d$n_cells_detected == 0, "observed_zero", "observed")); d
}
