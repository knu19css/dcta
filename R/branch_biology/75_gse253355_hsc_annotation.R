# =============================================================================
# 75_gse253355_hsc_annotation.R — GSE253355: HSC 부터 전체 세포 유형의 군집 잠정 후보 (배경 보정 module score, 두 DAG 마커 모델 따로) + 실제 R UMAP 그림 2 장
# =============================================================================
#
# 계약: Codex 인계서 codex_hsc_annotation_instruction.md (2026-10-20, 사용자 "PM이 이상하네. 근데, HSC부터의 주석을 붙여야함"). Ruling: claude_hsc_annotation_progress.md
#   - 마커: main repo node_markers_sign_review_2026-10-18.csv (sha256 고정, 읽기 전용; 결과에 원본 사본과 출처). use_in_current_panel = true · is_annotation_candidate = true ·
#     marker_sign = ±1 행만. bm_dag_broad · bm_dag_slide2 를 따로 계산하고 각각에 t_cell_broad_categories 를 더한다 (node · gene · sign 중복 제거, 같은 node 의 양 · 음
#     충돌은 보류). DAG edge · legacy stage · terminal destination 은 쓰지 않는다. 부호 = 현재 유형의 지지 / 반대 마커 (시간 증가 · 감소 아님). 검토 후보이지 검증 모델이 아니다.
#   - 점수: log1p(1e4 · count / 모든 유전자 library). 유전자 전체 세포 평균으로 24 개 같은 개수 bin (동률은 symbol 순서), 마커마다 같은 bin 의 비마커 유전자 50 개
#     (seed 20261020, 모자라면 전부). 유전자 보정값 = log 발현 − 그 대조 평균, 패널 = 측정 유전자 동일 가중. 군집마다 양성 · 음성 평균을 먼저 낸 뒤
#     score = pos − max(0, neg) (음성 패널 없으면 0, 양성 측정 유전자 없으면 NA). 0 = 미판정 대안. donor 동일 가중 = donor x 군집 >= 5 세포 donor 만.
#     p 값 · 확률 · 검증된 분류기가 아니다. bootstrap 없음.
#   - 두 모델 1 위가 같은 단일 유형이고 점수 > 0 일 때만 '잠정 후보', 아니면 '불일치 / 미판정' (두 후보 모두 표에). 동률은 이름 하나로 깨지 않는다.
#   - 기존 6 단계 점수 · 패널 · 후보는 바꾸지 않는다. 두 모델 모두 Neutrophil 인 군집에만 기존 최고 단계 후보를 보조 열로 보이고 maturity_status = unvalidated.
#   - 저자 라벨은 계산 고정 (calc_manifest.json) 뒤 표시 · 사후 대조에만 쓴다. 결과 · 공개 폴더가 비어 있지 않으면 쓰지 않는다 (exit 6).
#
# 실행 (W 에서, dcta 환경만)
#   run-r dcta R/branch_biology/75_gse253355_hsc_annotation.R --out-dir results/neutrophil_atlas_rna_2026-10-19/gse253355_hsc_annotation_01 [--public-dir …]
# exit code: 0 정상, 2 인자, 4 입력 불일치 · 검사 실패 · 계산 오류, 6 결과 폴더가 비어 있지 않음
# =============================================================================

GSH_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "23_nrv_reference.R", "24_nrv_structure.R",
                 "25_nrv_expression.R", "26_nrv_evaluate.R", "27_run_nrv.R", "34_mbv_umap_overlay.R", "36_mbv_readable_stage_zoom.R", "37_mbv_two_panel_compare.R",
                 "46_run_mbv_mono_quant.R", "50_mbv_mono_phase_view.R", "58_nso_inputs.R", "59_nso_evidence.R", "60_run_nso.R", "61_nso_view.R",
                 "62_nso_connectivity_audit.R", "63_atlas_rna_reference.R", "64_run_atlas_rna.R", "65_atlas_rna_view.R", "66_atlas_cluster_candidates.R",
                 "67_atlas_cluster_view.R", "68_gse253355_object.R", "69_gse253355_scores.R", "70_gse253355_graph.R", "71_gse253355_view.R", "72_run_gse253355.R",
                 "73_gse253355_umap.R", "74_gse253355_marker_umap.R", "75_gse253355_hsc_annotation.R")
GSH_MARKERS <- "/Users/edemhw/Honglab/CS/04_Claude/14_dcta/docs/cytoscape/rna_marker_review_2026-10-18/sign_revision_2026-10-18/node_markers_sign_review_2026-10-18.csv"   # main repo, 읽기 전용
GSH_UMAP <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_umap_03/umap_coordinates.rds"
GSH_SCORES_DIR <- "results/neutrophil_atlas_rna_2026-10-19/gse253355_02"
GSH_OBJECT <- "data/external/gse253355/S2_steady_GSE253355.local_copy.rds"
GSH_BASELINE <- "results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_hsc_annotation_baseline.json"
GSH_PUBLIC <- "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_hsc_annotation_01"
GSH_PINNED <- list(markers = "55d8f2056ffcdfb53f6dad9cf84af68d53e328ba8b82405d95b838524447f4c0", umap = "1f91bc928f06b736980bd56e18a595346ff1ddc430a96a021e7537149ab3001f",
                   lossless = "c423ecfe33be76951acfc3e71ea0173e949c908d289839a19828e48b3d6bb215", object = "83fd17576c48e710e5a03a60be485beedbc088f82c6a97fa79d7b80d5367186a",
                   calc_manifest = "5842fd04cf6e0fa943538ac671c11bbc778e4439db71662f38e4fd46871c6955")
GSH_MODELS <- c(broad = "bm_dag_broad", slide2 = "bm_dag_slide2")
GSH_TCELL <- "t_cell_broad_categories"
GSH_SEED <- 20261020L; GSH_N_BINS <- 24L; GSH_N_CTRL <- 50L; GSH_TIE_TOL <- 1e-12; GSH_MIN_DONOR_CELLS <- 5L; GSH_CHUNK <- 10000L; GSH_LABEL_MIN_FRAC <- 0.01
GSH_VERIFY <- list(c(model = "broad", type = "Neutrophil"), c(model = "slide2", type = "HSC"))                    # 세포 단위 재계산 검사 (양성만 · 양성 + 음성)
GSH_CAT_COLS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F3C8D", "#11A579", "#E73F74", "#3969AC", "#8C564B", "#80BA5A",
                  "#008695", "#CF1C90", "#F2B701", "#A5AA99", "#1F77B4", "#9467BD")
GSH_DISAGREE_COL <- "grey45"
GSH_UNDET_COL <- "grey82"


# ── 1. 마커 패널 ─────────────────────────────────────────────────────────────
gsh_read_markers <- function(path) {
  d <- utils::read.csv(path, colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = "UTF-8-BOM")
  need <- c("record_type", "network", "node_id", "display_label", "gene_symbol", "marker_sign", "use_in_current_panel", "is_annotation_candidate", "reason_ko",
            "approved_for_model", "RNA_specificity_validated")
  miss <- setdiff(need, names(d)); if (length(miss)) stop("마커 CSV 에 열이 없다: ", paste(miss, collapse = ", "))
  d$row_id <- seq_len(nrow(d)); d
}

#' 두 모델 패널: 각 DAG network + T 범주의 사용 행 (use · annotation true, 부호 ±1), node · gene · sign 중복 제거, 같은 node 양 · 음 충돌은 보류. 제외 행은 이유와 함께
gsh_panels <- function(mk, models = GSH_MODELS, tcell = GSH_TCELL) {
  sg <- suppressWarnings(as.integer(mk$marker_sign)); nets <- c(unname(models), tcell)
  ok <- mk$use_in_current_panel == "true" & mk$is_annotation_candidate == "true" & !is.na(sg) & sg %in% c(-1L, 1L) & nzchar(mk$gene_symbol) & mk$network %in% nets
  disp <- ifelse(mk$node_id == "DC", "cDC", mk$display_label)
  why <- function(i) { r <- character(0)
    if (mk$record_type[i] != "node_marker") r <- c(r, mk$record_type[i])
    if (is.na(sg[i]) || sg[i] == 0L) r <- c(r, "held: marker_sign 0")
    if (mk$use_in_current_panel[i] != "true") r <- c(r, "not in current panel")
    if (mk$is_annotation_candidate[i] != "true") r <- c(r, "not an annotation candidate")
    if (mk$node_id[i] == "Monocyte") r <- c(r, "Monocyte upper category (Classical and Nonclassical Monocyte are scored)")
    if (mk$node_id[i] == "T_cell_progenitor") r <- c(r, "T_cell_progenitor: no active annotation panel")
    if (!nzchar(mk$gene_symbol[i])) r <- c(r, "no gene symbol")
    if (!mk$network[i] %in% nets) r <- c(r, "other network")
    paste(r, collapse = "; ") }
  ex <- which(!ok)
  excluded <- data.frame(row_id = mk$row_id[ex], network = mk$network[ex], node_id = mk$node_id[ex], display_label = disp[ex], gene = mk$gene_symbol[ex], marker_sign = mk$marker_sign[ex],
                         record_type = mk$record_type[ex], use_in_current_panel = mk$use_in_current_panel[ex], is_annotation_candidate = mk$is_annotation_candidate[ex],
                         reason = vapply(ex, why, ""), reason_ko = mk$reason_ko[ex], stringsAsFactors = FALSE)
  used <- list(); held <- list(); ndup <- c()
  for (m in names(models)) { w <- which(ok & mk$network %in% c(models[[m]], tcell))
    z <- data.frame(model = m, node_id = mk$node_id[w], display_label = disp[w], gene = mk$gene_symbol[w], sign = sg[w], source_network = mk$network[w], row_id = mk$row_id[w], stringsAsFactors = FALSE)
    dup <- duplicated(paste(z$node_id, z$gene, z$sign, sep = "\r")); ndup[[m]] <- sum(dup); z <- z[!dup, , drop = FALSE]
    cf <- stats::ave(z$sign, z$node_id, z$gene, FUN = function(s) length(unique(s))) > 1
    if (any(cf)) held[[m]] <- data.frame(z[cf, c("model", "node_id", "gene", "sign", "source_network")], reason = "same node lists this gene as + and -: held", stringsAsFactors = FALSE)
    used[[m]] <- z[!cf, , drop = FALSE] }
  h0 <- data.frame(model = character(0), node_id = character(0), gene = character(0), sign = character(0), source_network = character(0), reason = character(0), stringsAsFactors = FALSE)   # 빈 표는 글자 열만 (CSV 왕복)
  out <- do.call(rbind, used); rownames(out) <- NULL
  list(used = out, excluded = excluded, held = if (length(held)) do.call(rbind, held) else h0, n_duplicates = unlist(ndup))
}


# ── 2. 유전자 집계 (원 count, 세포 덩어리 단위) ────────────────────────────────
#' log1p(1e4 · count / library) 를 세포 덩어리마다 만들어 유전자 전체 평균 · 군집 평균 · 군집 검출률 · donor x 군집 평균을 쌓는다 (세포 x 유전자 log 행렬을 한 번에 만들지 않는다)
gsh_aggregate <- function(X, lib, cluster, donor, clusters, chunk = GSH_CHUNK) {
  if (!methods::is(X, "dgCMatrix")) stop("count 행렬이 dgCMatrix 가 아니다.")
  n <- ncol(X); cluster <- as.character(cluster); donor <- as.character(donor); fc <- factor(cluster, levels = clusters); if (anyNA(fc)) stop("군집 목록에 없는 세포가 있다.")
  key <- paste(cluster, donor, sep = "\r"); gk <- sort(unique(key), method = "radix"); fk <- factor(key, levels = gk)
  G <- nrow(X); tot <- numeric(G); Cs <- matrix(0, G, length(clusters)); Ds <- Cs; Ks <- matrix(0, G, length(gk))
  for (s in seq(1L, n, by = chunk)) { idx <- s:min(n, s + chunk - 1L); sub <- X[, idx, drop = FALSE]
    sub@x <- log1p(1e4 * sub@x / rep(lib[idx], diff(sub@p))); det <- sub; det@x <- as.numeric(det@x > 0)
    Ic <- Matrix::sparseMatrix(i = seq_along(idx), j = as.integer(fc[idx]), x = 1, dims = c(length(idx), length(clusters)))
    Ik <- Matrix::sparseMatrix(i = seq_along(idx), j = as.integer(fk[idx]), x = 1, dims = c(length(idx), length(gk)))
    tot <- tot + Matrix::rowSums(sub); Cs <- Cs + as.matrix(sub %*% Ic); Ds <- Ds + as.matrix(det %*% Ic); Ks <- Ks + as.matrix(sub %*% Ik) }
  nc <- tabulate(fc, length(clusters)); nk <- tabulate(fk, length(gk))
  C <- sweep(Cs, 2, nc, "/"); D <- sweep(Ds, 2, nc, "/"); K <- sweep(Ks, 2, nk, "/"); dimnames(C) <- dimnames(D) <- list(rownames(X), clusters); dimnames(K) <- list(rownames(X), gk)
  list(gene_mean = stats::setNames(tot / n, rownames(X)), C = C, D = D, K = K, n_cells = stats::setNames(nc, clusters),
       groups = data.frame(group = gk, cluster = sub("\r.*$", "", gk), donor = sub("^.*\r", "", gk), n_cells = nk, stringsAsFactors = FALSE))
}


# ── 3. 배경 대조 유전자 ──────────────────────────────────────────────────────
#' 유전자 전체 평균으로 같은 개수 bin (동률은 symbol C locale 순서)
gsh_bins <- function(gene_mean, n_bins = GSH_N_BINS) {
  g <- names(gene_mean); o <- order(gene_mean, g, method = "radix"); b <- integer(length(g)); b[o] <- as.integer(ceiling(seq_along(o) * n_bins / length(o))); stats::setNames(b, g)
}

#' 마커 유전자마다 같은 bin 의 비마커 유전자 n_ctrl 개 (seed 한 번 고정, 마커 순서 = symbol C locale 순서; 모자라면 bin 의 비마커 전부)
gsh_controls <- function(bins, markers, exclude, n_ctrl = GSH_N_CTRL, seed = GSH_SEED) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  mk <- sort(unique(markers), method = "radix"); g <- names(bins); out <- list()
  tab <- do.call(rbind, lapply(mk, function(m) { cand <- sort(g[bins == bins[[m]] & !g %in% exclude], method = "radix")
    pick <- if (length(cand) > n_ctrl) cand[sort(sample.int(length(cand), n_ctrl))] else cand; out[[m]] <<- pick
    data.frame(gene = m, bin = unname(bins[[m]]), n_candidates_in_bin = length(cand), n_controls = length(pick), shortfall = n_ctrl - length(pick), controls = paste(pick, collapse = "|"),
               stringsAsFactors = FALSE) }))
  list(controls = out, table = tab)
}


# ── 4. 유형 점수 · 후보 ──────────────────────────────────────────────────────
#' 한 모델의 군집 x 유형 점수 (긴 표). agg = gsh_aggregate(), panel = 그 모델의 사용 행, ctrl = 마커별 대조 유전자
gsh_type_scores <- function(agg, panel, ctrl, min_cells = GSH_MIN_DONOR_CELLS) {
  meas <- unique(panel$gene); meas <- meas[meas %in% rownames(agg$C) & vapply(meas, function(g) length(ctrl[[g]]) > 0, TRUE)]
  em <- function(M) { cm <- matrix(vapply(meas, function(g) colMeans(M[ctrl[[g]], , drop = FALSE]), numeric(ncol(M))), nrow = ncol(M)); M[meas, , drop = FALSE] - t(cm) }
  EC <- em(agg$C); EK <- em(agg$K); grp <- agg$groups; okg <- grp$n_cells >= min_cells; cl <- colnames(agg$C); types <- unique(panel$node_id)
  cm2 <- function(E, g, k) if (length(g)) colMeans(E[g, , drop = FALSE]) else rep(NA_real_, k)
  d <- do.call(rbind, lapply(types, function(ty) { z <- panel[panel$node_id == ty, , drop = FALSE]
    pg <- intersect(z$gene[z$sign == 1L], meas); ng <- intersect(z$gene[z$sign == -1L], meas)
    pos <- cm2(EC, pg, length(cl)); neg <- cm2(EC, ng, length(cl)); pen <- if (length(ng)) pmax(0, neg) else rep(0, length(cl))
    kp <- cm2(EK, pg, ncol(EK)); kn <- cm2(EK, ng, ncol(EK))
    dep <- vapply(cl, function(k) { w <- grp$cluster == k & okg; if (any(w)) mean(kp[w]) else NA_real_ }, 0)
    den <- vapply(cl, function(k) { w <- grp$cluster == k & okg; if (any(w) && length(ng)) mean(kn[w]) else NA_real_ }, 0)
    det <- if (length(pg)) colMeans(agg$D[pg, , drop = FALSE]) else rep(NA_real_, length(cl))
    top <- vapply(cl, function(k) { if (!length(pg)) return(""); v <- EC[pg, k] / length(pg); o <- utils::head(order(-v, pg, method = "radix"), 3L)
      paste(sprintf("%s %+.3f", pg[o], v[o]), collapse = " | ") }, "")
    data.frame(cluster = cl, n_cells = as.integer(agg$n_cells[cl]), type = ty, display_label = z$display_label[1], n_pos_listed = sum(z$sign == 1L), n_pos_measured = length(pg),
               n_neg_listed = sum(z$sign == -1L), n_neg_measured = length(ng), missing_markers = paste(setdiff(z$gene, meas), collapse = "|"),
               pos_enrichment = unname(pos), neg_enrichment = unname(neg), penalty = unname(pen), score = unname(pos - pen), mean_pos_detection = unname(det),
               top_pos_contributors = unname(top), n_donors_eligible = vapply(cl, function(k) sum(grp$cluster == k & okg), 0L, USE.NAMES = FALSE),
               donor_equal_pos = unname(dep), donor_equal_neg = unname(den), donor_equal_score = unname(dep - if (length(ng)) pmax(0, den) else 0), stringsAsFactors = FALSE) }))
  d$rank_in_cluster <- as.integer(stats::ave(seq_len(nrow(d)), d$cluster, FUN = function(ix) { o <- order(-d$score[ix], match(d$type[ix], types), na.last = TRUE)
    r <- integer(length(ix)); r[o] <- seq_along(ix); r }))
  rownames(d) <- NULL; d
}

#' 군집마다 1 위 (동률 모두) · 2 위 (다음 유형과 미판정 0 중 큰 것) · 차이 · 상태. 1 위 <= 0 이면 미판정 (0 = 미판정 대안)
gsh_rank <- function(sc, col = "score", tol = GSH_TIE_TOL) {
  do.call(rbind, lapply(unique(sc$cluster), function(k) { z <- sc[sc$cluster == k & !is.na(sc[[col]]), , drop = FALSE]
    mk <- function(st, t1, l1, s1, t2, l2, s2, nt) data.frame(cluster = k, status = st, top1 = t1, top1_label = l1, top1_score = s1, top2 = t2, top2_label = l2, top2_score = s2,
                                                       gap = s1 - s2, n_top_tied = nt, stringsAsFactors = FALSE)
    if (!nrow(z)) return(mk("undetermined", "", "undetermined (no score)", NA_real_, "", "", NA_real_, 0L))
    z <- z[order(-z[[col]], seq_len(nrow(z))), , drop = FALSE]; v <- z[[col]]; top <- v[1]
    if (top <= 0) return(mk("undetermined", "", "undetermined", 0, z$type[1], z$display_label[1], top, 0L))
    tied <- v >= top - tol; nt <- sum(tied); nx <- if (nt < nrow(z)) v[nt + 1L] else -Inf
    if (nx > 0) mk(if (nt > 1L) "tied" else "candidate", paste(z$type[tied], collapse = "|"), paste(z$display_label[tied], collapse = " / "), top, z$type[nt + 1L], z$display_label[nt + 1L], nx, nt)
    else mk(if (nt > 1L) "tied" else "candidate", paste(z$type[tied], collapse = "|"), paste(z$display_label[tied], collapse = " / "), top, "", "undetermined", 0, nt) }))
}

#' 두 모델 합의: 둘 다 단일 1 위 · 같은 유형 · 점수 > 0 → provisional, 둘 다 미판정 → undetermined, 그 밖 → disagreement
gsh_consensus <- function(ra, rb) {
  st <- ifelse(ra$status == "candidate" & rb$status == "candidate" & ra$top1 == rb$top1, "provisional",
               ifelse(ra$status == "undetermined" & rb$status == "undetermined", "undetermined", "disagreement"))
  data.frame(cluster = ra$cluster, consensus_status = st, consensus_type = ifelse(st == "provisional", ra$top1, ""),
             consensus_label = ifelse(st == "provisional", ra$top1_label, ifelse(st == "undetermined", "undetermined", paste(ra$top1_label, rb$top1_label, sep = " vs "))),
             stringsAsFactors = FALSE)
}

#' 유형마다 점수가 가장 높은 군집 (동률 모두) · 2 위 군집 · 차이 (단계별 peak 와 같은 비교, 군집 안 1 위와 다르다)
gsh_type_peaks <- function(sc, tol = GSH_TIE_TOL) {
  key <- unique(sc[, c("model", "type", "display_label")])
  do.call(rbind, lapply(seq_len(nrow(key)), function(i) { z <- sc[sc$model == key$model[i] & sc$type == key$type[i] & !is.na(sc$score), , drop = FALSE]
    if (!nrow(z)) return(data.frame(key[i, ], peak_clusters = "", peak_score = NA_real_, second_cluster = "", second_score = NA_real_, gap = NA_real_, peak_positive = FALSE, stringsAsFactors = FALSE))
    z <- z[order(-z$score, seq_len(nrow(z))), , drop = FALSE]; tie <- z$score >= z$score[1] - tol; nt <- sum(tie)
    data.frame(key[i, ], peak_clusters = paste(z$cluster[tie], collapse = "|"), peak_score = z$score[1], second_cluster = if (nt < nrow(z)) z$cluster[nt + 1L] else "",
               second_score = if (nt < nrow(z)) z$score[nt + 1L] else NA_real_, gap = if (nt < nrow(z)) z$score[1] - z$score[nt + 1L] else NA_real_, peak_positive = z$score[1] > 0,
               stringsAsFactors = FALSE) }))
}

#' 세포 단위 재계산 (검사): 한 유형의 세포별 양성 · 음성 보정값을 세포 덩어리로 다시 계산해 군집 평균을 낸다
gsh_cell_level <- function(X, lib, cluster, clusters, genes_pos, genes_neg, ctrl, chunk = GSH_CHUNK) {
  mg <- unique(c(genes_pos, genes_neg)); cg <- unique(unlist(ctrl[mg])); ug <- unique(c(mg, cg)); Xs <- X[ug, , drop = FALSE]; n <- ncol(X)
  W <- Matrix::sparseMatrix(i = rep(seq_along(mg), lengths(ctrl[mg])), j = match(unlist(ctrl[mg]), ug), x = rep(1 / lengths(ctrl[mg]), lengths(ctrl[mg])), dims = c(length(mg), length(ug)))
  pos <- numeric(n); neg <- numeric(n)
  for (s in seq(1L, n, by = chunk)) { idx <- s:min(n, s + chunk - 1L); sub <- Xs[, idx, drop = FALSE]; sub@x <- log1p(1e4 * sub@x / rep(lib[idx], diff(sub@p)))
    E <- as.matrix(sub[match(mg, ug), , drop = FALSE] - W %*% sub); rownames(E) <- mg
    pos[idx] <- colMeans(E[genes_pos, , drop = FALSE]); neg[idx] <- if (length(genes_neg)) colMeans(E[genes_neg, , drop = FALSE]) else NA_real_ }
  f <- factor(as.character(cluster), levels = clusters)
  data.frame(cluster = clusters, pos_cell_mean = as.numeric(tapply(pos, f, mean)), neg_cell_mean = if (length(genes_neg)) as.numeric(tapply(neg, f, mean)) else NA_real_, stringsAsFactors = FALSE)
}


# ── 5. 그림 (R base graphics, R73 저자 색 · R74 배치) ───────────────────────────
#' 그림 1 · 2: 왼쪽 저자 주석, 오른쪽 군집 잠정 후보 (잠정 = 유형 색, 불일치 = 짙은 회색, 미판정 = 옅은 회색) + 번호. mask = NULL 이면 전체 세포.
#' 번호: 전체 = 모든 군집 번호, 부분집합 = 부분집합 세포 label_min_frac 이상 군집에 "번호 · 이름" 크게. 반환: 글자 표 · 그린 좌표 · 번호 군집
gsh_draw_pair <- function(path, U, author, acols, cluster, cons, ccol, mask = NULL, subset_labels = NULL, label_min_frac = 0, sub = "", ko = TRUE, fam = "") {
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note; cluster <- as.character(cluster); full <- is.null(mask); keep <- if (full) rep(TRUE, nrow(U)) else as.logical(mask)
  Us <- U[keep, , drop = FALSE]; au <- as.character(author[keep]); cl <- cluster[keep]; ns <- nrow(Us); fmt <- .nso_fmt_n; pad <- function(r) r + c(-1, 1) * 0.03 * diff(r)
  b <- if (full) mbvu_bounds(U) else list(xlim = pad(range(Us[, 1])), ylim = pad(range(Us[, 2])))
  ck <- stats::setNames(ifelse(cons$consensus_status == "provisional", cons$consensus_label, ifelse(cons$consensus_status == "disagreement", "disagreement", "undetermined")), cons$cluster)
  if (!all(cl %in% names(ck))) stop("합의 표에 없는 군집의 세포가 있다."); pcol <- unname(ccol[ck[cl]])
  nm_ko <- function(k) if (k == "disagreement") (if (ko) "불일치" else "disagreement") else if (k == "undetermined") (if (ko) "미판정" else "undetermined") else k
  cnt <- table(factor(cl, levels = cons$cluster)); shown <- names(cnt)[cnt > 0]
  numbered <- if (full) shown else shown[as.integer(cnt[shown]) >= label_min_frac * ns]; unn <- setdiff(shown, numbered)
  short <- function(k) { i <- match(k, cons$cluster); if (cons$consensus_status[i] == "provisional") cons$consensus_label[i]
    else if (cons$consensus_status[i] == "undetermined") nm_ko("undetermined") else sprintf("%s (%s)", nm_ko("disagreement"), gsub(" vs ", " / ", cons$consensus_label[i], fixed = TRUE)) }
  pcex <- if (full) 0.28 else 0.42
  grDevices::png(path, width = 4400, height = 2100, units = "px", res = 100, pointsize = 12)
  tryCatch({
    graphics::layout(matrix(1:4, 1), widths = c(3, 1.4, 3, 1.4)); graphics::par(oma = c(5.4, 0, 7.6, 0), family = fam)
    graphics::par(mar = c(0.5, 0.5, 5.6, 0.5)); mbvu_open_panel(b); o <- .gsu_order(au)
    graphics::points(Us[o, 1], Us[o, 2], pch = 16, cex = pcex, col = acols[au[o]]); graphics::box(col = "grey60")
    if (!full) { rp <- atlas_cl_rep_cells(Us, au); cx <- 1.7; lp <- .acl_label_layout(rp, cx, b, rp$n_cells); .acl_draw_labels(lp, lp$text, cx); for (k in lp$text) note("author", "label", k) }
    t1 <- if (ko) "저자 주석 (cluster_anno_l2)" else "Author annotation (cluster_anno_l2)"
    s1 <- sprintf(if (ko) "라벨 %d 개 · 사후 대조용 (정답 아님)" else "%d labels · post hoc comparison (not truth)", length(unique(au)))
    graphics::title(main = t1, cex.main = 2.15, line = 2.9); graphics::mtext(s1, side = 3, line = 0.8, cex = 1.4, col = "grey20"); note("author", "title", t1); note("author", "sub", s1)
    graphics::par(mar = c(0.5, 0, 5.6, 0.5)); graphics::plot.new(); labs <- if (full) names(acols)[names(acols) %in% au] else subset_labels; tb <- table(au)
    lg <- sprintf("%s (%s)", labs, fmt(as.integer(tb[labs])))
    graphics::legend("left", legend = lg, pch = 16, col = acols[labs], pt.cex = if (full) 2.4 else 3.0, cex = if (full) 1.7 else 1.95, bty = "n", y.intersp = if (full) 1.0 else 1.3,
                     title = if (ko) "저자 라벨 (세포 수)" else "author label (cells)", title.adj = 0)
    for (k in lg) note("author_legend", "item", k)
    graphics::par(mar = c(0.5, 0.5, 5.6, 0.5)); mbvu_open_panel(b); o2 <- .gsu_order(cl)
    graphics::points(Us[o2, 1], Us[o2, 2], pch = 16, cex = pcex, col = pcol[o2]); graphics::box(col = "grey60")
    rc <- atlas_cl_rep_cells(Us, cl); rc <- rc[rc$cluster %in% numbered, , drop = FALSE]; if (!full) rc$cluster <- vapply(rc$cluster, function(k) sprintf("%s · %s", k, short(k)), "")
    cx2 <- if (full) 1.8 else 2.0; lc <- .acl_label_layout(rc, cx2, b, rc$n_cells); .acl_draw_labels(lc, lc$text, cx2)
    for (k in lc$text) note("candidate", "label", k); if (any(!lc$free)) note("candidate", "labels_overlapping", paste(lc$text[!lc$free], collapse = " | "))
    np <- sum(cons$consensus_status[match(shown, cons$cluster)] == "provisional"); nd <- sum(cons$consensus_status[match(shown, cons$cluster)] == "disagreement")
    nu <- sum(cons$consensus_status[match(shown, cons$cluster)] == "undetermined")
    t2 <- if (ko) "군집 잠정 후보 (두 마커 모델이 같을 때만, 검증 전)" else "Provisional cluster candidates (only when both marker models agree; unvalidated)"
    s2 <- sprintf(if (ko) "잠정 %d · 불일치 %d · 미판정 %d 군집%s · 성숙 단계 (PM 등) 는 확정하지 않음" else "provisional %d · disagreement %d · undetermined %d clusters%s · maturity stage (e.g. PM) not assigned",
                  np, nd, nu, if (full) "" else sprintf(if (ko) " · 번호 = 이 세포의 %s%% 이상 군집 %d 개, 나머지 %d 군집 %s 세포는 색만" else " · numbered = %d clusters >= %s%%, other %d clusters (%s cells) colour only",
                                                          if (ko) format(100 * label_min_frac) else length(numbered), if (ko) length(numbered) else format(100 * label_min_frac), length(unn), fmt(sum(cnt[unn]))))
    graphics::title(main = t2, cex.main = 2.15, line = 2.9); graphics::mtext(s2, side = 3, line = 0.8, cex = if (full) 1.4 else 1.2, col = "grey20"); note("candidate", "title", t2); note("candidate", "sub", s2)
    graphics::par(mar = c(0.5, 0, 5.6, 0.5)); graphics::plot.new()
    keys <- unique(ck[shown]); keys <- c(setdiff(keys, c("disagreement", "undetermined")), intersect(c("disagreement", "undetermined"), keys))
    items <- vapply(keys, function(k) { kc <- shown[ck[shown] == k]; nn <- sum(cnt[kc]); ids <- if (full) kc else intersect(kc, numbered)
      idt <- if (!length(ids)) (if (ko) "번호 없음" else "unnumbered") else if (length(ids) <= 4L) paste0(if (ko) "군집 " else "clusters ", paste(ids, collapse = ","))
             else sprintf(if (ko) "%d 군집" else "%d clusters", length(ids))
      sprintf(if (ko) "%s  %s · %s 세포" else "%s  %s · %s cells", nm_ko(k), idt, fmt(nn)) }, "")
    graphics::legend("left", legend = items, pch = 15, col = ccol[keys], pt.cex = 3.6, cex = 1.6, bty = "n", y.intersp = 1.45, title = if (ko) "잠정 후보 (군집 · 세포)" else "provisional candidate (clusters · cells)", title.adj = 0)
    for (k in items) note("candidate_legend", "item", k)
    main <- if (full) sprintf(if (ko) "GSE253355 골수 전체 %s 세포 — 저자 주석과 HSC 부터 전체 세포 유형 군집 잠정 후보" else "GSE253355 marrow, all %s cells - author annotation vs provisional cell-type candidates from HSC", fmt(ns))
            else sprintf(if (ko) "HSC 부터 Neutrophil · Monocyte 까지 저자 라벨 %d 개 (%s 세포) — 같은 UMAP, 잠정 후보" else "Author labels from HSC to Neutrophil / Monocyte (%d labels, %s cells) - same UMAP, provisional candidates", length(subset_labels), fmt(ns))
    graphics::mtext(main, side = 3, outer = TRUE, line = 4.5, cex = 2.5, font = 2); graphics::mtext(sub, side = 3, outer = TRUE, line = 2.0, cex = 1.4, col = "grey20")
    note("figure", "title", main); note("figure", "sub", sub)
    foot <- if (ko) c("점수 = 배경 보정 module score (마커 log1p(CP10k) - 같은 발현 bin 비마커 50 유전자 평균), 군집 평균 양성 - max(0, 음성). 0 이하 = 미판정. 확률 · 검증된 분류가 아니다.",   # ASCII '-' (한글 글꼴에 U+2212 없음, hsc_annotation_01 검토)
                      "bm_dag_broad · bm_dag_slide2 (각각 + T 범주) 를 따로 계산해 1 위가 같은 단일 유형일 때만 잠정 후보. MSC · 내피 등은 패널이 없어 미판정이 이들을 배제한다는 뜻이 아니다.")
            else c("Score = background-corrected module score (marker log1p(CP10k) minus the mean of 50 same-bin non-marker genes), cluster mean positive - max(0, negative). <= 0 = undetermined. Not a probability or validated classifier.",
                   "bm_dag_broad and bm_dag_slide2 (each + T categories) scored separately; provisional only when both give the same single top. No MSC / endothelial panels, so undetermined does not exclude them.")
    for (k in seq_along(foot)) { graphics::mtext(foot[k], side = 1, outer = TRUE, line = 0.9 + (k - 1) * 1.75, cex = 1.3, col = "grey20"); note("footer", "footer", foot[k]) }
  }, finally = grDevices::dev.off())
  list(text = nf$table(), coords = Us, numbered = numbered, unnumbered = data.frame(cluster = unn, n_cells = as.integer(cnt[unn]), stringsAsFactors = FALSE),
       colours = data.frame(cluster = shown, key = unname(ck[shown]), colour = vapply(shown, function(k) pcol[cl == k][1], "", USE.NAMES = FALSE),
                            n_colours = vapply(shown, function(k) length(unique(pcol[cl == k])), 0L, USE.NAMES = FALSE), stringsAsFactors = FALSE))
}


# ── 6. 실행기 ────────────────────────────────────────────────────────────────
gsh_main <- function(argv, root = NULL, pins = GSH_PINNED) {
  keys <- c("--markers" = "markers", "--umap" = "umap", "--scores-dir" = "scores_dir", "--object" = "object", "--baseline" = "baseline", "--out-dir" = "out", "--public-dir" = "public",
            "--assay" = "assay"); a <- list(); i <- 1L
  while (i <= length(argv)) { k <- argv[i]; v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (!k %in% names(keys)) { cat("알 수 없는 인자:", k, "\n"); return(2L) }
    a[[keys[[k]]]] <- v; i <- i + 2L }
  if (is.null(a$out) || is.na(a$out) || !nzchar(a$out)) { cat("--out-dir 가 필요하다.\n"); return(2L) }
  if (is.null(root)) root <- mdu_find_root()
  if (is.null(a$markers)) a$markers <- GSH_MARKERS
  dflt <- list(umap = GSH_UMAP, scores_dir = GSH_SCORES_DIR, object = GSH_OBJECT, baseline = GSH_BASELINE, public = GSH_PUBLIC)
  for (k in names(dflt)) if (is.null(a[[k]])) a[[k]] <- file.path(root, dflt[[k]])
  if (is.null(a$assay)) a$assay <- "RNA"
  for (d in c(a$out, a$public)) if (dir.exists(d) && length(list.files(d, all.files = TRUE, no.. = TRUE))) { cat("결과 폴더가 비어 있지 않다(덮어쓰지 않는다):", d, "\n"); return(6L) }
  out <- a$out; dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), basename(out), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", GSH_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", length(GSH_SOURCES), "files"); ld <- mdu_load_snapshot(snap, files = GSH_SOURCES)
    get("gsh_build", envir = ld$env, inherits = FALSE)(a = a, out = out, snap = snap, log = log, argv = argv, pins = pins, root = root)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}

gsh_build <- function(a, out, snap, log, argv, pins, root) {
  ck <- .mbvq_checks(); t0 <- Sys.time(); stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"); rel <- function(x) substring(normalizePath(x), nchar(normalizePath(out)) + 2L)
  files <- character(0); crt <- logical(0); fmt <- .nso_fmt_n
  wc <- function(df, f) { p0 <- file.path(out, f); dir.create(dirname(p0), recursive = TRUE, showWarnings = FALSE); crt[[f]] <<- atlas_write_csv(df, p0); files <<- c(files, p0); invisible(p0) }
  rdc <- function(p) utils::read.csv(p, colClasses = "character", na.strings = character(0), check.names = FALSE)
  sd <- a$scores_dir; f_ls <- file.path(sd, "lossless.rds"); f_cm <- file.path(sd, "calc_manifest.json"); f_cc <- file.path(sd, "cluster_candidates.csv")
  # ── (1) 입력 고정 · 기준선 ──
  if (grepl("\\.part$", a$object)) stop("다운로드 중인 .part 입력은 쓰지 않는다: ", a$object)
  for (p in c(a$markers, a$umap, f_ls, f_cm, f_cc, a$object, a$baseline)) if (!file.exists(p)) stop("입력 파일이 없다: ", p)
  for (k in c("markers", "umap", "lossless", "object", "calc_manifest")) if (is.null(pins[[k]])) stop("hash 고정값이 없다: ", k)
  h0 <- c(markers = .nrv_sha(a$markers), umap = .nrv_sha(a$umap), lossless = .nrv_sha(f_ls), object = .nrv_sha(a$object), calc_manifest = .nrv_sha(f_cm), cluster_candidates = .nrv_sha(f_cc),
          baseline = .nrv_sha(a$baseline))
  for (k in c("markers", "umap", "lossless", "object", "calc_manifest")) if (!identical(h0[[k]], pins[[k]])) stop(k, " hash 가 고정값과 다르다: ", h0[[k]])
  if (!identical(h0[["cluster_candidates"]], jsonlite::fromJSON(f_cm)$files[["cluster_candidates.csv"]])) stop("cluster_candidates.csv hash 가 calc_manifest.json 기록과 다르다.")
  bl <- jsonlite::fromJSON(a$baseline); blf <- function() vapply(names(bl), function(f) { p <- file.path(root, f); if (file.exists(p)) .nrv_sha(p) else "missing" }, "")
  bl0 <- blf(); if (!identical(unname(bl0), unname(unlist(bl)))) stop("기준선 파일 hash 가 codex_hsc_annotation_baseline.json 과 다르다: ", paste(names(bl)[bl0 != unlist(bl)], collapse = ", "))
  ck$add("inputs_pinned", TRUE, sprintf("markers %s, umap %s, lossless %s, object %s, calc_manifest %s == pinned; cluster_candidates.csv == calc_manifest.json record",
                                        substr(h0[["markers"]], 1, 12), substr(h0[["umap"]], 1, 12), substr(h0[["lossless"]], 1, 12), substr(h0[["object"]], 1, 12), substr(h0[["calc_manifest"]], 1, 12)))
  # ── (2) 마커 패널 (원본 사본 · 출처) ──
  mk <- gsh_read_markers(a$markers); pn <- gsh_panels(mk); pu <- pn$used
  mdir <- file.path(out, "marker_input"); dir.create(mdir, showWarnings = FALSE); mcp <- file.path(mdir, basename(a$markers))
  if (!file.copy(a$markers, mcp, overwrite = FALSE) || !identical(.nrv_sha(mcp), h0[["markers"]])) stop("마커 CSV 사본이 원본과 다르다."); Sys.chmod(mcp, "0444"); files <- c(files, mcp)
  wc(data.frame(item = c("source_path", "sha256", "rows", "approved_for_model", "RNA_specificity_validated", "note"),
                value = c(a$markers, h0[["markers"]], as.character(nrow(mk)), paste(unique(mk$approved_for_model), collapse = "|"), paste(unique(mk$RNA_specificity_validated), collapse = "|"),
                          "review candidates, not a validated model; copied read-only without edits"), stringsAsFactors = FALSE), "marker_input/source.csv")
  # ── (3) 고정 좌표 · 군집 · 원 count (저자 라벨은 계산 뒤에만) ──
  um <- readRDS(a$umap); ls <- readRDS(f_ls); cells <- ls$cells; cl <- as.character(ls$clusters); donor <- as.character(ls$donor); n <- length(cells)
  U <- um$U; mask <- um$subset_mask; sl <- as.character(um$subset_labels); clusters <- unique(rdc(f_cc)$cluster)
  ro <- gsx_open_rds(a$object); obj <- ro$object; md <- gsx_meta(obj); cn <- gsx_counts(obj, a$assay); X <- cn$counts; rm(obj, ro); invisible(gc())
  ok_ids <- identical(um$cells, cells) && identical(as.character(um$new_cluster), cl) && identical(colnames(X), cells) && identical(rownames(md), cells) && length(clusters) == 26L &&
            setequal(clusters, unique(cl)) && identical(dim(U), c(n, 2L)) && identical(rownames(U), cells) && all(is.finite(U)) && length(mask) == n && length(donor) == n
  if (!ok_ids) stop("세포 ID · 순서 · 군집 · 좌표가 입력끼리 맞지 않는다.")
  ck$add("cells_clusters_unchanged", TRUE, sprintf("%d cell IDs identical in order (UMAP RDS, lossless, raw counts, meta.data); %d clusters == lossless clusters == cluster_candidates.csv; U %d x 2 finite, drawn as stored", n, length(clusters), n))
  lib <- Matrix::colSums(X); if (!identical(unname(as.numeric(lib)), ls$library)) stop("library 가 lossless 와 다르다.")
  ck$add("library_matches_lossless", TRUE, sprintf("all-gene library (colSums of %d genes x %d cells raw counts) identical to lossless.rds", nrow(X), ncol(X)))
  rm(md); invisible(gc()); log("inputs fixed:", n, "cells,", length(clusters), "clusters,", nrow(X), "genes | panels:", paste(names(GSH_MODELS), table(pu$model)[names(GSH_MODELS)], collapse = " "), "rows")
  # ── (4) 집계 · bin · 대조 · 점수 (저자 라벨 없음) ──
  agg <- gsh_aggregate(X, lib, cl, donor, clusters); bins <- gsh_bins(agg$gene_mean)
  allm <- unique(mk$gene_symbol[mk$use_in_current_panel == "true" & mk$is_annotation_candidate == "true" & mk$marker_sign %in% c("1", "-1") & nzchar(mk$gene_symbol)])
  meas_m <- intersect(unique(pu$gene), names(bins)); ct <- gsh_controls(bins, meas_m, exclude = allm)
  log("aggregated:", length(bins), "genes in", GSH_N_BINS, "bins; controls for", length(meas_m), "measured marker genes")
  sc <- do.call(rbind, lapply(names(GSH_MODELS), function(m) data.frame(model = m, gsh_type_scores(agg, pu[pu$model == m, , drop = FALSE], ct$controls), stringsAsFactors = FALSE)))
  rk <- lapply(stats::setNames(names(GSH_MODELS), names(GSH_MODELS)), function(m) gsh_rank(sc[sc$model == m, , drop = FALSE], "score"))
  rd <- lapply(stats::setNames(names(GSH_MODELS), names(GSH_MODELS)), function(m) gsh_rank(sc[sc$model == m, , drop = FALSE], "donor_equal_score"))
  cons <- gsh_consensus(rk[[1]], rk[[2]]); pk <- gsh_type_peaks(sc)
  # 세포 단위 재계산 (선형 집계 검사)
  vr <- do.call(rbind, lapply(GSH_VERIFY, function(v) { z <- pu[pu$model == v[["model"]] & pu$node_id == v[["type"]], , drop = FALSE]
    gp <- intersect(z$gene[z$sign == 1L], meas_m); gn <- intersect(z$gene[z$sign == -1L], meas_m); c0 <- gsh_cell_level(X, lib, cl, clusters, gp, gn, ct$controls)
    s0 <- sc[sc$model == v[["model"]] & sc$type == v[["type"]], , drop = FALSE]; s0 <- s0[match(clusters, s0$cluster), ]
    data.frame(model = v[["model"]], type = v[["type"]], cluster = clusters, pos_aggregated = s0$pos_enrichment, pos_cell_level = c0$pos_cell_mean, neg_aggregated = s0$neg_enrichment,
               neg_cell_level = c0$neg_cell_mean, stringsAsFactors = FALSE) }))
  vd <- max(abs(c(vr$pos_aggregated - vr$pos_cell_level, vr$neg_aggregated - vr$neg_cell_level)), na.rm = TRUE)
  rm(X); invisible(gc())
  # 검사: 패널 · 대조 · 표
  pchk <- vapply(names(GSH_MODELS), function(m) { z <- pu[pu$model == m, ]; !anyDuplicated(paste(z$node_id, z$gene, z$sign)) && all(z$sign %in% c(-1L, 1L)) &&
    !any(stats::ave(z$sign, z$node_id, z$gene, FUN = function(s) length(unique(s))) > 1) && !any(z$row_id %in% pn$excluded$row_id) }, TRUE)
  ck$add("panels_split_dedup_held", all(pchk) && all(nzchar(pn$excluded$reason)),
         sprintf("%s; duplicates removed %s; +/- conflicts held %d; %d excluded rows with reasons (held sign 0, empty panels, unresolved, Monocyte upper category, T_cell_progenitor)",
                 paste(vapply(names(GSH_MODELS), function(m) { z <- pu[pu$model == m, ]; sprintf("%s %d types (+%d / -%d rows)", m, length(unique(z$node_id)), sum(z$sign == 1L), sum(z$sign == -1L)) }, ""), collapse = "; "),
                 paste(names(pn$n_duplicates), pn$n_duplicates, sep = "=", collapse = " "), nrow(pn$held), nrow(pn$excluded)))
  cg <- unlist(ct$controls); bsz <- table(bins)
  ok_ct <- !any(cg %in% allm) && all(vapply(names(ct$controls), function(g) all(bins[ct$controls[[g]]] == bins[[g]]), TRUE)) && all(ct$table$n_controls == pmin(GSH_N_CTRL, ct$table$n_candidates_in_bin)) &&
           length(bsz) == GSH_N_BINS && (max(bsz) - min(bsz)) <= 1L
  ck$add("controls_exclude_markers", ok_ct, sprintf("%d bins of %d-%d genes; %d marker genes x up to %d same-bin controls (seed %d); no control is any of %d listed marker genes; shortfall total %d",
                                                    length(bsz), min(bsz), max(bsz), nrow(ct$table), GSH_N_CTRL, GSH_SEED, length(allm), sum(ct$table$shortfall)))
  npos <- sc$n_pos_measured > 0; nty <- tapply(sc$type, sc$model, function(z) length(unique(z)))
  ok_tab <- all(vapply(names(GSH_MODELS), function(m) sum(sc$model == m) == length(clusters) * nty[[m]], TRUE)) && all(is.finite(sc$score[npos])) && all(is.na(sc$score[!npos])) &&
            all(tapply(sc$n_cells, sc$model, function(z) sum(z)) == n * nty[names(GSH_MODELS)]) && nrow(cons) == length(clusters)
  ck$add("score_tables_complete", ok_tab, sprintf("%s rows = %d clusters x types; scores finite where >= 1 positive marker measured (NA otherwise: %d rows); cluster sizes sum to %d cells per type",
                                                  paste(sprintf("%s %d", names(GSH_MODELS), as.integer(table(sc$model)[names(GSH_MODELS)])), collapse = ", "), length(clusters), sum(!npos), n))
  ck$add("aggregate_equals_cell_level", is.finite(vd) && vd <= 1e-10, sprintf("cluster means of per-cell enrichments (%s) == aggregated values; max |diff| %.3g <= 1e-10",
                                                                              paste(vapply(GSH_VERIFY, function(v) paste(v, collapse = ":"), ""), collapse = ", "), vd))
  # 계산 표 고정
  cand <- data.frame(cluster = clusters, n_cells = as.integer(agg$n_cells[clusters]), n_subset_cells = as.integer(table(factor(cl[mask], levels = clusters))), stringsAsFactors = FALSE)
  for (m in names(GSH_MODELS)) { r <- rk[[m]]; e <- rd[[m]]
    cand[[paste0(m, "_status")]] <- r$status; cand[[paste0(m, "_top1")]] <- r$top1_label; cand[[paste0(m, "_top1_score")]] <- r$top1_score; cand[[paste0(m, "_top2")]] <- r$top2_label
    cand[[paste0(m, "_top2_score")]] <- r$top2_score; cand[[paste0(m, "_gap")]] <- r$gap; cand[[paste0(m, "_donor_equal_top1")]] <- e$top1_label; cand[[paste0(m, "_donor_equal_status")]] <- e$status
    cand[[paste0(m, "_donor_equal_agrees")]] <- e$status == r$status & e$top1 == r$top1 }
  cand <- cbind(cand, cons[, c("consensus_status", "consensus_type", "consensus_label")])
  pu_out <- pu; pu_out$measured <- pu_out$gene %in% meas_m; pu_out$bin <- ifelse(pu_out$measured, unname(bins[pu_out$gene]), NA_integer_)
  pu_out$n_controls <- ifelse(pu_out$measured, ct$table$n_controls[match(pu_out$gene, ct$table$gene)], NA_integer_)
  wc(sc, "cluster_type_scores.csv"); wc(pk, "type_top_clusters.csv"); wc(pu_out, "panel_markers_used.csv"); wc(pn$excluded, "panel_markers_excluded.csv"); wc(pn$held, "panel_markers_held_conflicts.csv")
  bm <- data.frame(bin = seq_len(GSH_N_BINS), n_genes = as.integer(bsz), min_mean = as.numeric(tapply(agg$gene_mean, bins, min)), max_mean = as.numeric(tapply(agg$gene_mean, bins, max)), stringsAsFactors = FALSE)
  wc(bm, "gene_bins.csv"); ctab <- ct$table; ctab$mean_log1p_cp10k <- unname(agg$gene_mean[ctab$gene]); wc(ctab, "control_genes.csv"); wc(vr, "check_cell_level_recompute.csv")
  wc(cand, "cluster_candidates_calc.csv")
  ug <- unique(c(meas_m, cg)); saveRDS(list(gene_mean = agg$gene_mean, bins = bins, controls = ct$controls, cluster_mean = agg$C[ug, , drop = FALSE], cluster_detection = agg$D[ug, , drop = FALSE],
                                            donor_cluster_mean = agg$K[ug, , drop = FALSE], groups = agg$groups), file.path(out, "aggregates.rds"), version = 3); files <- c(files, file.path(out, "aggregates.rds"))
  cmf <- file.path(out, "calc_manifest.json"); calc_files <- files
  writeLines(jsonlite::toJSON(list(phase = "type scores, candidates and consensus fixed before author labels are read", written_at = stamp(), inputs = as.list(h0),
                                   files = as.list(stats::setNames(vapply(files, .nrv_sha, ""), rel(files)))), auto_unbox = TRUE, pretty = TRUE), cmf)
  log("calculation fixed:", length(calc_files), "files | consensus:", paste(names(table(cons$consensus_status)), table(cons$consensus_status), collapse = " "))
  # ── (5) 사후: 저자 라벨 · 기존 성숙 후보 (Neutrophil 잠정 군집만) ──
  au <- as.character(um$author_label); comp <- gsm_author_composition(cl, au, clusters)
  atop <- vapply(clusters, function(k) comp$author_label[comp$cluster == k][1], ""); afr <- vapply(clusters, function(k) comp$frac_of_cluster[comp$cluster == k][1], 0)
  old <- rdc(f_cc); neu <- cand$consensus_status == "provisional" & cand$consensus_type == "Neutrophil"
  fin <- cand; fin$existing_maturity_top_candidate <- ifelse(neu, old$score_top_candidate[match(fin$cluster, old$cluster)], "")
  fin$maturity_status <- ifelse(neu, "unvalidated", "not_applicable"); fin$author_top_label <- unname(atop); fin$author_top_frac <- unname(afr)
  fin$note <- rep("provisional = both marker models give the same single top with score > 0; not validated; author labels post hoc only; maturity candidate shown only for Neutrophil-provisional clusters and not confirmed", nrow(fin))
  wc(fin, "cluster_provisional_annotation.csv")
  pk$peak_cluster_author_top <- vapply(pk$peak_clusters, function(z) if (!nzchar(z)) "" else paste(atop[strsplit(z, "|", fixed = TRUE)[[1]]], collapse = "|"), "")
  pk$peak_cluster_consensus <- vapply(pk$peak_clusters, function(z) if (!nzchar(z)) "" else paste(cons$consensus_label[match(strsplit(z, "|", fixed = TRUE)[[1]], cons$cluster)], collapse = "|"), "")
  wc(pk, "posthoc_type_top_clusters.csv")
  ck$add("author_labels_after_calculation", identical(unname(vapply(calc_files, .nrv_sha, "")), unname(unlist(jsonlite::fromJSON(cmf)$files[rel(calc_files)]))),
         "author labels were read only after calc_manifest.json (scores, candidates, consensus); they set no cutoff or feature; fixed calculation files unchanged afterwards")
  # ── (6) 그림 ──
  ko <- mbvr_hangul_font(); fam <- if (ko) MBVR_TITLE_FAMILY else ""; acols <- .gsu_cols(au, sl); fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
  ptypes <- unique(pu$display_label); pkeys <- intersect(ptypes, unique(cons$consensus_label[cons$consensus_status == "provisional"]))
  if (length(pkeys) > length(GSH_CAT_COLS)) stop("잠정 유형 수가 준비한 색보다 많다: ", length(pkeys))
  ccol <- c(stats::setNames(GSH_CAT_COLS[seq_along(pkeys)], pkeys), disagreement = GSH_DISAGREE_COL, undetermined = GSH_UNDET_COL)
  sub <- sprintf(if (ko) "R %s base graphics · 좌표 = %s (uwot %s, 재계산 없음) · 마커 = %s (검토 후보, 검증 전)" else "R %s base graphics · coordinates = %s (uwot %s, not recomputed) · markers = %s (review candidates, unvalidated)",
                 as.character(getRversion()), basename(dirname(a$umap)), if (is.null(um$uwot_version)) "?" else um$uwot_version, basename(a$markers))
  f1 <- file.path(fig, "fig1_umap_full_author_vs_provisional.png"); r1 <- gsh_draw_pair(f1, U, au, acols, cl, cons, ccol, sub = sub, ko = ko, fam = fam)
  f2 <- file.path(fig, "fig2_umap_lineage_author_vs_provisional.png")
  r2 <- gsh_draw_pair(f2, U, au, acols, cl, cons, ccol, mask = mask, subset_labels = sl, label_min_frac = GSH_LABEL_MIN_FRAC, sub = sub, ko = ko, fam = fam)
  files <- c(files, f1, f2); ft <- rbind(r1$text, r2$text); wc(ft, "figure_text.csv"); log("figures written")
  same_xy <- function(A, B) identical(as.vector(A), as.vector(B)) && identical(dim(A), dim(B)) && identical(dimnames(A), dimnames(B))
  ck$add("figures_same_coordinates_and_colours", same_xy(r1$coords, U) && same_xy(r2$coords, U[mask, , drop = FALSE]) && all(r1$colours$n_colours == 1L) && all(r2$colours$n_colours == 1L) &&
           identical(unname(ccol[r1$colours$key]), r1$colours$colour) && !any(ft$text %in% c("NA", "NaN", "Inf")) && !any(grepl("→|->", ft$text)),
         sprintf("figure 1 all %d cells, figure 2 all %d subset cells at the stored coordinates (values and IDs identical); one colour per cluster = its consensus; no NA text or arrows", n, sum(mask)))
  # ── (7) 기록 · 요약 · 공개 ──
  allinp <- c(markers = a$markers, umap = a$umap, lossless = f_ls, calc_manifest = f_cm, cluster_candidates = f_cc, object = a$object, baseline = a$baseline)
  h1 <- vapply(allinp, .nrv_sha, ""); bl1 <- blf()
  ck$add("existing_inputs_unchanged", identical(unname(h1), unname(h0[names(allinp)])) && identical(unname(bl1), unname(unlist(bl))),
         sprintf("%d inputs re-hashed after the run (markers, coordinates, lossless, calc_manifest, cluster_candidates, object, baseline json) and all %d baseline files (R05-74, tests) unchanged", length(allinp), length(bl)))
  wc(data.frame(input = names(allinp), path = unname(allinp), sha256 = unname(h1), stringsAsFactors = FALSE), "input_sha256.csv"); wc(snap, "source_snapshot.csv")
  pkg <- c("Matrix", "digest", "jsonlite")
  wc(data.frame(item = c("R_executable", "R", "platform", "libPaths", "direct_packages", "base_packages", "rng", "seed", "graphics", "title_font", "command", "working_dir", "started_at", "written_at"),
                value = c(file.path(R.home("bin"), "R"), R.version.string, R.version$platform, paste(.libPaths(), collapse = ";"),
                          paste(vapply(pkg, function(x) paste0(x, "=", as.character(utils::packageVersion(x))), ""), collapse = " "), "stats, utils, methods, grDevices, graphics",
                          "Mersenne-Twister / Inversion / Rejection", as.character(GSH_SEED), "R base graphics (grDevices::png, res 100)", if (ko) MBVR_TITLE_FAMILY else "default",
                          paste(c("75_gse253355_hsc_annotation.R", argv), collapse = " "), getwd(), format(t0, "%Y-%m-%dT%H:%M:%S%z"), stamp()), stringsAsFactors = FALSE), "environment.csv")
  # 요약 (20 줄 이하) · 방법
  ct2 <- table(factor(cons$consensus_status, levels = c("provisional", "disagreement", "undetermined")))
  plist <- fin[fin$consensus_status == "provisional", , drop = FALSE]; ptxt <- if (nrow(plist)) paste(sprintf("%s %s", plist$cluster, plist$consensus_label), collapse = " · ") else "없음"
  ht <- function(ty) { w <- vapply(names(GSH_MODELS), function(m) { r <- rk[[m]]; k <- r$cluster[r$status == "candidate" & r$top1 == ty]; if (length(k)) paste(k, collapse = ",") else "없음" }, "")
    p <- pk[pk$type == ty, ]; sprintf("%s — 1 위 군집 broad %s · slide2 %s; 점수 최고 군집 %s", ty, w[["broad"]], w[["slide2"]],
                                      paste(sprintf("%s %s (%.2f, 저자 %s)", p$model, p$peak_clusters, p$peak_score, p$peak_cluster_author_top), collapse = " · ")) }
  main7 <- r2$numbered; m7 <- paste(vapply(main7, function(k) { i <- match(k, fin$cluster); sprintf("%s %s (broad %s · slide2 %s)", k, if (fin$consensus_status[i] == "provisional") fin$consensus_label[i]
                                                                  else if (fin$consensus_status[i] == "undetermined") "미판정" else "불일치", fin$broad_top1[i], fin$slide2_top1[i]) }, ""), collapse = " · ")
  mat <- fin[neu, , drop = FALSE]; mtxt <- if (nrow(mat)) paste(sprintf("군집 %s 기존 최고 단계 후보 %s (unvalidated)", mat$cluster, mat$existing_maturity_top_candidate), collapse = " · ") else "Neutrophil 잠정 군집 없음"
  sm <- c("# GSE253355 — HSC 부터 전체 세포 유형 군집 잠정 후보 (요약)", "",
          sprintf("- 사용: run-r dcta R/branch_biology/75_gse253355_hsc_annotation.R --out-dir <새 폴더>. 입력 sha256: 마커 %s… · 좌표 %s… · lossless %s… · 객체 %s….",
                  substr(h0[["markers"]], 1, 12), substr(h0[["umap"]], 1, 12), substr(h0[["lossless"]], 1, 12), substr(h0[["object"]], 1, 12)),
          sprintf("- 마커: 검토 후보 CSV (approved_for_model = false) 의 사용 행만. %s. 제외 %d 행은 이유와 함께 panel_markers_excluded.csv.",
                  paste(vapply(names(GSH_MODELS), function(m) { z <- pu[pu$model == m, ]; sprintf("%s %d 유형 (양성 %d · 음성 %d)", m, length(unique(z$node_id)), sum(z$sign == 1L), sum(z$sign == -1L)) }, ""), collapse = ", "), nrow(pn$excluded)),
          sprintf("- 점수: 배경 보정 module score (log1p CP10k, 24 bin, 마커당 대조 %d 유전자, seed %d), 군집 평균 양성 − max(0, 음성). 0 이하 = 미판정. 확률 · 검증된 분류가 아니다.", GSH_N_CTRL, GSH_SEED),
          sprintf("- 결과 26 군집: 잠정 %d · 불일치 %d · 미판정 %d. 잠정: %s.", ct2[["provisional"]], ct2[["disagreement"]], ct2[["undetermined"]], ptxt),
          paste0("- HSC: ", ht("HSC"), "."), paste0("- MPP: ", ht("MPP"), "."),
          paste0("- 그림 2 주요 군집: ", m7, "."),
          paste0("- PM: 기존 6 단계 점수 · 후보는 바꾸지 않았다. ", mtxt, ". PM 문제는 여전히 해결되지 않았고 성숙 단계는 확정하지 않는다."),
          "- 한계: MSC · 내피 · 골모세포 등은 패널이 없어 이 군집들도 혈액 유형 중 가장 높은 것이 나오거나 미판정이 된다. 미판정이 이 계통을 배제한다는 뜻이 아니다.",
          "- 한계: 마커 목록 (검토 후보) · 군집 해상도 · 대조 표본에 의존한다. bootstrap 없음. 저자 라벨은 계산 고정 뒤 사후 대조에만 썼다.",
          "- 표: cluster_provisional_annotation · cluster_type_scores · type_top_clusters · panel_markers_used / excluded · control_genes · gene_bins. 그림: figures/fig1 · fig2.")
  mt <- c("# 방법 — 군집 잠정 세포 유형 후보 (R75)", "",
          "1. 마커: node_markers_sign_review_2026-10-18.csv 원본 사본 (marker_input/). use_in_current_panel = true, is_annotation_candidate = true, marker_sign = ±1 행만 쓴다.",
          "2. 모델: bm_dag_broad 와 bm_dag_slide2 를 따로 쓰고, 각각에 t_cell_broad_categories 를 더한다. node · gene · sign 중복은 지우고, 같은 node 의 양 · 음 충돌은 보류한다. DAG edge 는 쓰지 않는다.",
          "3. 발현: log1p(1e4 · count / 모든 유전자 library). 유전자 전체 세포 평균으로 24 개 같은 개수 bin 을 만든다(동률은 symbol 순서).",
          sprintf("4. 대조: 마커 유전자마다 같은 bin 의 비마커 유전자 %d 개를 뽑는다(seed %d, Mersenne-Twister / Rejection, 마커 symbol 순서, 모자라면 전부). 목록은 control_genes.csv.", GSH_N_CTRL, GSH_SEED),
          "5. 유전자 보정값 = log 발현 − 대조 평균. 양성 · 음성 패널은 측정 유전자 동일 가중 평균이고, 없는 마커는 missing_markers 에 적는다.",
          "6. 군집 점수 = 군집 평균 양성 − max(0, 군집 평균 음성). 음성 패널이 없으면 0 을 빼고, 양성 측정 유전자가 없으면 NA.",
          "   선형인 부분은 양성 · 음성 패널 평균의 군집 집계까지다: 그래서 유전자별 군집 평균으로 계산했고 세포 단위 재계산과 같다. 음성 벌점 max(0, ·) 는 비선형이며 세포마다가 아니라 군집 평균 음성에 적용한다.",
          "7. donor 동일 가중: donor x 군집 >= 5 세포 donor 의 양성 · 음성 평균을 같은 무게로 평균한 뒤 같은 식을 쓴다. 해당 donor 가 없으면 미판정이다.",
          "8. 후보: 1 위가 0 이하이면 미판정이다. 2 위 = 다음 유형과 미판정 (0) 중 큰 것이다. 동률은 이름 하나로 깨지 않는다.",
          "9. 합의: 두 모델 1 위가 같은 단일 유형이고 점수 > 0 → 잠정 후보. 둘 다 미판정 → 미판정. 그 밖 → 불일치.",
          "10. 성숙 단계: 기존 6 단계 점수는 그대로다. Neutrophil 잠정 군집에만 기존 최고 단계 후보를 보조로 적고 unvalidated 로 둔다.")
  pub_list <- c("summary_ko.md", "methods_ko.md", "cluster_provisional_annotation.csv", "cluster_type_scores.csv", "posthoc_type_top_clusters.csv", "panel_markers_used.csv",
                "panel_markers_excluded.csv", "panel_markers_held_conflicts.csv", "control_genes.csv", "gene_bins.csv", "check_cell_level_recompute.csv", "checks.csv", "figure_text.csv",
                file.path("marker_input", c(basename(a$markers), "source.csv")), file.path("figures", basename(c(f1, f2))))
  ck$add("tables_roundtrip_ok", all(crt), sprintf("%d CSVs re-read with identical dimensions, names and cell text", length(crt)))
  ck$add("public_copy_aggregate_only", !any(grepl("lossless|coordinates|aggregates|\\.rds$|\\.gz$", pub_list)), sprintf("%d public files: aggregate tables, PNG, input marker copy, methods, checks; no cell-level data", length(pub_list)))
  cks <- ck$table(); sm <- c(sm, sprintf("- 검사 %d / %d 통과 (checks.csv).", sum(cks$pass), nrow(cks)))
  writeLines(sm, file.path(out, "summary_ko.md")); writeLines(mt, file.path(out, "methods_ko.md")); files <- c(files, file.path(out, c("summary_ko.md", "methods_ko.md"))); wc(cks, "checks.csv")
  allf <- c(files, cmf)
  writeLines(jsonlite::toJSON(list(output = basename(out), written_at = stamp(), inputs = as.list(h0), files = as.list(stats::setNames(vapply(allf, .nrv_sha, ""), rel(allf)))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(out, "manifest.json"))
  if (!all(cks$pass)) stop("확인 실패: ", paste(cks$check[!cks$pass], collapse = ", "))
  pub <- a$public; for (dd in unique(dirname(file.path(pub, pub_list)))) dir.create(dd, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(file.path(out, pub_list), file.path(pub, pub_list), overwrite = FALSE))) stop("공개 사본 복사 실패")
  writeLines(jsonlite::toJSON(list(source_output = normalizePath(out), inputs = as.list(h0), files = as.list(stats::setNames(vapply(file.path(pub, pub_list), .nrv_sha, ""), pub_list))),
                              auto_unbox = TRUE, pretty = TRUE), file.path(pub, "public_manifest.json"))
  log("done | checks:", paste(cks$check, cks$pass, collapse = ", "))
}


if (!interactive() && sys.nframe() == 0L) {
  .gsh_bb <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
  for (.f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")) source(file.path(.gsh_bb, .f))
  quit(status = gsh_main(commandArgs(trailingOnly = TRUE)))
}
