# =============================================================================
# 63_atlas_rna_reference.R — atlas RNA (GSE294330) source 전용 인접 단계 대비 · 단계 양성 패널 (target 을 읽지 않는다)
# =============================================================================
#
# 계약: results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_instruction.md, docs/experiments/neutrophil_atlas_rna_2026-10-19/design.md, plan Task 1.
#   - source: GSE294330 mRNA TPM 40 열 + GEO full SOFT (열 ↔ GSM · 단계 · code · 세부 분획 · 병합 표기) + HGNC (정확한 ENSG 일대일만). counts 는 행 · 열 대조에만 쓴다.
#   - code x 단계마다 TPM 중앙값 → log2(TPM + 1). code(A–E)는 nominal collection group 이다. 논문 4 인 + 반복 1 과의 동일인 관계는 확인되지 않았으므로
#     독립 공여자 검증 · p 값 · 정확도 · bootstrap CI 를 내지 않는다. 같은 code 의 재시퀀싱 · 세부 분획 열은 그 code x 단계 하나로 묶는다(독립 표본 아님).
#   - 인접 대비 5 개 (뒤 − 앞): 두 단계가 모두 있는 code 만. 후보 = 쌍 ≥ 3 · 높은 쪽 TPM ≥ 1 인 code ≥ 3 · |중앙값| ≥ 0.5 · 중앙값 부호 일치율 ≥ 0.75.
#     부호별 |중앙값| 순 top 20 (동률 ENSG 순). 모든 유전자 x 모든 대비 표도 보존한다. TPM 차이는 서술적 효과이며 DE 검정이 아니다.
#   - 단계 양성 패널: 있는 모든 인접 단계보다 높음(각 대비 같은 후보 기준)의 교집합. 끝 단계는 이웃 하나. 최소 인접 효과 순 top 10 (동률 ENSG).
#     후보가 없으면 빈 패널(판정불가)로 두고 억지로 채우지 않는다. DCTA 실험 규칙이며 논문 저자의 공식 패널이 아니다.
#   - 선택 목록(변화 · 패널)은 HGNC 일대일 protein-coding 만 쓴다. 그 밖의 후보는 제외 표에 이유와 함께 남긴다. 혈액 PMN 은 골수 6 단계와 분리한다.
#   - 파라미터 · 순위는 target RNA · 라벨 · 연결을 보기 전에 고정한다.
# =============================================================================

ATLAS_STAGES <- c("Myeloblast", "PM", "MY", "MM", "BN", "SN")                                 # 골수 6 단계 (기존 순서)
ATLAS_CELLTYPE <- c(Myeloblast = "Myeloblast", Promyelocyte = "PM", Myelocyte = "MY", Metamyelocyte = "MM", Band = "BN", Segmented = "SN", Neutrophil = "PB_PMN")
ATLAS_PARAMS <- list(min_codes = 3L, min_tpm = 1, min_codes_tpm = 3L, min_abs_median = 0.5, min_sign_agree = 0.75, top_transition = 20L, top_panel = 10L)   # design 고정값
ATLAS_FILES <- c(tpm = "GSE294330_exp_mrn_tpm.csv.gz", counts = "GSE294330_exp_mrn_counts.csv.gz", soft = "GSE294330_full.soft", hgnc = "hgnc_complete_set.txt")
ATLAS_INDEPENDENCE_NOTE <- paste("Paper: 4 donors + 1 biological repeat; GEO: 5 nominal collection codes (A-E). Same-person mapping not confirmed:",
                                 "codes are nominal collection groups, not independent donors; no donor-level CV, p-value, accuracy or bootstrap CI.")


# ── 1. source 읽기 ───────────────────────────────────────────────────────────
#' SOFT 의 표본 블록에서 GSM · title · cell type · 행렬 열 이름을 읽는다 (열 이름이 있는 표본만)
atlas_parse_soft <- function(path) {
  x <- readLines(path, warn = FALSE); st <- which(startsWith(x, "^SAMPLE = "))
  if (!length(st)) stop("SOFT 에 SAMPLE 블록이 없다: ", path)
  en <- c(st[-1] - 1L, length(x))
  d <- do.call(rbind, lapply(seq_along(st), function(k) { b <- x[st[k]:en[k]]
    val <- function(prefix) substring(b[startsWith(b, prefix)], nchar(prefix) + 1L)
    de <- val("!Sample_description = "); col <- de[startsWith(de, "Column name in raw_counts_mrn.csv: ")]; col <- substring(col, nchar("Column name in raw_counts_mrn.csv: ") + 1L)
    ch <- val("!Sample_characteristics_ch1 = "); ct <- ch[startsWith(ch, "cell type: ")]; ct <- substring(ct, nchar("cell type: ") + 1L)
    data.frame(gsm = substring(x[st[k]], nchar("^SAMPLE = ") + 1L), title = c(val("!Sample_title = "), "")[1], cell_type = c(ct, "")[1], column = c(col, "")[1],
               stringsAsFactors = FALSE) }))
  d[nzchar(d$column), , drop = FALSE]
}

#' HGNC 대응: source ENSG 마다 HGNC 행 수(정확한 ensembl_gene_id). 1 행 = one_to_one, 2 행 이상 = multi_hgnc, 0 행 = unmapped,
#' 일대일이어도 같은 symbol 이 source ENSG 여럿에 붙으면 shared_symbol. 패널 사용 = one_to_one 이고 protein-coding gene.
atlas_map_hgnc <- function(ensg, hg) {
  miss <- setdiff(c("symbol", "ensembl_gene_id", "locus_group"), names(hg)); if (length(miss)) stop("HGNC 표에 필요한 열이 없다: ", paste(miss, collapse = ", "))
  h <- hg[nzchar(hg$ensembl_gene_id), , drop = FALSE]; nr <- as.integer(table(factor(h$ensembl_gene_id, levels = unique(ensg)))[ensg]); j <- match(ensg, h$ensembl_gene_id)
  st <- ifelse(nr == 0L, "unmapped", ifelse(nr > 1L, "multi_hgnc", "one_to_one")); sym <- ifelse(st == "one_to_one", h$symbol[j], "")
  dup <- st == "one_to_one" & (duplicated(sym) | duplicated(sym, fromLast = TRUE)); st[dup] <- "shared_symbol"; sym[dup] <- ""
  lg <- ifelse(st == "one_to_one", h$locus_group[j], "")
  data.frame(ensg = ensg, symbol = sym, hgnc_id = if ("hgnc_id" %in% names(h)) ifelse(st == "one_to_one", h$hgnc_id[j], "") else "", n_hgnc_rows = nr,
             locus_group = lg, map_status = st, use_for_panels = st == "one_to_one" & lg == "protein-coding gene", row.names = NULL, stringsAsFactors = FALSE)
}

#' source 읽기: TPM (유전자 x 40 열), 열 metadata (SOFT 로 대응), HGNC 대응, 파일 hash. 대응이 어긋나면 멈춘다.
atlas_read_source <- function(data_dir) {
  p <- stats::setNames(file.path(data_dir, ATLAS_FILES), names(ATLAS_FILES))
  miss <- names(p)[!file.exists(p)]; if (length(miss)) stop("source 파일이 없다: ", paste(basename(p[miss]), collapse = ", "))
  rd <- function(f) utils::read.csv(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  tp <- rd(p[["tpm"]]); ct <- rd(p[["counts"]])
  if (!"ENSG" %in% names(tp)) stop("TPM 표에 ENSG 열이 없다.")
  if (!identical(names(ct), names(tp)) || !identical(ct$ENSG, tp$ENSG)) stop("counts 와 TPM 의 행 · 열이 다르다.")
  if (anyDuplicated(tp$ENSG)) stop("TPM 표의 ENSG 가 중복이다.")
  cols <- setdiff(names(tp), "ENSG"); if (anyDuplicated(cols)) stop("TPM 열 이름이 중복이다.")
  tpm <- as.matrix(tp[, cols, drop = FALSE]); rownames(tpm) <- tp$ENSG
  if (!is.numeric(tpm) || anyNA(tpm) || any(tpm < 0)) stop("TPM 값이 음수 · 결측 · 비숫자다.")
  sf <- atlas_parse_soft(p[["soft"]]); sf <- sf[startsWith(sf$title, "mRNA."), , drop = FALSE]
  if (anyDuplicated(sf$column)) stop("SOFT 에서 같은 열 이름이 두 mRNA 표본에 있다.")
  miss <- setdiff(cols, sf$column); if (length(miss)) stop("행렬 열이 SOFT mRNA 표본에 없다: ", paste(miss, collapse = ", "))
  extra <- setdiff(sf$column, cols); if (length(extra)) stop("SOFT mRNA 열이 행렬에 없다: ", paste(extra, collapse = ", "))
  sf <- sf[match(cols, sf$column), , drop = FALSE]
  stage <- unname(ATLAS_CELLTYPE[sf$cell_type]); if (anyNA(stage)) stop("모르는 cell type: ", paste(unique(sf$cell_type[is.na(stage)]), collapse = ", "))
  code <- sub("^.*_donor_([A-Za-z])(_.*)?$", "\\1", sf$title); bad <- !grepl("^[A-Z]$", code)
  if (any(bad)) stop("title 에서 donor code 를 읽지 못했다: ", paste(sf$title[bad], collapse = ", "))
  md <- data.frame(column = cols, gsm = sf$gsm, title = sf$title, cell_type = sf$cell_type, stage = stage, code = code,
                   subfraction = ifelse(grepl("\\.early_donor", sf$title), "early", ifelse(grepl("\\.late_donor", sf$title), "late", "")),
                   merge_note = ifelse(grepl("_mergedInFinalWith_", sf$title, fixed = TRUE), sub("^.*_mergedInFinalWith_", "", sf$title), ""),
                   in_bm = stage %in% ATLAS_STAGES, row.names = NULL, stringsAsFactors = FALSE)
  hg <- utils::read.delim(p[["hgnc"]], stringsAsFactors = FALSE, quote = "", colClasses = "character")
  list(tpm = tpm, metadata = md, map = atlas_map_hgnc(rownames(tpm), hg),
       files = data.frame(role = names(p), file = basename(p), path = unname(p), sha256 = vapply(p, function(f) digest::digest(f, algo = "sha256", file = TRUE), "", USE.NAMES = FALSE),
                          stringsAsFactors = FALSE))
}


# ── 2. code x 단계 요약 ──────────────────────────────────────────────────────
#' code x 단계마다 TPM 중앙값 (같은 code x 단계의 재시퀀싱 · 세부 분획 열을 묶는다). 반환: tpm_med (유전자 x group), groups (group · code · stage · 열)
atlas_code_stage_tpm <- function(tpm, md) {
  if (!identical(colnames(tpm), md$column)) stop("TPM 열과 metadata 열 순서가 다르다.")
  key <- paste(md$code, md$stage, sep = "|"); fk <- !duplicated(key)
  uk <- key[fk][order(match(md$stage[fk], c(ATLAS_STAGES, "PB_PMN")), md$code[fk], method = "radix")]
  med <- matrix(NA_real_, nrow(tpm), length(uk), dimnames = list(rownames(tpm), uk))
  for (k in uk) { j <- which(key == k); med[, k] <- if (length(j) == 1L) tpm[, j] else apply(tpm[, j, drop = FALSE], 1, stats::median) }
  grp <- data.frame(group = uk, code = sub("\\|.*$", "", uk), stage = sub("^.*\\|", "", uk), n_columns = as.integer(table(factor(key, levels = uk))),
                    columns = vapply(uk, function(k) paste(md$column[key == k], collapse = "|"), ""),
                    merge_notes = vapply(uk, function(k) paste(unique(md$merge_note[key == k & nzchar(md$merge_note)]), collapse = "|"), ""),
                    subfractions = vapply(uk, function(k) paste(unique(md$subfraction[key == k & nzchar(md$subfraction)]), collapse = "|"), ""),
                    row.names = NULL, stringsAsFactors = FALSE)
  list(tpm_med = med, groups = grp)
}


# ── 3. 인접 대비 · 선택 변화 · 단계 패널 ─────────────────────────────────────
.atlas_exclusion_reason <- function(map_status) ifelse(map_status == "unmapped", "map_unmapped", ifelse(map_status == "multi_hgnc", "map_multi_hgnc",
                                                       ifelse(map_status == "shared_symbol", "map_shared_symbol", "not_protein_coding")))

#' source 전용 reference. 반환: params · stages · groups · all_contrasts · transitions · panels · panel_status · exclusions · independence_note
atlas_build_reference <- function(source, params = ATLAS_PARAMS, stages = ATLAS_STAGES) {
  md <- source$metadata; keep <- md$stage %in% stages; md <- md[keep, , drop = FALSE]; tpm <- source$tpm[, keep, drop = FALSE]
  g <- atlas_code_stage_tpm(tpm, md); L <- log2(g$tpm_med + 1); codes <- sort(unique(g$groups$code))
  mp <- source$map[match(rownames(tpm), source$map$ensg), , drop = FALSE]; if (anyNA(mp$ensg)) stop("HGNC 대응 표에 없는 source 유전자가 있다.")
  pr <- data.frame(from = stages[-length(stages)], to = stages[-1], stringsAsFactors = FALSE); pr$pair <- paste(pr$from, pr$to, sep = ">")
  allc <- vector("list", nrow(pr))
  for (r in seq_len(nrow(pr))) { a <- pr$from[r]; b <- pr$to[r]
    cc <- sort(intersect(g$groups$code[g$groups$stage == a], g$groups$code[g$groups$stage == b])); n <- length(cc)
    ga <- paste(cc, a, sep = "|"); gb <- paste(cc, b, sep = "|"); D <- L[, gb, drop = FALSE] - L[, ga, drop = FALSE]
    med <- if (n) apply(D, 1, stats::median) else rep(NA_real_, nrow(L))
    agree <- if (n) rowMeans(sign(D) == sign(med)) else rep(NA_real_, nrow(L))
    hi <- if (n) ifelse(matrix(med > 0, nrow(L), n), g$tpm_med[, gb, drop = FALSE], g$tpm_med[, ga, drop = FALSE]) else matrix(0, nrow(L), 0)
    nhi <- as.integer(rowSums(hi >= params$min_tpm))
    cand <- n >= params$min_codes & !is.na(med) & med != 0 & abs(med) >= params$min_abs_median & !is.na(agree) & agree >= params$min_sign_agree & nhi >= params$min_codes_tpm
    d <- data.frame(pair = pr$pair[r], from = a, to = b, ensg = rownames(L), symbol = mp$symbol, map_status = mp$map_status, use_for_panels = mp$use_for_panels, n_pairs = n,
                    codes = paste(cc, collapse = "|"), median_diff = unname(med), sign_agree = unname(agree), n_codes_high_tpm = nhi, candidate = unname(cand),
                    direction = ifelse(cand, ifelse(med > 0, "up", "down"), ""), stringsAsFactors = FALSE)
    for (cd in codes) { ok <- cd %in% cc
      d[[paste0("tpm_from_", cd)]] <- if (ok) unname(g$tpm_med[, paste(cd, a, sep = "|")]) else NA_real_
      d[[paste0("tpm_to_", cd)]] <- if (ok) unname(g$tpm_med[, paste(cd, b, sep = "|")]) else NA_real_
      d[[paste0("diff_", cd)]] <- if (ok) unname(D[, match(cd, cc)]) else NA_real_ }
    allc[[r]] <- d }
  allc <- do.call(rbind, allc); rownames(allc) <- NULL
  # 선택 변화: 부호별 |중앙값| 순 top N (동률 ENSG), 일대일 protein-coding 만
  tr <- list(); ex <- list()
  for (p in pr$pair) for (dr in c("up", "down")) { x <- allc[allc$pair == p & allc$candidate & allc$direction == dr, , drop = FALSE]
    if (any(!x$use_for_panels)) { y <- x[!x$use_for_panels, ]
      ex[[length(ex) + 1L]] <- data.frame(kind = "transition", stage_or_pair = p, direction = dr, ensg = y$ensg, symbol = y$symbol, reason = .atlas_exclusion_reason(y$map_status),
                                          stringsAsFactors = FALSE) }
    x <- x[x$use_for_panels, , drop = FALSE]; x <- x[order(-abs(x$median_diff), x$ensg, method = "radix"), , drop = FALSE]; x <- utils::head(x, params$top_transition)
    if (nrow(x)) tr[[length(tr) + 1L]] <- data.frame(pair = p, direction = dr, rank = seq_len(nrow(x)), ensg = x$ensg, symbol = x$symbol, median_diff = x$median_diff,
                                                     n_pairs = x$n_pairs, sign_agree = x$sign_agree, n_codes_high_tpm = x$n_codes_high_tpm, stringsAsFactors = FALSE) }
  # 단계 양성 패널: 있는 모든 인접 단계보다 높음의 교집합, 최소 인접 효과 순 top N
  pn <- list(); ps <- list()
  for (k in seq_along(stages)) { s <- stages[k]; nb <- character(0); conds <- list()
    if (k > 1L) { nb <- c(nb, stages[k - 1L]); conds[["prev"]] <- allc[allc$pair == paste(stages[k - 1L], s, sep = ">"), , drop = FALSE] }
    if (k < length(stages)) { nb <- c(nb, stages[k + 1L]); conds[["next"]] <- allc[allc$pair == paste(s, stages[k + 1L], sep = ">"), , drop = FALSE] }
    ok <- rep(TRUE, nrow(L)); eff <- matrix(NA_real_, nrow(L), 2, dimnames = list(NULL, c("prev", "next")))
    for (cn in names(conds)) { z <- conds[[cn]]; z <- z[match(rownames(L), z$ensg), ]
      ok <- ok & z$candidate & z$direction == if (cn == "prev") "up" else "down"; eff[, cn] <- z$median_diff }
    me <- apply(abs(eff[, names(conds), drop = FALSE]), 1, min)
    x <- data.frame(ensg = rownames(L), symbol = mp$symbol, map_status = mp$map_status, use = mp$use_for_panels, min_effect = me, effect_vs_prev = eff[, "prev"], effect_vs_next = eff[, "next"],
                    stringsAsFactors = FALSE)[ok, , drop = FALSE]
    if (any(!x$use)) { y <- x[!x$use, ]
      ex[[length(ex) + 1L]] <- data.frame(kind = "panel", stage_or_pair = s, direction = "positive", ensg = y$ensg, symbol = y$symbol, reason = .atlas_exclusion_reason(y$map_status),
                                          stringsAsFactors = FALSE) }
    y <- x[x$use, , drop = FALSE]; y <- y[order(-y$min_effect, y$ensg, method = "radix"), , drop = FALSE]; y <- utils::head(y, params$top_panel)
    if (nrow(y)) pn[[length(pn) + 1L]] <- data.frame(stage = s, rank = seq_len(nrow(y)), ensg = y$ensg, symbol = y$symbol, min_adjacent_effect = y$min_effect,
                                                     effect_vs_prev = y$effect_vs_prev, effect_vs_next = y$effect_vs_next, stringsAsFactors = FALSE)
    ps[[k]] <- data.frame(stage = s, neighbors = paste(nb, collapse = "|"), n_candidates = nrow(x), n_excluded_mapping = sum(!x$use), n_genes = nrow(y),
                          status = if (nrow(y)) "panel_listed" else if (nrow(x)) "empty_after_mapping_filter" else "empty_no_candidate", stringsAsFactors = FALSE) }
  bind <- function(z, tmpl) { if (!length(z)) return(tmpl); out <- do.call(rbind, z); rownames(out) <- NULL; out }
  list(params = params, stages = stages, codes = codes, groups = g$groups, all_contrasts = allc,
       transitions = bind(tr, data.frame(pair = character(0), direction = character(0), rank = integer(0), ensg = character(0), symbol = character(0), median_diff = numeric(0),
                                         n_pairs = integer(0), sign_agree = numeric(0), n_codes_high_tpm = integer(0), stringsAsFactors = FALSE)),
       panels = bind(pn, data.frame(stage = character(0), rank = integer(0), ensg = character(0), symbol = character(0), min_adjacent_effect = numeric(0), effect_vs_prev = numeric(0),
                                    effect_vs_next = numeric(0), stringsAsFactors = FALSE)),
       panel_status = bind(ps, NULL),
       exclusions = bind(ex, data.frame(kind = character(0), stage_or_pair = character(0), direction = character(0), ensg = character(0), symbol = character(0), reason = character(0),
                                        stringsAsFactors = FALSE)),
       independence_note = ATLAS_INDEPENDENCE_NOTE)
}
