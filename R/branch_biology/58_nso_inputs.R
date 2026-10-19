# =============================================================================
# 58_nso_inputs.R — 호중구 성숙 단계 관측(NSO): reference · 입력 검사 · 계산용 메타데이터 분리
# =============================================================================
#
# 계약: results/claude_bridge/neut_stage_observation_2026-10-19/codex_instruction.md (2026-10-19), 계획 Task 1. Ruling: 같은 폴더 claude_progress.md
#   - 형태 단계 6 개(Myeloblast → PM → MY → MM → BN → SN)의 순서를 그대로 고정한다. SN → PB_PMN 은 골수 → 혈액 이동이라 따로 둔다.
#   - 참고표 marker 는 문헌 기반 RNA 프로그램 요약이다. 단계 정답 · 배타 marker 가 아니다. 여러 단계가 같이 쓰는 marker 는 지우지 않고 표시한다.
#   - 단계 사이 변화(E05–E09)는 core / conditional, up / down 을 나눠 남긴다. 빈칸 = 근거 미지정 (변화 없음이 아니다).
#   - 계산 입력에는 원저자 라벨과 UMAP 을 넣지 않는다. 들어오면 조용히 버리지 않고 멈춘다.
#   R05 · R17(ID · 열 검사), R23(참고표 분해), R24(donor block) 를 수정 없이 재사용한다.
#
# 주요 함수
#   nso_reference(nodes_path, edges_path)                              6 단계 정의 · marker 긴 표 · 변화 긴 표 · 공유 marker · PB 표 · scope
#   nso_check_inputs(cells92, cells51, raw_ids, graph_ids, reference)  네 ID 목록 · 군집 · block 검사 → canonical cell_id
#   nso_calculation_metadata(cells)                                    cell_id · cluster92 · cluster51 · donor · study · block 만
# =============================================================================

NSO_STAGES <- c("Myeloblast", "PM", "MY", "MM", "BN", "SN")                          # 기존 형태 순서 (바꾸지 않는다)
NSO_PB_NODE <- "PB_PMN"
NSO_NODE_COLS <- c("node_id", "stage_order", "display_name", "full_name", "compartment", "marker_genes", "program_ko", "state_note_ko", "reference_doi",
                   "adult_evidence_scope")
NSO_EDGE_COLS <- c("edge_id", "from", "to", "edge_type", "up_genes", "down_genes", "conditional_up_genes", "conditional_down_genes", "note_ko", "direction_basis",
                   "adult_evidence_scope", "reference_doi")
NSO_CHANGE_COLS <- list(c("up_genes", "up", "core"), c("down_genes", "down", "core"), c("conditional_up_genes", "up", "conditional"),
                        c("conditional_down_genes", "down", "conditional"))
NSO_META_COLS <- c("cell_id", "cluster92", "cluster51", "donor", "study", "block")


# ── 1. reference ─────────────────────────────────────────────────────────────
#' 참고표 하나를 문자 그대로 읽는다 (R23 과 같은 CSV 옵션, "NA" 를 결측으로 바꾸지 않는다). 0 행 · 필수 열 누락은 멈춘다.
.nso_read_table <- function(path, need, what) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE, colClasses = "character", na.strings = character(0))
  miss <- setdiff(need, names(d)); if (length(miss)) stop("참고표 ", what, " 표에 필수 열이 없다: ", paste(miss, collapse = ", "), " (", path, ")")
  if (!nrow(d)) stop("참고표 ", what, " 표가 0 행이다: ", path)
  id <- if (what == "노드") d$node_id else d$edge_id
  if (anyDuplicated(id)) stop("참고표 ", what, " ID 가 중복이다: ", paste(unique(id[duplicated(id)]), collapse = ", "))
  d
}

#' 형태 6 단계 정의 · marker 긴 표 · 단계 사이 변화 긴 표 · 단계 사이 연결 표 · 공유 marker · PB 표 · scope · 파일 hash
nso_reference <- function(nodes_path, edges_path) {
  nodes <- .nso_read_table(nodes_path, NSO_NODE_COLS, "노드"); edges <- .nso_read_table(edges_path, NSO_EDGE_COLS, "edge")
  miss <- setdiff(NSO_STAGES, nodes$node_id); if (length(miss)) stop("형태 단계가 참고표에 없다: ", paste(miss, collapse = ", "))
  st <- nodes[match(NSO_STAGES, nodes$node_id), , drop = FALSE]; so <- suppressWarnings(as.integer(st$stage_order))
  if (anyNA(so) || any(diff(so) <= 0)) stop("형태 단계 순서가 Myeloblast → PM → MY → MM → BN → SN 이 아니다 (stage_order: ",
                                            paste(st$node_id, st$stage_order, sep = "=", collapse = ", "), ").")
  gl <- lapply(st$marker_genes, .nrv_split); names(gl) <- NSO_STAGES
  for (s in NSO_STAGES) if (anyDuplicated(gl[[s]])) stop("한 단계의 marker 목록에 같은 유전자가 두 번 있다: ", s)
  pbm <- if (NSO_PB_NODE %in% nodes$node_id) .nrv_split(nodes$marker_genes[nodes$node_id == NSO_PB_NODE]) else character(0)
  stage_of <- rep(NSO_STAGES, lengths(gl)); gene <- unlist(gl, use.names = FALSE)
  mk <- data.frame(stage_id = stage_of, stage_rank = match(stage_of, NSO_STAGES), gene = gene, gene_order = unlist(lapply(gl, seq_along), use.names = FALSE),
                   shared_with = vapply(seq_along(gene), function(i) paste(setdiff(stage_of[gene == gene[i]], stage_of[i]), collapse = "|"), ""),
                   pb_shared = gene %in% pbm, marker_kind = "rna_program_summary", evidence_scope = st$adult_evidence_scope[match(stage_of, NSO_STAGES)],
                   reference_doi = st$reference_doi[match(stage_of, NSO_STAGES)], stringsAsFactors = FALSE)
  stages <- data.frame(stage_id = NSO_STAGES, stage_rank = seq_along(NSO_STAGES), source_stage_order = so, display_name = st$display_name, full_name = st$full_name,
                       compartment = st$compartment, program_ko = st$program_ko, state_note_ko = st$state_note_ko, reference_doi = st$reference_doi,
                       adult_evidence_scope = st$adult_evidence_scope, n_markers = as.integer(lengths(gl)), marker_genes = vapply(gl, paste, "", collapse = "|"),
                       marker_status = ifelse(lengths(gl) > 0, "markers_listed", "no_reference_markers"), row.names = NULL, stringsAsFactors = FALSE)
  ug <- unique(c(gene, pbm[pbm %in% gene]))
  sh <- do.call(rbind, lapply(ug, function(g) { s <- stage_of[gene == g]; r <- match(s, NSO_STAGES)
    data.frame(gene = g, stages = paste(s, collapse = "|"), n_stages = length(s), adjacent_only = all(diff(r) == 1L), pb_shared = g %in% pbm, first_rank = min(r),
               stringsAsFactors = FALSE) }))
  sh <- sh[sh$n_stages >= 2L | sh$pb_shared, , drop = FALSE]; sh <- sh[order(sh$first_rank, match(sh$gene, gene), method = "radix"), setdiff(names(sh), "first_rank"), drop = FALSE]
  rownames(sh) <- NULL
  # 단계 사이 변화: 양끝이 모두 형태 단계인 edge 는 바로 다음 단계로 가야 한다 (건너뜀 · 역방향은 오류)
  inS <- function(x) x %in% NSO_STAGES; tr_e <- edges[inS(edges$from) & inS(edges$to), , drop = FALSE]
  bad <- tr_e$edge_id[match(tr_e$to, NSO_STAGES) != match(tr_e$from, NSO_STAGES) + 1L]
  if (length(bad)) stop("단계 사이 edge 가 바로 다음 단계로 가지 않는다: ", paste(bad, "(", tr_e$from[tr_e$edge_id %in% bad], "→", tr_e$to[tr_e$edge_id %in% bad], ")", collapse = ", "))
  long <- function(E) do.call(rbind, lapply(seq_len(nrow(E)), function(i) do.call(rbind, lapply(NSO_CHANGE_COLS, function(cc) {
    g <- .nrv_split(E[[cc[1]]][i]); if (!length(g)) return(NULL)
    data.frame(edge_id = E$edge_id[i], from = E$from[i], to = E$to[i], gene = g, direction = cc[2], evidence = cc[3], gene_order = seq_along(g),
               direction_basis = E$direction_basis[i], reference_support = nrv_reference_support(E$direction_basis[i]),
               reference_scope = if (is.na(NRV_EDGE_LAYER[E$edge_id[i]])) "unknown_layer" else unname(NRV_EDGE_LAYER[E$edge_id[i]]),
               adult_evidence_scope = E$adult_evidence_scope[i], reference_doi = E$reference_doi[i], note_ko = E$note_ko[i], stringsAsFactors = FALSE) }))))
  ch <- long(tr_e); if (is.null(ch)) ch <- data.frame(edge_id = character(0), from = character(0), to = character(0), gene = character(0), direction = character(0),
                                                     evidence = character(0), gene_order = integer(0), stringsAsFactors = FALSE)
  rownames(ch) <- NULL
  tr <- do.call(rbind, lapply(seq_len(length(NSO_STAGES) - 1L), function(k) { a <- NSO_STAGES[k]; b <- NSO_STAGES[k + 1L]; id <- tr_e$edge_id[tr_e$from == a & tr_e$to == b]
    n <- function(dr, ev) as.integer(sum(ch$from == a & ch$direction == dr & ch$evidence == ev))
    data.frame(from = a, to = b, edge_id = if (length(id)) id else "", n_core_up = n("up", "core"), n_core_down = n("down", "core"), n_conditional_up = n("up", "conditional"),
               n_conditional_down = n("down", "conditional"), status = if (!length(id)) "no_reference_edge" else if (any(ch$from == a)) "change_rules_listed" else "edge_without_genes",
               stringsAsFactors = FALSE) }))
  # PB (골수 → 혈액) 는 따로: PB_PMN node · marker, marrow_to_blood edge 변화
  pb_e <- edges[edges$to == NSO_PB_NODE | edges$edge_type == "marrow_to_blood", , drop = FALSE]
  pb <- list(node = nodes[nodes$node_id == NSO_PB_NODE, , drop = FALSE], markers = data.frame(node_id = rep(NSO_PB_NODE, length(pbm)), gene = pbm, stringsAsFactors = FALSE),
             changes = if (nrow(pb_e)) long(pb_e) else NULL)
  # scope: 모든 node · edge 를 행으로 남기고 역할을 붙인다
  nrole <- ifelse(inS(nodes$node_id), "morphology_stage", ifelse(nodes$node_id == NSO_PB_NODE | nodes$compartment == "peripheral_blood", "marrow_to_blood_target", "upstream_context"))
  erole <- ifelse(inS(edges$from) & inS(edges$to), "stage_transition", ifelse(edges$to == NSO_PB_NODE | edges$edge_type == "marrow_to_blood", "marrow_to_blood",
                  ifelse(inS(edges$to), "entry_projection", ifelse(inS(edges$from), "exit_projection", "upstream_projection"))))
  scope <- rbind(data.frame(kind = "node", id = nodes$node_id, from = "", to = "", role = nrole, used_for_stage_evidence = nrole == "morphology_stage", stringsAsFactors = FALSE),
                 data.frame(kind = "edge", id = edges$edge_id, from = edges$from, to = edges$to, role = erole, used_for_stage_evidence = erole == "stage_transition",
                            stringsAsFactors = FALSE))
  files <- data.frame(file = basename(c(nodes_path, edges_path)), path = c(nodes_path, edges_path),
                      sha256 = vapply(c(nodes_path, edges_path), function(p) digest::digest(p, algo = "sha256", file = TRUE), "", USE.NAMES = FALSE), stringsAsFactors = FALSE)
  list(stages = stages, stage_markers = mk, changes = ch, transitions = tr, shared_markers = sh, pb = pb, scope = scope, files = files)
}


# ── 2. 입력 검사 ─────────────────────────────────────────────────────────────
#' UMAP · 표시 좌표 열이 계산 입력에 있으면 멈춘다 (저자 · DAG 열은 R17 mdu_check_fit_columns 가 막는다)
.nso_check_display_columns <- function(cols) {
  bad <- cols[grepl("umap|tsne", cols, ignore.case = TRUE)]
  if (length(bad)) stop("UMAP · 표시 좌표 열은 계산 입력에 넣지 않는다: ", paste(bad, collapse = ", "))
  invisible(TRUE)
}

#' reference 의 형태 단계가 고정 순서 그대로인지
.nso_check_reference_order <- function(reference) {
  st <- reference$stages
  if (is.null(st) || !identical(as.character(st$stage_id), NSO_STAGES) || !identical(as.integer(st$stage_rank), seq_along(NSO_STAGES)) || any(diff(st$source_stage_order) <= 0))
    stop("reference 의 형태 단계 순서가 Myeloblast → PM → MY → MM → BN → SN 이 아니다.")
  invisible(TRUE)
}

#' p92 · p51 세포표(NRV cells_lossless$cells), raw obs ID, graph ID 를 대조한다. 재정렬하지 않고, 어긋나면 멈춘다.
#' 반환: list(checks = data.frame(check, pass, detail), cell_id = canonical ID, cells = data.frame(cell_id, cluster92, cluster51, block, study))
nso_check_inputs <- function(cells92, cells51, raw_ids, graph_ids, reference) {
  rows <- list(); add <- function(check, detail) rows[[length(rows) + 1L]] <<- data.frame(check = check, pass = TRUE, detail = detail, stringsAsFactors = FALSE)
  need <- c("cell", "cell_id", "cluster", "block", "study")
  for (p in c("p92", "p51")) { d <- if (p == "p92") cells92 else cells51
    mdu_check_fit_columns(names(d)); .nso_check_display_columns(names(d))
    miss <- setdiff(need, names(d)); if (length(miss)) stop(p, " 세포표에 필수 열이 없다: ", paste(miss, collapse = ", "))
    if (!identical(as.integer(d$cell), seq_len(nrow(d)))) stop(p, " 세포표의 cell 번호가 1..n 순서가 아니다.")
    if (anyNA(d$cluster) || any(!nzchar(as.character(d$cluster)))) stop(p, " 세포표의 군집에 결측이 있다.") }
  add("no_author_or_display_columns", "p92 / p51 cell tables carry no author, label, DAG, ontology or UMAP columns")
  ids <- as.character(cells92$cell_id)
  ct_check_ids(ids, cells51$cell_id); ct_check_ids(ids, raw_ids); ct_check_ids(ids, graph_ids)
  add("cell_ids_identical_order", sprintf("p92, p51, raw obs and graph IDs: %d cells, no duplicates, identical order (no reordering)", length(ids)))
  add("cell_index_sequential", "cell column equals 1..n in both partitions")
  if (!identical(as.character(cells92$block), as.character(cells51$block))) stop("p92 와 p51 세포표의 block 이 다르다.")
  if (!identical(as.character(cells92$study), as.character(cells51$study))) stop("p92 와 p51 세포표의 study 가 다르다.")
  add("block_study_match_partitions", "block and study columns identical between p92 and p51 cell tables")
  add("clusters_complete", sprintf("no missing cluster: p92 %d clusters, p51 %d clusters", length(unique(cells92$cluster)), length(unique(cells51$cluster))))
  .nso_check_reference_order(reference); add("reference_stage_order", paste(NSO_STAGES, collapse = " > "))
  list(checks = do.call(rbind, rows), cell_id = ids,
       cells = data.frame(cell_id = ids, cluster92 = as.character(cells92$cluster), cluster51 = as.character(cells51$cluster), block = as.character(cells92$block),
                          study = as.character(cells92$study), stringsAsFactors = FALSE))
}


# ── 3. 계산용 메타데이터 ─────────────────────────────────────────────────────
#' 계산 함수에 넘기는 세포 정보: cell_id · cluster92 · cluster51 · donor · study · block 만. 저자 라벨 · UMAP 열이 들어오면 멈춘다.
#' block 은 R24 nrv_blocks(donor) 매핑(반복 표본을 같은 block 으로)과 같아야 한다.
nso_calculation_metadata <- function(cells) {
  mdu_check_fit_columns(names(cells)); .nso_check_display_columns(names(cells))
  miss <- setdiff(NSO_META_COLS, names(cells)); if (length(miss)) stop("계산용 메타데이터에 필수 열이 없다: ", paste(miss, collapse = ", "))
  out <- as.data.frame(lapply(cells[NSO_META_COLS], as.character), stringsAsFactors = FALSE)
  if (anyDuplicated(out$cell_id)) stop("세포 ID 에 중복이 있다.")
  if (any(vapply(out, function(x) anyNA(x) || any(!nzchar(x)), logical(1)))) stop("계산용 메타데이터에 결측이 있다.")
  b <- nrv_blocks(out$donor, require_all = FALSE)
  if (!identical(b, out$block)) stop("block 이 R24 donor → block 매핑과 다르다: ", paste(utils::head(unique(out$block[b != out$block]), 5), collapse = ", "))
  out
}
