# =============================================================================
# build_rna_dag_review.R — RNA 중심 조혈 DAG 검토안: 원본 76 노드 / 89 연결 대응표와 제안 view (사용자 승인 대기)
# =============================================================================
# 입력 (--source, 기본 data/priors/rna_dag_review/source_v1_1_0): 문헌 기반 v1.1.0 의 작은 표 복사본
# 출력 (--out, 기본 data/priors/rna_dag_review/proposal)
#   source_node_map.tsv  원본 76 노드 -> 처리 (keep / split / reference_only / held) · 제안 ID · 이유 · RNA 구별 상태 · 근거 ID
#   proposed_nodes.tsv   제안 77 항목 (RNA 상태 후보 63 · assay 정의 7 · 산물 1 · 보류 2 · taxonomy 4). HSC_candidate · MPP_candidate 는 별도 행
#   source_edge_map.tsv  원본 89 연결 -> 처리 (keep / unresolved / reference_only / held) · 제안 연결 ID · 이유
#   proposed_edges.tsv   RNA 상태 후보 view 에 그리는 연결만 (유지 제안 + 출발점 미해결)
#   --checklist (기본 docs/rna_dag_review/approval_checklist.tsv): 사용자 확인 3 항목. 선택 · 메모 칸은 비워 둔다
# 원칙
#   - 모든 user_approval 은 pending. HSC · MPP 를 별도 후보로 둔다는 사용자 결정만 split_decision 에 적는다 (연결 · 마커 승인이 아님).
#   - HSPC_continuum 출발 19 개 coarse 관계는 출발점 미해결로 두고 두 후보에 복제하지 않는다.
#   - assay 정의 · 산물 · taxonomy 는 참고 항목이다. 추론용 연결에 섞지 않는다. 원본 layer · adult_support 는 덮어쓰지 않는다.
#   - 새 마커 · 생물학적 근거를 만들지 않는다. RNA 구별은 평가하지 않았으므로 not_assessed (HSC · MPP 는 unresolved).
# 사용: Rscript R/rna_dag_review/build_rna_dag_review.R [--source DIR] [--out DIR] [--checklist FILE]
# =============================================================================

RDR_SOURCE_SHA256 <- c(nodes.tsv = "438fd18aa47a5c201f993f3f09ee018eb980e3ee06f118438c567fa718a81eca",
                       edges.tsv = "ccc67f11019d7107adfa1a1c967f1d6bca1a5b038cfe19253f55f43c82fc8004")
RDR_SOURCE_FILES <- c(nodes = "nodes.tsv", edges = "edges.tsv", edge_evidence = "edge_evidence.tsv", taxonomy = "taxonomy.tsv",
                      gate_membership = "gate_membership.tsv", evidence = "evidence_public_subset.tsv")
RDR_SPLIT_SOURCE <- "HSPC_continuum"
RDR_SPLIT_IDS <- c("HSC_candidate", "MPP_candidate")
RDR_UNRESOLVED_ORIGIN <- "unresolved_HSC_or_MPP"
RDR_BADGE_PROTEIN_FUNCTION <- c("EoP", "eNeP", "cMoP", "MCp", "NKP", "ILCP", "NK_CD56bright", "NK_CD56dim")   # 계획에 명시된 단백질 · 기능 정의 의존 후보
RDR_RNA_KINDS <- c("precursor_state", "stage_state", "mature_state", "activation_state")
RDR_NODE_KINDS <- c(RDR_RNA_KINDS, "stem_state", "assay_gate", "product", "held_candidate", "taxonomy_category")
RDR_REFERENCE_ROLES <- c("assay_defined_reference", "product_reference", "taxonomy_reference")
RDR_NODE_COLS <- c("node_id", "label", "module", "compartment", "node_kind", "layer", "layout_x", "layout_y", "definition_ko")
RDR_EDGE_COLS <- c("edge_id", "source", "target", "relation", "directness", "support_class", "layer", "adult_support", "evidence_basis", "endpoint_mapping")


# ── 읽기 · 쓰기 ──────────────────────────────────────────────────────────────
rdr_read_tsv <- function(path) utils::read.delim(path, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
                                                 colClasses = "character", na.strings = character(0), fileEncoding = "UTF-8")

#' NA · 탭 · 줄바꿈이 없는 문자 표만 쓴다 (빈 칸은 ""). 줄 끝은 "\n", UTF-8.
rdr_write_tsv <- function(x, path) {
  cells <- unlist(x, use.names = FALSE)
  if (anyNA(cells) || any(grepl("[\t\n\r]", c(cells, names(x))))) stop("cannot write ", basename(path), ": NA, tab or newline inside a cell")
  lines <- c(paste(names(x), collapse = "\t"), if (nrow(x)) do.call(paste, c(unname(as.list(x)), sep = "\t")))
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
  invisible(path)
}

#' 원본 복사본을 읽는다. nodes · edges 는 v1.1.0 SHA-256 과 같아야 한다.
rdr_read_source <- function(dir) {
  miss <- RDR_SOURCE_FILES[!file.exists(file.path(dir, RDR_SOURCE_FILES))]
  if (length(miss)) stop("source copy lacks: ", paste(miss, collapse = ", "), " (run R/rna_dag_review/copy_source_v1_1_0.R first)")
  for (f in names(RDR_SOURCE_SHA256)) {
    h <- digest::digest(file = file.path(dir, f), algo = "sha256")
    if (!identical(h, RDR_SOURCE_SHA256[[f]])) stop(f, " SHA-256 differs from v1.1.0: ", h)
  }
  lapply(RDR_SOURCE_FILES, function(f) rdr_read_tsv(file.path(dir, f)))
}

.rdr_join <- function(x) { x <- unique(x[nzchar(x)]); if (length(x)) paste(x, collapse = ";") else "none" }
.rdr_split <- function(s) if (identical(s, "none") || !nzchar(s)) character(0) else strsplit(s, ";", fixed = TRUE)[[1]]


# ── 입력 검증 (조용히 넘어가지 않는다) ─────────────────────────────────────────
rdr_validate_source <- function(src) {
  n <- src$nodes; e <- src$edges
  lack <- c(setdiff(RDR_NODE_COLS, names(n)), setdiff(RDR_EDGE_COLS, names(e)))
  if (length(lack)) stop("source tables lack columns: ", paste(lack, collapse = ", "))
  d <- unique(n$node_id[duplicated(n$node_id)]); if (length(d)) stop("duplicate node_id: ", paste(d, collapse = ", "))
  d <- unique(e$edge_id[duplicated(e$edge_id)]); if (length(d)) stop("duplicate edge_id: ", paste(d, collapse = ", "))
  k <- !n$node_kind %in% RDR_NODE_KINDS
  if (any(k)) stop("unknown node_kind: ", paste(sprintf("%s (%s)", n$node_kind[k], n$node_id[k]), collapse = ", "))
  if (!identical(n$node_id[n$node_kind == "stem_state"], RDR_SPLIT_SOURCE)) stop("expected exactly one stem_state node named ", RDR_SPLIT_SOURCE)
  bad <- !(e$source %in% n$node_id & e$target %in% n$node_id)
  if (any(bad)) stop("edges reference unknown nodes: ", paste(sprintf("%s %s->%s", e$edge_id[bad], e$source[bad], e$target[bad]), collapse = ", "))
  s <- e$source == e$target; if (any(s)) stop("self-loop edge: ", paste(e$edge_id[s], collapse = ", "))
  if (any(e$target == RDR_SPLIT_SOURCE)) stop("edge into ", RDR_SPLIT_SOURCE, ": ", paste(e$edge_id[e$target == RDR_SPLIT_SOURCE], collapse = ", "))
  g <- igraph::graph_from_data_frame(e[, c("source", "target")], vertices = n["node_id"])
  if (!igraph::is_dag(g)) stop("source graph has a cycle")
  ee <- src$edge_evidence
  if (!all(ee$edge_id %in% e$edge_id)) stop("edge_evidence references unknown edges: ", paste(setdiff(ee$edge_id, e$edge_id), collapse = ", "))
  if (!all(ee$evidence_id %in% src$evidence$evidence_id)) stop("edge_evidence references unknown evidence: ", paste(setdiff(ee$evidence_id, src$evidence$evidence_id), collapse = ", "))
  tx <- src$taxonomy; if (!all(c(tx$child, tx$parent) %in% n$node_id)) stop("taxonomy references unknown nodes")
  gm <- src$gate_membership; if (!all(c(gm$subset_node, gm$gate_node) %in% n$node_id)) stop("gate_membership references unknown nodes")
  invisible(TRUE)
}


# ── 노드 ─────────────────────────────────────────────────────────────────────
#' 원본 노드의 제안 역할. assay 정의 = assay_gate 이거나 원본 이름에 "(assay-defined)" 가 있는 노드 (EryP).
rdr_node_roles <- function(nodes) {
  k <- nodes$node_kind
  assay <- k == "assay_gate" | grepl("(assay-defined)", nodes$label, fixed = TRUE)
  role <- ifelse(k == "stem_state", "split_source", ifelse(assay, "assay_defined_reference", ifelse(k == "product", "product_reference",
          ifelse(k == "held_candidate", "held_candidate", ifelse(k == "taxonomy_category", "taxonomy_reference", "rna_state_candidate")))))
  stats::setNames(role, nodes$node_id)
}

.RDR_NODE_TEXT <- list(
  treatment = c(split_source = "split", rna_state_candidate = "keep", assay_defined_reference = "reference_only", product_reference = "reference_only",
                taxonomy_reference = "reference_only", held_candidate = "held"),
  reason = c(split_source = "HSC · MPP 를 별도 후보로 두기로 한 사용자 결정(2026-10-16)에 따라 HSC_candidate · MPP_candidate 로 나눔. 원본 상자는 provenance 로 보존",
             rna_state_candidate = "원본 상태 노드를 RNA 상태 후보로 유지 (원본 정의 · ID 그대로)",
             assay_defined_reference = "assay 정의 항목 (FACS gate 또는 colony assay). 삭제하지 않고 RNA 후보 본체와 분리한 참고 panel 에 표시",
             product_reference = "생산 산물 (produces). 세포 상태 후보와 분리해 참고로 표시",
             taxonomy_reference = "분류 상위어. 분화 연결 탐색에서 빼고 참고 panel 에 표시",
             held_candidate = "원본에서 보류(held) 된 후보. 보류 유지"),
  rna = c(split_source = "unresolved", rna_state_candidate = "not_assessed", assay_defined_reference = "not_assessed", product_reference = "not_assessed",
          taxonomy_reference = "not_applicable", held_candidate = "not_assessed"),
  badge = c(split_source = "split_from_HSPC_continuum", rna_state_candidate = "none", assay_defined_reference = "assay_defined",
            product_reference = "production_product", taxonomy_reference = "taxonomy_category", held_candidate = "held_in_source"),
  uncertainty = c(rna_state_candidate = "RNA 상태 구별은 이번 검토에서 평가하지 않음 (원본 정의 유지)",
                  protein_or_function_defined = "원본 정의가 단백질 표면항원 · 기능/배양 측정에 의존. RNA 만으로 구별되는지는 평가하지 않음",
                  split = "HSC · MPP 를 RNA 로 구별할 수 있는지와 하류 연결의 출발점이 미해결. 마커 · 연결은 승인되지 않음",
                  assay_defined_reference = "assay 로 정의된 범주. 참고 항목이며 RNA 상태 구별은 평가하지 않음",
                  product_reference = "무핵 산물. 세포 상태 전환과 구별",
                  taxonomy_reference = "분류 상위어. 세포 상태가 아님",
                  held_candidate = "원본에서 정의 · 기원이 확정되지 않아 보류"))

.rdr_definition_evidence <- function(def, ids) .rdr_join(ids[vapply(ids, function(i) grepl(paste0("\\b", i, "\\b"), def, perl = TRUE), NA)])

.rdr_incident_evidence <- function(node, edges, ee) {
  e <- edges$edge_id[edges$source == node | edges$target == node]
  .rdr_join(sort(unique(ee$evidence_id[ee$edge_id %in% e]), method = "radix"))
}

rdr_source_node_map <- function(src, role) {
  n <- src$nodes; ids <- src$evidence$evidence_id; r <- unname(role[n$node_id]); tt <- .RDR_NODE_TEXT
  badge <- tt$badge[r]; badge[r == "rna_state_candidate" & n$node_id %in% RDR_BADGE_PROTEIN_FUNCTION] <- "protein_or_function_defined"
  data.frame(source_node_id = n$node_id, source_label = n$label, source_module = n$module, source_node_kind = n$node_kind, source_layer = n$layer,
             source_definition_ko = n$definition_ko, treatment = unname(tt$treatment[r]),
             proposed_ids = ifelse(r == "split_source", paste(RDR_SPLIT_IDS, collapse = ";"), n$node_id),
             view_role = ifelse(r == "split_source", "rna_state_candidate", r), reason_ko = unname(tt$reason[r]),
             rna_distinction_status = unname(tt$rna[r]), badges = unname(badge),
             evidence_ids_definition = vapply(n$definition_ko, .rdr_definition_evidence, "", ids = ids, USE.NAMES = FALSE),
             evidence_ids_incident_edges = vapply(n$node_id, .rdr_incident_evidence, "", edges = src$edges, ee = src$edge_evidence, USE.NAMES = FALSE),
             user_approval = "pending", stringsAsFactors = FALSE)
}

rdr_proposed_nodes <- function(src, nm) {
  n <- src$nodes[match(nm$source_node_id, src$nodes$node_id), ]; tt <- .RDR_NODE_TEXT
  kept <- nm$treatment != "split"
  unc <- ifelse(nm$view_role == "rna_state_candidate", tt$uncertainty[["rna_state_candidate"]], tt$uncertainty[nm$view_role])
  unc[nm$badges == "protein_or_function_defined"] <- tt$uncertainty[["protein_or_function_defined"]]
  base <- data.frame(proposed_id = nm$source_node_id, label = n$label, view_role = nm$view_role, node_kind = n$node_kind, module = n$module,
                     compartment = n$compartment, source_node_ids = nm$source_node_id, source_layer = n$layer, definition_ko = n$definition_ko,
                     badges = nm$badges, rna_distinction_status = nm$rna_distinction_status, uncertainty_ko = unname(unc),
                     evidence_ids_source_definition = nm$evidence_ids_definition, split_decision = "not_applicable", user_approval = "pending",
                     stringsAsFactors = FALSE)[kept, ]
  h <- n[!kept, ]; hm <- nm[!kept, ]
  split <- data.frame(proposed_id = RDR_SPLIT_IDS, label = c("HSC candidate", "MPP candidate"), view_role = "rna_state_candidate",
                      node_kind = "hspc_candidate_state", module = h$module, compartment = h$compartment, source_node_ids = h$node_id,
                      source_layer = h$layer,
                      definition_ko = sprintf("원본 %s(HSC 와 MPP 를 하나로 묶은 요약 상자)에서 나눈 %s 후보. 이번 개정은 후보를 나누는 것만 반영하며 RNA 상태 구별 · 하류 연결 · 마커는 정하지 않는다",
                                              h$node_id, c("HSC", "MPP")),
                      badges = hm$badges, rna_distinction_status = "unresolved", uncertainty_ko = tt$uncertainty[["split"]],
                      evidence_ids_source_definition = hm$evidence_ids_definition, split_decision = "confirmed_separate_candidates",
                      user_approval = "pending", stringsAsFactors = FALSE)
  out <- rbind(split, base); rownames(out) <- NULL; out
}


# ── 연결 ─────────────────────────────────────────────────────────────────────
rdr_source_edge_map <- function(src, role) {
  e <- src$edges; ee <- src$edge_evidence; rs <- unname(role[e$source]); rt <- unname(role[e$target])
  tr <- ifelse(e$support_class == "held", "held", ifelse(e$source == RDR_SPLIT_SOURCE, "unresolved",
        ifelse(rs %in% RDR_REFERENCE_ROLES | rt %in% RDR_REFERENCE_ROLES, "reference_only", "keep")))
  odd <- tr == "keep" & !(rs == "rna_state_candidate" & rt == "rna_state_candidate")
  if (any(odd)) stop("non-held edges touching held or unexpected nodes: ", paste(e$edge_id[odd], collapse = ", "))
  odd <- tr == "unresolved" & !rt %in% c("rna_state_candidate", "assay_defined_reference")
  if (any(odd)) stop("HSPC-start edges to unexpected targets: ", paste(e$edge_id[odd], collapse = ", "))
  product <- tr == "reference_only" & (rs == "product_reference" | rt == "product_reference")
  panel <- ifelse(tr == "held", "held", ifelse(tr == "reference_only", ifelse(product, "product_reference", "assay_reference"),
           ifelse(tr == "unresolved" & rt == "assay_defined_reference", "assay_reference", "rna_state")))
  ref_ends <- vapply(seq_len(nrow(e)), function(i) paste(c(e$source[i], e$target[i])[c(rs[i], rt[i]) %in% RDR_REFERENCE_ROLES], collapse = " · "), "")
  reason <- ifelse(tr == "keep", "두 끝점 모두 RNA 상태 후보. 원본 방향 근거를 그대로 둔 제안 연결 (사용자 승인 대기)",
            ifelse(tr == "held", "원본 support_class = held. 보류 유지; 추론용 아님",
            ifelse(tr == "unresolved", ifelse(panel == "rna_state",
                   "원본 HSPC_continuum 출발 coarse 관계. HSC_candidate · MPP_candidate 중 출발점을 정하지 않았고 두 후보에 복제하지 않음",
                   "원본 HSPC_continuum 출발 coarse 관계. 출발점 미해결이며 도착점이 assay 정의 항목이라 참고 panel 에만 표시"),
            ifelse(product, "Platelet 생산(produces) 관계. 세포 상태 전환과 구별해 참고로만 표시; 추론용 아님",
                   paste0("끝점 ", ref_ends, " 이(가) assay 정의 항목. 삭제하지 않고 참고 panel 에만 표시; 추론용 아님")))))
  drawn <- tr == "keep" | (tr == "unresolved" & panel == "rna_state")
  out <- data.frame(source_edge_id = e$edge_id, e[, setdiff(RDR_EDGE_COLS, "edge_id")],
                    source_evidence_ids = vapply(e$edge_id, function(id) .rdr_join(ee$evidence_id[ee$edge_id == id]), "", USE.NAMES = FALSE),
                    treatment = tr, view_panel = panel, origin_attachment = ifelse(tr == "unresolved", RDR_UNRESOLVED_ORIGIN, "unchanged_source_ids"),
                    proposed_edge_ids = ifelse(drawn, paste0("PE_", e$edge_id), "none"),
                    inference_use = ifelse(tr == "keep", "candidate_pending_approval", "not_for_inference"), reason_ko = reason,
                    user_approval = "pending", stringsAsFactors = FALSE, check.names = FALSE)
  rownames(out) <- NULL; out
}

#' 한 연결의 근거 종류를 direction_evidence_ids 와 같은 순서 · 같은 개수로 적는다 (한 문헌에 종류가 여럿이면 "+" 로 잇는다).
.rdr_aligned_kinds <- function(rows) {
  ids <- unique(rows$evidence_id)
  if (!length(ids)) return("none")
  paste(vapply(ids, function(i) paste(unique(rows$evidence_kind[rows$evidence_id == i]), collapse = "+"), "", USE.NAMES = FALSE), collapse = ";")
}

rdr_proposed_edges <- function(src, em) {
  ee <- src$edge_evidence; d <- em[em$proposed_edge_ids != "none", ]; u <- d$treatment == "unresolved"
  out <- data.frame(proposed_edge_id = d$proposed_edge_ids, from_id = ifelse(u, RDR_UNRESOLVED_ORIGIN, d$source),
                    from_candidates = ifelse(u, paste(RDR_SPLIT_IDS, collapse = ";"), d$source), to_id = d$target, source_edge_ids = d$source_edge_id,
                    edge_view = ifelse(u, "unresolved_origin", "state_relation_candidate"),
                    source_relation = d$relation, source_directness = d$directness, source_support_class = d$support_class, source_layer = d$layer,
                    source_adult_support = d$adult_support, source_evidence_basis = d$evidence_basis,
                    direction_evidence_ids = d$source_evidence_ids,
                    direction_evidence_kinds = vapply(d$source_edge_id, function(id) .rdr_aligned_kinds(ee[ee$edge_id == id, ]), "", USE.NAMES = FALSE),
                    rna_distinction_evidence = ifelse(u, "unresolved", "not_assessed"),
                    review_status = ifelse(u, "unresolved_origin", "pending_user_review"), inference_use = d$inference_use,
                    user_approval = "pending", stringsAsFactors = FALSE)
  rownames(out) <- NULL; out
}


# ── 사용자 확인 3 항목 ───────────────────────────────────────────────────────
rdr_checklist <- function(pn, em) {
  r <- table(factor(pn$view_role, c("rna_state_candidate", "assay_defined_reference", "product_reference", "held_candidate", "taxonomy_reference")))
  t <- table(factor(em$treatment, c("keep", "unresolved", "reference_only", "held")))
  hspc <- em$source_edge_id[em$treatment == "unresolved"]
  ref <- pn$proposed_id[pn$view_role %in% RDR_REFERENCE_ROLES]
  data.frame(
    item_id = c("A1", "A2", "A3"),
    question_ko = c("RNA 중심 DAG 개정안 전체를 받아들일까요?",
                    sprintf("HSC/MPP 에서 출발하는 원본 연결 %d 개의 출발점을 어떻게 둘까요?", length(hspc)),
                    "assay 정의 · 산물 · taxonomy 항목을 RNA 상태 후보 본체와 분리해 둘까요?"),
    proposal_ko = c(sprintf("제안 %d 항목 (RNA 상태 후보 %d · assay 정의 %d · 산물 %d · 보류 %d · taxonomy %d). 원본 연결 %d 개는 유지 제안 %d · 출발점 미해결 %d · 참고 %d · 보류 %d 로 대응",
                            nrow(pn), r[[1]], r[[2]], r[[3]], r[[4]], r[[5]], nrow(em), t[["keep"]], t[["unresolved"]], t[["reference_only"]], t[["held"]]),
                    "모두 출발점 미해결로 두고 HSC_candidate · MPP_candidate 어느 쪽에도 붙이지 않는다. 원본에 없는 HSC-MPP 사이 연결은 그리지 않는다",
                    "6 개 gate 와 EryP(BFU-E / CFU-E) 는 삭제하지 않고 참고 panel 로 분리한다. Platelet 은 산물, taxonomy 4 개는 분류 상위어로 표시한다"),
    options_ko = c("(a) 그대로 승인 / (b) 수정 후 승인: 고칠 노드 · 연결 ID 와 내용을 user_notes_ko 에 적기 / (c) 승인하지 않음: 원본 v1.1.0 유지",
                   "(a) 미해결 유지 / (b) 수정: 연결별 출발점(HSC_candidate 또는 MPP_candidate) 지정이나 HSC-MPP 연결 추가를 edge ID 와 함께 user_notes_ko 에 적기 / (c) 원본처럼 하나의 HSPC 출발점으로 되돌림",
                   "(a) 분리 유지 / (b) 수정: RNA 상태 후보로 되돌릴 ID 를 user_notes_ko 에 적기 / (c) 원본 배치로 되돌림"),
    scope_ids = c("all", paste(hspc, collapse = ";"), paste(ref, collapse = ";")),
    user_choice = "", user_notes_ko = "", user_approval = "pending", stringsAsFactors = FALSE)
}


# ── 조립 · 불변식 검사 ───────────────────────────────────────────────────────
rdr_check_review <- function(rv, src) {
  nm <- rv$source_node_map; pn <- rv$proposed_nodes; em <- rv$source_edge_map; pe <- rv$proposed_edges
  if (!setequal(nm$source_node_id, src$nodes$node_id) || nrow(nm) != nrow(src$nodes)) stop("source_node_map does not trace every source node once")
  if (!setequal(em$source_edge_id, src$edges$edge_id) || nrow(em) != nrow(src$edges)) stop("source_edge_map does not trace every source edge once")
  if (anyDuplicated(pn$proposed_id)) stop("duplicate proposed_id")
  ids <- unique(unlist(lapply(nm$proposed_ids, .rdr_split)))
  if (!setequal(ids, pn$proposed_id)) stop("node map and proposed nodes disagree: ", paste(union(setdiff(ids, pn$proposed_id), setdiff(pn$proposed_id, ids)), collapse = ", "))
  r <- pe$from_id != RDR_UNRESOLVED_ORIGIN
  if (!all(c(pe$from_id[r], pe$to_id) %in% pn$proposed_id)) stop("proposed edges reference unknown proposed nodes")
  if (any(c(pe$from_id, pe$to_id) %in% RDR_SPLIT_IDS)) stop("an edge was attached to HSC_candidate or MPP_candidate without approval")
  if (anyDuplicated(pe$source_edge_ids)) stop("a source edge was copied into more than one proposed edge")
  g <- igraph::graph_from_data_frame(pe[r, c("from_id", "to_id")], vertices = data.frame(id = pn$proposed_id))
  if (!igraph::is_dag(g)) stop("proposed view has a cycle")
  if (!all(c(nm$user_approval, pn$user_approval, em$user_approval, pe$user_approval, rv$checklist$user_approval) == "pending")) stop("user_approval must stay pending")
  invisible(TRUE)
}

#' 원본 표 목록 -> 대응표 · 제안 view · 확인 항목.
rdr_build_review <- function(src) {
  rdr_validate_source(src)
  role <- rdr_node_roles(src$nodes)
  nm <- rdr_source_node_map(src, role); pn <- rdr_proposed_nodes(src, nm)
  em <- rdr_source_edge_map(src, role); pe <- rdr_proposed_edges(src, em)
  rv <- list(source_node_map = nm, proposed_nodes = pn, source_edge_map = em, proposed_edges = pe, checklist = rdr_checklist(pn, em))
  rdr_check_review(rv, src)
  rv
}

RDR_REVIEW_TABLES <- c("source_node_map", "proposed_nodes", "source_edge_map", "proposed_edges")

#' 이미 있는 출력 파일의 사용자 검토 기록: user_choice · user_notes_ko 가 비어 있지 않거나 user_approval 이 pending 이 아닌 열 이름.
#' 열 구성이 이번 출력과 다른 파일은 무엇인지 확인할 수 없으므로 기록이 있는 것으로 본다.
.rdr_review_entries <- function(path, expected_cols) {
  if (!file.exists(path)) return(character(0))
  x <- tryCatch(rdr_read_tsv(path), error = function(e) NULL)
  if (is.null(x) || !identical(names(x), expected_cols)) return(paste0(basename(path), " (unexpected layout)"))
  hit <- c(intersect(c("user_choice", "user_notes_ko"), names(x))[vapply(intersect(c("user_choice", "user_notes_ko"), names(x)), function(k) any(nzchar(x[[k]])), NA)],
           if ("user_approval" %in% names(x) && any(x$user_approval != "pending")) "user_approval")
  if (length(hit)) paste0(basename(path), " (", paste(hit, collapse = ", "), ")") else character(0)
}

#' 대응표 4 개와 확인 항목을 쓴다. 기존 출력에 사용자 검토 기록이 있으면 하나도 쓰지 않고 멈춘다 (자동 승인 · 초기화 없음).
rdr_write_review <- function(rv, out_dir, checklist_path) {
  paths <- c(file.path(out_dir, paste0(RDR_REVIEW_TABLES, ".tsv")), checklist_path)
  tabs <- c(rv[RDR_REVIEW_TABLES], list(rv$checklist))
  found <- unlist(Map(function(p, x) .rdr_review_entries(p, names(x)), paths, tabs))
  if (length(found)) stop("existing outputs hold user review entries: ", paste(found, collapse = "; "),
                          ". Nothing was written; keep that review and use a new --out / --checklist path.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(checklist_path), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(paths)) rdr_write_tsv(tabs[[i]], paths[i])
  invisible(paths)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  arg <- function(k, default) { i <- match(k, a); if (is.na(i)) default else a[i + 1L] }
  src <- rdr_read_source(arg("--source", file.path("data", "priors", "rna_dag_review", "source_v1_1_0")))
  rv <- rdr_build_review(src)
  out <- rdr_write_review(rv, arg("--out", file.path("data", "priors", "rna_dag_review", "proposal")),
                          arg("--checklist", file.path("docs", "rna_dag_review", "approval_checklist.tsv")))
  cat("proposed items by view_role:\n"); print(table(rv$proposed_nodes$view_role))
  cat("source edges by treatment:\n"); print(table(rv$source_edge_map$treatment))
  cat("written:\n"); cat(paste0("  ", out), sep = "\n")
}
