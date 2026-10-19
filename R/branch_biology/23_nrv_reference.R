# =============================================================================
# 23_nrv_reference.R — 참고표를 읽고, 층·근거·유전자 역할을 실행 전에 고정한다
# =============================================================================
#
# 무엇을 하나
#   사용자가 제공한 참고 CSV(02_nodes, 03_edges)를 읽어 다음을 만든다.
#   - 문헌 edge 의 유전자를 core / conditional, up / down 으로 나눈 긴 표 (빈칸 = 근거 미지정)
#   - node 의 marker 표
#   - 층(성인 확인 / 연령 미확인 / 범위 밖)과 근거 종류(고전 단계 투영 / 프로그램 파동 / 직접 인접 비교)
#   - 역할 배정 A/B/C: anchor 유전자와 그 밖의 판정용 유전자를 나눈다.
#     anchor 에 쓴 유전자의 증감은 판정에 쓰지 않는다(anchor_circular).
#   새 marker 규칙을 만들지 않는다. 참고표 원본은 읽기만 한다.
#
# 주요 함수
#   nrv_read_reference(dir)        참고표 읽기
#   nrv_gene_roles(ref, allocation) 역할 배정
#   nrv_programs(roles)            판정용 프로그램(부호 포함)
#   nrv_check_no_leak(roles, used) anchor 유전자가 판정에 쓰이면 멈춘다
#   nrv_stress_audit()             스트레스·염증 근거 점검 결과
# =============================================================================

NRV_REF_DIR <- "/Users/edemhw/.codex/visualizations/2026/10/03/01a0ffdf-e36a-70b0-b52b-3389d008cd25/neutrophil_expression/cytoscape"

# 층: 성인 확인(E01–E03), 연령 미확인 사람 근거(E04–E09), 범위 밖(E10: 골수 → 혈액 이동)
NRV_EDGE_LAYER <- c(E01 = "adult_confirmed", E02 = "adult_confirmed", E03 = "adult_confirmed", E04 = "human_age_unconfirmed",
                    E05 = "human_age_unconfirmed", E06 = "human_age_unconfirmed", E07 = "human_age_unconfirmed",
                    E08 = "human_age_unconfirmed", E09 = "human_age_unconfirmed", E10 = "out_of_scope")

# 근거 종류: 참고표 direction_basis 를 그대로 대응한다 (모르는 값은 멈춘다)
NRV_SUPPORT_OF_BASIS <- c(adult_human_RNA_programs_projected_to_classical_edge = "classical_stage_projection",
                          candidate_mapping_between_progenitor_and_morphology = "classical_stage_projection",
                          granule_wave_with_conditional_decline = "program_wave", granule_wave_literature = "program_wave",
                          maturation_program_literature = "program_wave", granule_wave_and_candidate_receptor_trend = "program_wave",
                          reported_decline_and_candidate_receptor_trend = "direct_adjacent_comparison",
                          marrow_blood_transcriptome_comparison = "direct_adjacent_comparison")

# 역할 배정 (실행 전 고정, 셋 다 보고하고 합치지 않는다)
NRV_ALLOCATIONS <- list(
  A = list(S = c("HLF", "HOPX", "CRHBP", "MLLT3", "MECOM"), T_early = c("MPO", "ELANE", "PRTN3"), T_terminal = c("FCGR3B", "FPR1", "CXCR2")),
  B = list(S = c("HOPX", "MECOM"), T_early = c("MPO", "ELANE", "PRTN3"), T_terminal = c("FCGR3B", "FPR1", "CXCR2")),
  C = list(S = c("HLF", "HOPX", "CRHBP", "MLLT3", "MECOM"), T_early = c("CEBPA", "GFI1"), T_terminal = c("FCGR3B", "FPR1", "CXCR2")))

# 기준 상수: 성격은 statistical / mathematical / operational 중 하나 (생물학적 기준이 아니다)
NRV_THRESHOLDS <- data.frame(
  name = c("set_min_cells", "anchor_block_min_cells", "anchor_min_blocks", "anchor_candidate_freq", "program_min_genes",
           "gene_min_detected_cells", "n_for_strictest_holm_threshold", "alpha", "n_boot", "n_gene_null", "num_tol", "tie_tol",
           "yen_k", "n_folds"),
  value = c("5", "20", "3", "0.10", "2", "5", "ceiling(log2(2*m/alpha)) (설명 열, 거르지 않음)", "0.05", "2000", "1000",
            "MDU_NUM_TOL", "1e-12", "10", "5"),
  nature = c("operational", "operational", "operational", "operational", "operational", "operational", "mathematical",
             "statistical", "statistical", "statistical", "operational", "operational", "operational", "operational"),
  stringsAsFactors = FALSE)

.nrv_split <- function(x) { x <- if (is.na(x)) "" else as.character(x); v <- trimws(strsplit(x, ";", fixed = TRUE)[[1]]); v[nzchar(v)] }

#' direction_basis → 근거 종류
nrv_reference_support <- function(basis) {
  s <- NRV_SUPPORT_OF_BASIS[basis]
  if (anyNA(s)) stop("모르는 direction_basis: ", paste(basis[is.na(s)], collapse = ", "))
  unname(s)
}

#' 참고표 읽기. 반환: list(nodes, edges, gene_long(edge_id, gene, direction, evidence), markers(node_id, gene))
nrv_read_reference <- function(dir) {
  rd <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  nodes <- rd("02_nodes.csv"); edges <- rd("03_edges.csv")
  cols <- list(c("up_genes", "up", "core"), c("down_genes", "down", "core"),
               c("conditional_up_genes", "up", "conditional"), c("conditional_down_genes", "down", "conditional"))
  gl <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) do.call(rbind, lapply(cols, function(cc) {
    g <- .nrv_split(edges[[cc[1]]][i]); if (!length(g)) return(NULL)
    data.frame(edge_id = edges$edge_id[i], gene = g, direction = cc[2], evidence = cc[3], stringsAsFactors = FALSE) }))))
  mk <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
    g <- .nrv_split(nodes$marker_genes[i]); if (!length(g)) return(NULL)
    data.frame(node_id = nodes$node_id[i], gene = g, stringsAsFactors = FALSE) }))
  edges$reference_support <- nrv_reference_support(edges$direction_basis)
  list(nodes = nodes, edges = edges, gene_long = gl, markers = mk)
}

#' 배정의 anchor 유전자 목록
nrv_anchor_genes <- function(allocation) {
  a <- NRV_ALLOCATIONS[[allocation]]; if (is.null(a)) stop("모르는 배정: ", allocation); a
}

#' 역할 배정: anchor 유전자의 증감 행은 anchor_circular, 나머지는 evaluation
nrv_gene_roles <- function(ref, allocation) {
  an <- unique(unlist(nrv_anchor_genes(allocation)))
  r <- ref$gene_long; r$role <- ifelse(r$gene %in% an, "anchor_circular", "evaluation"); r$allocation <- allocation; r
}

#' 판정에 쓴 유전자가 선택에 쓴 유전자와 겹치면 멈춘다 (Codex 재현 수정).
#'  - 배정의 anchor 유전자 전부를 본다. edge 유전자 표에 없는 node marker 전용 유전자(배정 A 의 HOPX, MECOM)도 포함한다.
#'  - selection_genes: 그 조건의 전역 선택 유전자(masked panel, anchor·node 대조 유전자 등)
nrv_check_no_leak <- function(roles, used_genes, selection_genes = character(0)) {
  al <- unique(roles$allocation); if (length(al) != 1L) stop("역할 표의 배정이 하나가 아니다: ", paste(al, collapse = ", "))
  bad <- intersect(used_genes, union(unlist(nrv_anchor_genes(al)), roles$gene[roles$role == "anchor_circular"]))
  if (length(bad)) stop("anchor 유전자가 판정에 쓰였다: ", paste(bad, collapse = ", "))
  bad <- intersect(used_genes, selection_genes)
  if (length(bad)) stop("선택에 쓴 유전자가 판정에 쓰였다: ", paste(bad, collapse = ", "))
  invisible(TRUE)
}

#' 판정용 프로그램: (edge, core/conditional, up/down) 마다 하나. 부호는 up = +1, down = -1. AND 조건은 없다.
nrv_programs <- function(roles) {
  ev <- roles[roles$role == "evaluation", ]
  data.frame(program_id = paste(ev$edge_id, ev$evidence, ev$direction, sep = "_"), edge_id = ev$edge_id, evidence = ev$evidence,
             direction = ev$direction, gene = ev$gene, sign = ifelse(ev$direction == "up", 1, -1), stringsAsFactors = FALSE)
}

#' 스트레스·염증 근거 점검 결과 (2026-10-13 계획 단계 점검; 새 규칙을 만들지 않는다)
nrv_stress_audit <- function() {
  data.frame(
    topic = c("stress", "inflammation", "stress", "inflammation", "gene_set_database", rep("diagnostic_metadata", 6)),
    item = c("reference_csv_6_files", "reference_csv_6_files", "data/priors + differentiation_dag_v1", "R/branch_biology/programs_template.tsv",
             "msigdbr/msigdbdf/GSVA/fgsea in dcta renv", "S.Score", "G2M.Score", "CyclePhase", "nCount_RNA", "nFeature_RNA", "scrublet_scores"),
    status = c(rep("no_approved_program_found", 4), "not_installed_install_forbidden", rep("available_author_metadata_diagnostic_only", 6)),
    stringsAsFactors = FALSE)
}
