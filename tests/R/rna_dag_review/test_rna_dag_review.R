# =============================================================================
# test_rna_dag_review.R — RNA 중심 DAG 검토안 (원본 76 노드 / 89 연결 대응표 · 제안 view · 그림)
# 실행: Rscript tests/R/rna_dag_review/run_tests.R   (저장소 루트 기준 경로만 쓴다)
# =============================================================================

root <- Sys.getenv("RDR_REPO_ROOT", ".")
source(file.path(root, "R", "rna_dag_review", "build_rna_dag_review.R"))
source(file.path(root, "R", "rna_dag_review", "render_rna_dag_review.R"))

src_dir <- file.path(root, "data", "priors", "rna_dag_review", "source_v1_1_0")
src <- rdr_read_source(src_dir)
rv <- rdr_build_review(src)

HSPC_START <- c(sprintf("E%03d", 2:15), sprintf("E%03d", 201:205))
STATUS_COLS <- c("treatment", "review_status", "rna_distinction_status", "rna_distinction_evidence", "inference_use", "user_approval", "user_choice")
OVERCLAIM <- "valid|verif|confirm|approved|accepted|확정|검증됨|승인됨"
PRIVATE <- "/Users/|/home/|claude|codex|session|worktree|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"

split_ids <- function(s) unlist(lapply(s, function(z) if (identical(z, "none")) character(0) else strsplit(z, ";", fixed = TRUE)[[1]]))

# ── 원본 복사본 ──────────────────────────────────────────────────────────────
test_that("source copy matches the recorded and planned hashes", {
  sums <- rdr_read_tsv(file.path(src_dir, "SHA256SUMS.tsv"))
  expect_true(all(c("nodes.tsv", "edges.tsv", "edge_evidence.tsv", "taxonomy.tsv", "gate_membership.tsv", "evidence_public_subset.tsv") %in% sums$file))
  for (i in seq_len(nrow(sums))) if (sums$copy_mode[i] != "withheld")
    expect_identical(digest::digest(file = file.path(src_dir, sums$file[i]), algo = "sha256"), sums$sha256[i], info = sums$file[i])
  same <- sums$copy_mode == "byte_identical"
  expect_identical(sums$sha256[same], sums$source_sha256[same])
  expect_identical(sums$sha256[match(names(RDR_SOURCE_SHA256), sums$file)], unname(RDR_SOURCE_SHA256))
})

test_that("public source copy carries no internal reviewer tags or private paths", {
  for (f in list.files(src_dir, pattern = "\\.tsv$")) {
    x <- readLines(file.path(src_dir, f), warn = FALSE, encoding = "UTF-8")
    expect_false(any(grepl(PRIVATE, x, ignore.case = TRUE, perl = TRUE)), info = f)
  }
  ev <- rdr_read_tsv(file.path(src_dir, "evidence_public_subset.tsv"))
  expect_false("reviewed_by" %in% names(ev))
  expect_identical(nrow(ev), 42L)
})

# ── 76 / 89 추적 ────────────────────────────────────────────────────────────
test_that("all 76 source nodes are traced exactly once", {
  m <- rv$source_node_map
  expect_identical(nrow(m), 76L)
  expect_false(anyDuplicated(m$source_node_id) > 0)
  expect_setequal(m$source_node_id, src$nodes$node_id)
})

test_that("all 89 source edges are traced exactly once", {
  m <- rv$source_edge_map
  expect_identical(nrow(m), 89L)
  expect_false(anyDuplicated(m$source_edge_id) > 0)
  expect_setequal(m$source_edge_id, src$edges$edge_id)
})

# ── 제안 view ───────────────────────────────────────────────────────────────
test_that("proposed view has 77 unique items with the planned role counts", {
  p <- rv$proposed_nodes
  expect_identical(nrow(p), 77L)
  expect_false(anyDuplicated(p$proposed_id) > 0)
  roles <- table(p$view_role)
  expect_identical(as.integer(roles[c("rna_state_candidate", "assay_defined_reference", "product_reference", "held_candidate", "taxonomy_reference")]),
                   c(63L, 7L, 1L, 2L, 4L))
  expect_false("HSPC_continuum" %in% p$proposed_id)
  expect_setequal(p$proposed_id[p$view_role == "assay_defined_reference"], c("CMP_gate", "GMP_gate", "MEP_gate", "CLP_gate", "MDP_gate", "CDP_gate", "EryP"))
})

test_that("HSC and MPP are separate candidate rows split from HSPC_continuum", {
  p <- rv$proposed_nodes
  hm <- p[p$proposed_id %in% c("HSC_candidate", "MPP_candidate"), ]
  expect_identical(nrow(hm), 2L)
  expect_true(all(hm$source_node_ids == "HSPC_continuum"))
  expect_true(all(hm$view_role == "rna_state_candidate"))
  expect_true(all(hm$rna_distinction_status == "unresolved"))
  expect_true(all(hm$split_decision == "confirmed_separate_candidates"))
  expect_true(all(p$split_decision[!p$proposed_id %in% hm$proposed_id] == "not_applicable"))
  h <- rv$source_node_map[rv$source_node_map$source_node_id == "HSPC_continuum", ]
  expect_identical(h$treatment, "split")
  expect_identical(h$proposed_ids, "HSC_candidate;MPP_candidate")
})

test_that("every node mapping resolves in both directions", {
  m <- rv$source_node_map; p <- rv$proposed_nodes
  expect_true(all(split_ids(m$proposed_ids) %in% p$proposed_id))
  expect_true(all(split_ids(p$source_node_ids) %in% src$nodes$node_id))
  expect_setequal(unique(split_ids(m$proposed_ids)), p$proposed_id)
})

test_that("node and edge treatments match the first proposal", {
  nt <- table(rv$source_node_map$treatment)
  expect_identical(as.integer(nt[c("keep", "split", "reference_only", "held")]), c(61L, 1L, 12L, 2L))
  et <- table(rv$source_edge_map$treatment)
  expect_identical(as.integer(et[c("keep", "unresolved", "reference_only", "held")]), c(49L, 19L, 13L, 8L))
  expect_identical(nrow(rv$proposed_edges), 62L)
})

# ── HSPC 출발 19 개 ──────────────────────────────────────────────────────────
test_that("HSPC-start relations stay unresolved and are never copied to both candidates", {
  em <- rv$source_edge_map; pe <- rv$proposed_edges
  expect_setequal(em$source_edge_id[em$source == "HSPC_continuum"], HSPC_START)
  h <- em[em$source_edge_id %in% HSPC_START, ]
  expect_true(all(h$treatment == "unresolved"))
  expect_true(all(h$origin_attachment == "unresolved_HSC_or_MPP"))
  expect_true(all(h$inference_use == "not_for_inference"))
  expect_false(any(c(pe$from_id, pe$to_id) %in% c("HSC_candidate", "MPP_candidate")))
  expect_false(anyDuplicated(split_ids(pe$source_edge_ids)) > 0)
  u <- pe[pe$from_id == "unresolved_HSC_or_MPP", ]
  expect_identical(nrow(u), 13L)
  expect_true(all(u$from_candidates == "HSC_candidate;MPP_candidate"))
  expect_true(all(u$review_status == "unresolved_origin"))
  expect_true(all(split_ids(u$source_edge_ids) %in% HSPC_START))
})

# ── 참고 · 보류 연결 분리 ────────────────────────────────────────────────────
test_that("reference-only and held relations are not mixed into inference edges", {
  em <- rv$source_edge_map; pe <- rv$proposed_edges; p <- rv$proposed_nodes
  role <- stats::setNames(p$view_role, p$proposed_id)
  cand <- pe[pe$inference_use == "candidate_pending_approval", ]
  expect_identical(nrow(cand), 49L)
  expect_true(all(role[cand$from_id] == "rna_state_candidate" & role[cand$to_id] == "rna_state_candidate"))
  side <- em[em$treatment %in% c("reference_only", "held"), ]
  expect_true(all(side$inference_use == "not_for_inference"))
  expect_true(all(side$proposed_edge_ids == "none"))
  expect_false(any(split_ids(pe$source_edge_ids) %in% side$source_edge_id))
  expect_true(all(em$source_edge_id[em$support_class == "held"] %in% em$source_edge_id[em$treatment == "held"]))
})

test_that("proposed edges reference existing nodes and every graph is acyclic", {
  pe <- rv$proposed_edges; p <- rv$proposed_nodes
  r <- pe$from_id != "unresolved_HSC_or_MPP"
  expect_true(all(pe$from_id[r] %in% p$proposed_id))
  expect_true(all(pe$to_id %in% p$proposed_id))
  expect_true(all(split_ids(pe$source_edge_ids) %in% src$edges$edge_id))
  g_src <- igraph::graph_from_data_frame(src$edges[, c("source", "target")], vertices = src$nodes["node_id"])
  expect_true(igraph::is_dag(g_src))
  g_pro <- igraph::graph_from_data_frame(pe[r, c("from_id", "to_id")], vertices = data.frame(id = p$proposed_id))
  expect_true(igraph::is_dag(g_pro))
  side <- rv$source_edge_map[rv$source_edge_map$treatment %in% c("reference_only", "held"), c("source", "target")]
  g_all <- igraph::graph_from_data_frame(rbind(stats::setNames(pe[r, c("from_id", "to_id")], c("source", "target")), side), vertices = data.frame(id = p$proposed_id))
  expect_true(igraph::is_dag(g_all))
})

# ── 승인 · 과장 방지 ─────────────────────────────────────────────────────────
test_that("pending items are never reported as validated", {
  tabs <- list(rv$source_node_map, rv$proposed_nodes, rv$source_edge_map, rv$proposed_edges, rv$checklist)
  for (x in tabs) {
    for (col in intersect(STATUS_COLS, names(x))) expect_false(any(grepl(OVERCLAIM, x[[col]], ignore.case = TRUE)), info = col)
    if ("user_approval" %in% names(x)) expect_true(all(x$user_approval == "pending"))
  }
  expect_true(all(rv$proposed_nodes$rna_distinction_status %in% c("not_assessed", "unresolved", "not_applicable")))
  expect_true(all(rv$source_node_map$rna_distinction_status %in% c("not_assessed", "unresolved", "not_applicable")))
  expect_true(all(rv$proposed_edges$review_status %in% c("pending_user_review", "unresolved_origin")))
})

test_that("source attributes are preserved and not overwritten", {
  em <- rv$source_edge_map; e <- src$edges[match(em$source_edge_id, src$edges$edge_id), ]
  for (col in c("source", "target", "relation", "directness", "support_class", "layer", "adult_support", "evidence_basis"))
    expect_identical(em[[col]], e[[col]], info = col)
  nm <- rv$source_node_map; n <- src$nodes[match(nm$source_node_id, src$nodes$node_id), ]
  expect_identical(nm$source_layer, n$layer)
  expect_identical(nm$source_definition_ko, n$definition_ko)
  pe <- rv$proposed_edges; k <- src$edges[match(pe$source_edge_ids, src$edges$edge_id), ]
  expect_identical(pe$source_layer, k$layer)
  expect_identical(pe$source_adult_support, k$adult_support)
})

test_that("protein- or function-defined candidates carry a badge", {
  p <- rv$proposed_nodes
  b <- p$badges[match(RDR_BADGE_PROTEIN_FUNCTION, p$proposed_id)]
  expect_true(all(grepl("protein_or_function_defined", b, fixed = TRUE)))
  expect_setequal(p$proposed_id[grepl("protein_or_function_defined", p$badges, fixed = TRUE)], RDR_BADGE_PROTEIN_FUNCTION)
})

test_that("direction evidence and RNA distinction evidence are kept separate", {
  pe <- rv$proposed_edges
  ee <- src$edge_evidence
  for (i in seq_len(nrow(pe))) {
    want <- unique(ee$evidence_id[ee$edge_id == pe$source_edge_ids[i]])
    expect_identical(split_ids(pe$direction_evidence_ids[i]), want, info = pe$proposed_edge_id[i])
    kinds <- vapply(want, function(id) paste(unique(ee$evidence_kind[ee$edge_id == pe$source_edge_ids[i] & ee$evidence_id == id]), collapse = "+"), "", USE.NAMES = FALSE)
    expect_identical(split_ids(pe$direction_evidence_kinds[i]), kinds, info = pe$proposed_edge_id[i])
  }
  expect_true(all(pe$rna_distinction_evidence %in% c("not_assessed", "unresolved")))
})

# ── 승인 체크리스트 ──────────────────────────────────────────────────────────
test_that("approval checklist offers at most three concrete choices and is not pre-filled", {
  ck <- rv$checklist
  expect_true(nrow(ck) >= 1L && nrow(ck) <= 3L)
  expect_true(all(ck$user_choice == ""))
  expect_true(all(ck$user_notes_ko == ""))
  expect_true(all(ck$user_approval == "pending"))
  expect_true(all(nzchar(ck$options_ko)))
  expect_true(any(grepl("수정", ck$options_ko)))
})

# ── 쓰기 · 재현 · 실패 ───────────────────────────────────────────────────────
test_that("writing is deterministic and round-trips without NA", {
  a <- tempfile("rdr_a"); b <- tempfile("rdr_b")
  rdr_write_review(rv, a, file.path(a, "approval_checklist.tsv"))
  rdr_write_review(rdr_build_review(rdr_read_source(src_dir)), b, file.path(b, "approval_checklist.tsv"))
  fa <- sort(list.files(a)); expect_identical(fa, sort(list.files(b)))
  for (f in fa) expect_identical(digest::digest(file = file.path(a, f), algo = "sha256"), digest::digest(file = file.path(b, f), algo = "sha256"), info = f)
  back <- rdr_read_tsv(file.path(a, "proposed_edges.tsv"))
  expect_identical(back, rv$proposed_edges)
  for (f in fa) expect_false(any(grepl("(^|\t)NA(\t|$)", readLines(file.path(a, f), encoding = "UTF-8"))), info = f)
})

test_that("existing outputs with user review entries are never overwritten", {
  sha_all <- function(d) vapply(sort(list.files(d, full.names = TRUE)), function(f) digest::digest(file = f, algo = "sha256"), "")
  d <- tempfile("rdr_keep"); ck <- file.path(d, "approval_checklist.tsv")
  rdr_write_review(rv, d, ck)
  x <- rdr_read_tsv(ck); x$user_choice[1] <- "(b)"; x$user_notes_ko[1] <- "E078 재검토"; rdr_write_tsv(x, ck)
  before <- sha_all(d)
  expect_error(rdr_write_review(rv, d, ck), "user review entries")
  expect_identical(sha_all(d), before)
  d2 <- tempfile("rdr_keep2"); ck2 <- file.path(d2, "approval_checklist.tsv")
  rdr_write_review(rv, d2, ck2)
  p <- rdr_read_tsv(file.path(d2, "proposed_edges.tsv")); p$user_approval[1] <- "changed_by_user"; rdr_write_tsv(p, file.path(d2, "proposed_edges.tsv"))
  before2 <- sha_all(d2)
  expect_error(rdr_write_review(rv, d2, ck2), "user review entries")
  expect_identical(sha_all(d2), before2)
  d3 <- tempfile("rdr_new")
  expect_silent(rdr_write_review(rv, d3, file.path(d3, "approval_checklist.tsv")))
  expect_silent(rdr_write_review(rv, d3, file.path(d3, "approval_checklist.tsv")))
})

test_that("builder fails loudly on broken input", {
  bad <- src; bad$edges$target[1] <- "Missing_node"
  expect_error(rdr_build_review(bad), "Missing_node")
  bad <- src; bad$nodes <- rbind(bad$nodes, bad$nodes[2, ])
  expect_error(rdr_build_review(bad), "duplicate")
  bad <- src; bad$nodes$node_kind[bad$nodes$node_id == "MkP"] <- "new_kind"
  expect_error(rdr_build_review(bad), "new_kind")
  bad <- src; bad$edges$source[bad$edges$edge_id == "E002"] <- "MkP"; bad$edges$target[bad$edges$edge_id == "E002"] <- "MkP"
  expect_error(rdr_build_review(bad), "cycle|self")
})

# ── 그림 ─────────────────────────────────────────────────────────────────────
test_that("renderer writes the comparison, core and lineage figures", {
  out <- tempfile("rdr_fig")
  files <- rdr_render_all(rv, src, out)
  expect_true(length(files) >= 16L)
  expect_true(all(file.exists(files)))
  expect_true(all(file.size(files) > 5000))
  expect_true(all(c("rna_dag_original_vs_proposed.png", "rna_dag_core.png", "rna_dag_core.pdf") %in% basename(files)))
  expect_identical(sum(grepl("^rna_dag_lineage_.*\\.png$", basename(files))), length(RDR_LINEAGE_GROUPS))
})

test_that("every source edge appears in at least one lineage figure", {
  shown <- unique(unlist(lapply(names(RDR_LINEAGE_GROUPS), function(g) rdr_lineage_view(rv, src, g)$edges$source_edge_id)))
  expect_setequal(intersect(shown, src$edges$edge_id), src$edges$edge_id)
})

get_specs <- function() c(list(core = rdr_core_spec(rv, src)), lapply(stats::setNames(names(RDR_LINEAGE_GROUPS), names(RDR_LINEAGE_GROUPS)), function(g) rdr_lineage_spec(rv, src, g)))
inside <- function(p, b, pad = 0) any(abs(p[, 1] - b$x) < b$w / 2 - pad & abs(p[, 2] - b$y) < b$h / 2 - pad)

test_that("rendering raises no warnings", {
  expect_no_warning(rdr_render_all(rv, src, tempfile("rdr_fig_w")))
})

test_that("every legend and the whole layout fit inside the canvas", { specs <- get_specs()
  for (nm in names(specs)) {
    s <- specs[[nm]]; b <- s$legend_bbox
    expect_true(b[["x0"]] >= 0 && b[["x1"]] <= s$W && b[["y0"]] >= 0 && b[["y1"]] <= s$H, info = nm)
    bx <- s$boxes
    expect_true(all(bx$x - bx$w / 2 >= 0 & bx$x + bx$w / 2 <= s$W & bx$y - bx$h / 2 >= 0 & bx$y + bx$h / 2 <= s$H), info = nm)
    expect_true(all(bx$y - bx$h / 2 > b[["y1"]]), info = nm)
  }
})

test_that("unresolved HSPC-start lines leave the bracket and never touch a candidate box", { specs <- get_specs()
  for (nm in names(specs)) {
    s <- specs[[nm]]
    expect_true(length(s$bus) > 0, info = nm)
    for (k in names(s$bus)) {
      p <- s$bus[[k]]
      expect_equal(unname(p[1, ]), unname(c(s$bracket[["xc"]], s$bracket[["y0"]])), info = paste(nm, k))
      for (i in seq_len(nrow(s$candidates))) expect_false(inside(p, s$candidates[i, ]), info = paste(nm, k))
    }
  }
  expect_setequal(names(specs$core$bus), rv$source_edge_map$source_edge_id[rv$source_edge_map$treatment == "unresolved" & rv$source_edge_map$view_panel == "rna_state"])
  lin <- unique(unlist(lapply(specs[names(RDR_LINEAGE_GROUPS)], function(s) names(s$bus))))
  expect_setequal(lin, HSPC_START)
})

test_that("drawn relations never pass through a box other than their own ends", { specs <- get_specs()
  for (nm in names(specs)) {
    s <- specs[[nm]]
    for (r in seq_along(s$routes)) {
      q <- s$routes[[r]]; mid <- .rdr_dense(q$pts, 0.02)
      others <- s$boxes[!s$boxes$id %in% c(q$from, q$to), , drop = FALSE]
      for (i in seq_len(nrow(others))) expect_false(inside(mid, others[i, ], pad = -0.01), info = paste(nm, q$id, "through", others$id[i]))
    }
    for (k in names(s$bus)) {
      mid <- .rdr_dense(s$bus[[k]], 0.02)
      others <- s$boxes[s$boxes$id != s$bus_target[[k]], , drop = FALSE]
      for (i in seq_len(nrow(others))) expect_false(inside(mid, others[i, ], pad = -0.01), info = paste(nm, k, "through", others$id[i]))
    }
  }
})

# ── 공개 산출물 ──────────────────────────────────────────────────────────────
test_that("public outputs contain no private paths or internal tags", {
  pub <- c(list.files(file.path(root, "data", "priors", "rna_dag_review", "proposal"), full.names = TRUE),
           list.files(file.path(root, "docs", "rna_dag_review"), pattern = "\\.(md|tsv)$", full.names = TRUE),
           list.files(file.path(root, "R", "rna_dag_review"), full.names = TRUE),
           list.files(file.path(root, "tests", "R", "rna_dag_review"), full.names = TRUE))
  expect_true(length(pub) >= 8L)
  for (f in pub) {
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    if (basename(f) == "test_rna_dag_review.R") x <- x[!grepl("^(PRIVATE|OVERCLAIM) <-", x)]
    expect_false(any(grepl(PRIVATE, x, ignore.case = TRUE, perl = TRUE)), info = basename(f))
  }
})
