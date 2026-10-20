# User follow-up: restrict candidates to HSC -> neutrophil, 10 stages.
# Candidate set restriction ONLY: R77 expression scores / backgrounds stay fixed.
# Author labels select the displayed lineage (no monocytes); they do not select
# markers, change scores, define cutoffs, or establish lineage ground truth.
# run-r dcta R/branch_biology/78_gse253355_lineage_only.R

G78_BASE <- "neutrophil_atlas_rna_2026-10-19/gse253355_lineage_only_01"
G78_PLOT_LABELS <- c("HSC", "MPP", "Cycling HSPC", "GMP", "Early Myeloid Progenitor", "Late Myeloid", "Neutrophil")
G78_STAGES <- c("HSC", "MPP", "CMP", "GMP", "MB", "PM", "MY", "MM", "BN", "SN")

g78_draw <- function(path, U, author, cluster, ranks, fam) {
  bounds <- mbvu_bounds(U); ac <- .gsu_cols(author, G78_PLOT_LABELS)
  stagecols <- setNames(grDevices::hcl.colors(10, "Spectral", rev = TRUE), G78_STAGES)
  stagecols <- c(stagecols, undetermined = "grey80")
  labels <- ranks$top1_label[match(cluster, ranks$cluster)]
  stopifnot(!anyNA(labels), all(labels %in% names(stagecols)))
  counts <- table(cluster); numbered <- names(counts)[counts >= 0.01 * nrow(U)]
  grDevices::png(path, width = 2800, height = 1900, res = 110, pointsize = 12)
  on.exit(grDevices::dev.off())
  graphics::layout(matrix(c(1, 2, 3, 4), 2, byrow = TRUE), heights = c(4, 0.9))
  graphics::par(oma = c(4.3, 0, 4.5, 0), family = fam)
  for (j in 1:2) {
    x <- if (j == 1L) author else labels; cols <- if (j == 1L) ac else stagecols
    graphics::par(mar = c(0.4, 0.4, 3.1, 0.4)); mbvu_open_panel(bounds)
    o <- .gsu_order(if (j == 1L) author else cluster)
    graphics::points(U[o, 1], U[o, 2], pch = 16, cex = 0.42, col = cols[x[o]])
    graphics::box(col = "grey65")
    if (j == 1L) {
      rp <- atlas_cl_rep_cells(U, author)
    } else {
      rp <- atlas_cl_rep_cells(U, cluster); rp <- rp[rp$cluster %in% numbered, , drop = FALSE]
      rp$cluster <- vapply(rp$cluster, function(k) paste0(k, " · ", labels[match(k, cluster)]), "")
    }
    ly <- .acl_label_layout(rp, 1.25, bounds, rp$n_cells); .acl_draw_labels(ly, ly$text, 1.25)
    graphics::title(main = if (j == 1L) "저자 주석 · 참고용" else "10단계 중 최고 점수 후보", cex.main = 1.6, line = 1)
  }
  for (j in 1:2) {
    x <- if (j == 1L) author else labels; cols <- if (j == 1L) ac else stagecols
    tb <- table(x); order <- if (j == 1L) G78_PLOT_LABELS else c(G78_STAGES, "undetermined")
    present <- order[order %in% names(tb)]
    items <- sprintf("%s (%s)", present, format(as.integer(tb[present]), big.mark = ",", trim = TRUE))
    graphics::par(mar = c(0, 0.5, 0, 0.5)); graphics::plot.new()
    graphics::legend("top", legend = items, col = cols[present], pch = 16, pt.cex = 1.7, cex = 1.05, ncol = 2, bty = "n", title = "표시된 세포 수")
  }
  graphics::mtext("HSC부터 호중구까지 · 분류 후보를 10단계로 제한", 3, outer = TRUE, line = 2.3, cex = 1.9, font = 2)
  graphics::mtext(sprintf("%s 세포 표시 · 단핵구 제외 · 기존 R UMAP과 군집 유지", format(nrow(U), big.mark = ",")), 3, outer = TRUE, line = 0.5, cex = 1.1)
  graphics::mtext("후보: HSC / MPP / CMP / GMP / MB / PM / MY / MM / BN / SN", 1, outer = TRUE, line = 0.9, cex = 1.1)
  graphics::mtext("점수는 원래 군집 전체에서 계산한 R77 값을 유지. 저자 라벨은 표시 범위에만 사용. 잠정 후보이며 단계 확정 아님.", 1, outer = TRUE, line = 2.5, cex = 0.95)
  list(coords = U, numbered = numbered, colours = data.frame(stage = names(stagecols), colour = unname(stagecols)))
}

g78_main <- function() {
  root <- getwd(); out <- file.path(root, "results", G78_BASE); pub <- file.path(root, "docs/experiments", G78_BASE)
  for (p in c(out, pub)) if (dir.exists(p) && length(list.files(p, all.files = TRUE, no.. = TRUE))) stop("Refusing overwrite: ", p)
  dir.create(out, recursive = TRUE); dir.create(pub, recursive = TRUE)
  src <- c(GSP_SOURCES, "77_gse253355_meta_stage.R", "78_gse253355_lineage_only.R")
  snap <- mdu_snapshot_sources(file.path(root, "R/branch_biology", src), file.path(out, "source_snapshot"))
  ld <- mdu_load_snapshot(snap, files = src)
  get("g78_build", ld$env, inherits = FALSE)(root, out, pub)
}

g78_build <- function(root, out, pub) {
  prior <- file.path(root, "docs/experiments/neutrophil_atlas_rna_2026-10-19/gse253355_meta_stage_01")
  pin <- jsonlite::fromJSON(file.path(root, "results/claude_bridge/neutrophil_atlas_rna_2026-10-19/codex_lineage_only_baseline.json"))
  phash <- function() vapply(names(pin), function(p) .nrv_sha(if (startsWith(p, "/")) p else file.path(root, p)), "")
  stopifnot(identical(unname(phash()), unname(unlist(pin))))
  rd <- function(f) utils::read.csv(file.path(prior, f), stringsAsFactors = FALSE, check.names = FALSE)
  old <- rd("cluster_type_scores.csv"); sc <- old[old$arm == "meta_stage" & old$type %in% G78_STAGES, ]
  sc$cluster <- as.character(sc$cluster)  # IDs (including "0"), never numeric table positions.
  stopifnot(nrow(sc) == 260L, !anyDuplicated(paste(sc$cluster, sc$type)), setequal(sc$type, G78_STAGES))
  rank <- gsh_rank(sc); donor_rank <- gsh_rank(sc, "donor_equal_score"); fixed <- rd("conditional_stage_ranks.csv")
  stopifnot(identical(as.character(rank$cluster), as.character(fixed$cluster)), identical(rank$top1_label, fixed$top1_label), max(abs(rank$top1_score - fixed$top1_score)) < 1e-12)
  # Independent maxima and next-score checks; no imposed stage order or endpoint.
  for (i in seq_len(nrow(rank))) {
    v <- sc$score[as.character(sc$cluster) == rank$cluster[i]]; v <- sort(v[is.finite(v)], decreasing = TRUE)
    stopifnot(abs(rank$top1_score[i] - max(0, v[1])) < 1e-12)
    if (sum(v >= v[1] - 1e-12) == 1L && v[1] > 0) stopifnot(abs(rank$gap[i] - (v[1] - max(0, v[2]))) < 1e-12)
  }
  umpath <- file.path(root, GSH_UMAP); stopifnot(.nrv_sha(umpath) == GSH_PINNED$umap)
  um <- readRDS(umpath); cl <- as.character(um$new_cluster); author <- as.character(um$author_label)
  mask <- author %in% G78_PLOT_LABELS; stopifnot(!any(author[mask] == "Monocyte"), all(mask <= um$subset_mask))
  counts <- table(cl[mask]); shown <- rank$cluster[rank$cluster %in% names(counts)]
  side <- rank[match(shown, rank$cluster), ]; side$original_cluster_cells_used_for_score <- sc$n_cells[match(shown, sc$cluster)]
  side$lineage_cells_shown <- as.integer(counts[shown]); side$donor_equal_candidate <- donor_rank$top1_label[match(shown, donor_rank$cluster)]
  stopifnot(!anyNA(side$lineage_cells_shown), sum(side$lineage_cells_shown) == sum(mask),
            all(side$lineage_cells_shown == vapply(shown, function(k) sum(cl[mask] == k), 0L)))
  side$donor_equal_gap <- donor_rank$gap[match(shown, donor_rank$cluster)]
  side$major_display_cluster <- side$lineage_cells_shown >= 0.01 * sum(mask)
  side$author_top_in_display <- vapply(shown, function(k) names(sort(table(author[mask & cl == k]), decreasing = TRUE))[1], "")
  side$note <- "Candidate restriction only; scores from full original clusters, author labels restrict display only; provisional within-lineage candidates"
  wc <- function(x, f) { stopifnot(atlas_write_csv(x, file.path(out, f))) }
  wc(side, "cluster_lineage_candidates.csv")
  wc(sc[sc$cluster %in% shown, ], "lineage_cluster_stage_scores.csv")
  wc(data.frame(author_label = G78_PLOT_LABELS, n_cells = as.integer(table(factor(author[mask], levels = G78_PLOT_LABELS)))), "display_cell_counts.csv")
  panel <- rd("stage_marker_panel.csv"); wc(panel, "stage_marker_panel.csv")
  types <- data.frame(stage = G78_STAGES, n_major_clusters_top1 = vapply(G78_STAGES, function(s) sum(side$top1 == s & side$major_display_cluster), 0L))
  types$note <- "No highest-scoring major cluster does not establish absence of this stage"; wc(types, "stage_assignments_summary.csv")
  saveRDS(list(cells = um$cells[mask], U = um$U[mask, , drop = FALSE], cluster = cl[mask], candidate = side$top1_label[match(cl[mask], side$cluster)]), file.path(out, "displayed_cells.rds"))
  fam <- if (mbvr_hangul_font()) MBVR_TITLE_FAMILY else ""
  fig <- "fig1_hsc_neutrophil_only.png"
  view <- g78_draw(file.path(out, fig), um$U[mask, , drop = FALSE], author[mask], cl[mask], rank, fam)
  stopifnot(identical(view$coords, um$U[mask, , drop = FALSE]))
  wc(view$colours, "stage_colours.csv")
  major <- side[side$major_display_cluster, ]; major <- major[match(c("8", "17", "16", "7", "12", "15", "6"), major$cluster), ]
  sm <- c("# HSC→호중구 10단계만 적용", "", "- 분류 후보만 HSC/MPP/CMP/GMP/MB/PM/MY/MM/BN/SN으로 제한. R77 군집 점수·배경 유전자 그대로 유지.",
    sprintf("- 그림: 저자 HSC→Neutrophil 7개 라벨 %s세포. Monocyte 제외. 원래 UMAP·군집 유지, 군집을 다시 나누거나 발현 평균을 재추정하지 않음.", format(sum(mask), big.mark = ",")),
    "", "|군집|비교용 저자 최다 라벨|최고 후보|차점|점수 차이|", "|---|---|---|---|---|")
  sm <- c(sm, vapply(seq_len(nrow(major)), function(i) sprintf("|%s|%s|%s|%s|%.3f|", major$cluster[i], major$author_top_in_display[i], major$top1_label[i], major$top2_label[i], major$gap[i]), ""),
    "- 8번 HSC, 12·15번 MM, 6번 SN. 17번 MPP와 7번 GMP가 SN으로 남아 후보 제한만으로 초기 단계 혼동이 해결되지는 않음.",
    "- 점수 차이는 확률·신뢰도가 아님. 16번 MY 후보는 증식 신호와 단계 신호를 구별해야 함.",
    "- 이는 계통 안의 조건부 분류이며 계통 자체를 판별하거나 형태학적 성숙 단계를 검증한 결과가 아님.")
  writeLines(sm, file.path(out, "summary_ko.md"))
  stopifnot(length(sm) <= 20L, identical(unname(phash()), unname(unlist(pin))))
  verification <- list(score_rows_checked = nrow(sc), original_clusters = length(unique(sc$cluster)), displayed_cells = sum(mask), displayed_clusters = length(shown),
    nonlineage_candidates_removed = TRUE, frozen_R77_scores_unchanged = TRUE, ranks_equal_prior_conditional = TRUE, independent_maxima_and_gaps_checked = TRUE,
    coordinates_identical = TRUE, monocytes_displayed = sum(author[mask] == "Monocyte"), protected_hashes_unchanged = length(pin),
    selection_note = "Author labels select displayed cells only. Entire original clusters still contribute to scores; no lineage-only expression refit.")
  jsonlite::write_json(verification, file.path(out, "verification.json"), auto_unbox = TRUE, pretty = TRUE)
  wc(data.frame(item = c("R", "R_executable", "libPaths", "packages", "score_source", "umap_sha256", "command"), value = c(R.version.string,
    file.path(R.home("bin"), "R"), paste(.libPaths(), collapse = ";"), paste(vapply(c("jsonlite", "digest"), function(p) paste0(p, "=", utils::packageVersion(p)), ""), collapse = ";"),
    .nrv_sha(file.path(prior, "cluster_type_scores.csv")), .nrv_sha(umpath), "run-r dcta R/branch_biology/78_gse253355_lineage_only.R")), "environment.csv")
  pubfiles <- list.files(out, pattern = "\\.(csv|md|json|png)$")
  stopifnot(all(file.copy(file.path(out, pubfiles), file.path(pub, pubfiles), overwrite = FALSE)))
  jsonlite::write_json(list(files = as.list(setNames(vapply(file.path(pub, pubfiles), .nrv_sha, ""), pubfiles))), file.path(pub, "public_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("Completed:", sum(mask), "displayed cells,", length(shown), "displayed clusters;", nrow(sc), "fixed scores\n")
}

if (!interactive() && sys.nframe() == 0L) {
  bb <- file.path(getwd(), "R/branch_biology")
  for (f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "76_gse253355_manual_panel.R")) source(file.path(bb, f))
  tryCatch(g78_main(), error = function(e) { message(conditionMessage(e)); quit(status = 4L) })
}
