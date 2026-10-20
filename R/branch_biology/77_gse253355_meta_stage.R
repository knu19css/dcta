# Report-derived stage peak markers on fixed GSE253355 clusters and R UMAP.
# User authorized Codex R implementation for this experiment on 2026-10-20.
# R75 scoring is reused; only five panels are replaced. The source report's
# transition signs and functional panel remain evidence, not classifier features.
# Run from the worktree: run-r dcta R/branch_biology/77_gse253355_meta_stage.R
# Optional: --out-dir <empty directory> --public-dir <empty directory>.

G77_STAGES <- c("HSC", "MPP", "CMP", "GMP", "MB", "PM", "MY", "MM", "BN", "SN")
G77_REMOVE <- c("HSC", "MPP", "CMP", "GMP", "Neutrophil")
G77_ROOT <- "/Users/edemhw/Honglab/CS/05_Neutrophil classification"
G77_BASE <- "neutrophil_atlas_rna_2026-10-19/gse253355_meta_stage_01"
G77_BRIDGE <- "results/claude_bridge/neutrophil_atlas_rna_2026-10-19"
G77_ARMS <- c("baseline_current", "meta_stage")

g77_main <- function(argv) {
  root <- getwd()
  a <- list(out = file.path(root, "results", G77_BASE), public = file.path(root, "docs/experiments", G77_BASE))
  if (length(argv) %% 2L) stop("Arguments require values")
  if (length(argv)) for (i in seq(1L, length(argv), 2L)) {
    key <- c("--out-dir" = "out", "--public-dir" = "public")[[argv[i]]]
    if (is.null(key)) stop("Unknown argument: ", argv[i])
    a[[key]] <- argv[i + 1L]
  }
  for (p in a) if (dir.exists(p) && length(list.files(p, all.files = TRUE, no.. = TRUE))) stop("Refusing to overwrite: ", p)
  dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) {
    line <- paste(format(Sys.time(), "%H:%M:%S"), ...)
    cat(line, "\n"); cat(line, "\n", file = file.path(a$out, "run.log"), append = TRUE)
  }
  src <- c(GSP_SOURCES, "77_gse253355_meta_stage.R")
  snap <- mdu_snapshot_sources(file.path(root, "R/branch_biology", src), file.path(a$out, "source_snapshot"))
  ld <- mdu_load_snapshot(snap, files = src)
  get("g77_build", ld$env, inherits = FALSE)(root, a$out, a$public, snap, log, argv)
}

# Same UMAP values and cell IDs are used for every panel; no inferred arrows.
g77_draw <- function(path, um, cluster, ranks, palette, subset = FALSE, fam = "") {
  keep <- if (subset) um$subset_mask else rep(TRUE, nrow(um$U))
  U <- um$U[keep, , drop = FALSE]; cl <- cluster[keep]; au <- as.character(um$author_label[keep])
  bounds <- mbvu_bounds(U); ac <- .gsu_cols(um$author_label, um$subset_labels)
  cnt <- table(cl); numbered <- names(cnt)[if (subset) cnt >= 0.01 * nrow(U) else cnt > 0]
  labels <- list(author = au)
  for (arm in G77_ARMS) labels[[arm]] <- ranks[[arm]]$top1_label[match(cl, ranks[[arm]]$cluster)]
  titles <- c("저자 주석 · 사후 대조", "기존 패널 · 30유형", "보고서 정점 패널 · 35유형")
  nf <- .mbvp_note_factory(basename(path)); note <- nf$note
  grDevices::png(path, width = 5400, height = 2600, res = 100, pointsize = 12)
  on.exit(grDevices::dev.off())
  graphics::layout(matrix(1:6, nrow = 2, byrow = TRUE), heights = c(3.2, 1.3))
  graphics::par(oma = c(5, 0, 6, 0), family = fam)
  for (j in seq_along(labels)) {
    lab <- labels[[j]]; col <- if (j == 1L) ac else palette
    stopifnot(!anyNA(col[lab]))
    graphics::par(mar = c(0.8, 0.7, 4.5, 0.7)); mbvu_open_panel(bounds)
    o <- .gsu_order(if (j == 1L) au else cl)
    graphics::points(U[o, 1], U[o, 2], pch = 16, cex = if (subset) 0.43 else 0.27, col = col[lab[o]])
    graphics::box(col = "grey60"); graphics::title(main = titles[j], cex.main = 2.2, line = 2)
    note(names(labels)[j], "title", titles[j])
    if (j > 1L) {
      rc <- atlas_cl_rep_cells(U, cl); rc <- rc[rc$cluster %in% numbered, , drop = FALSE]
      if (subset) rc$cluster <- vapply(rc$cluster, function(k) paste0(k, " · ", lab[match(k, cl)]), "")
      ly <- .acl_label_layout(rc, if (subset) 1.8 else 1.6, bounds, rc$n_cells)
      .acl_draw_labels(ly, ly$text, if (subset) 1.8 else 1.6)
      for (v in ly$text) note(names(labels)[j], "label", v)
    } else if (subset) {
      rc <- atlas_cl_rep_cells(U, au); ly <- .acl_label_layout(rc, 1.65, bounds, rc$n_cells)
      .acl_draw_labels(ly, ly$text, 1.65)
    }
  }
  for (j in seq_along(labels)) {
    lab <- labels[[j]]; col <- if (j == 1L) ac else palette
    tb <- sort(table(lab), decreasing = TRUE)
    items <- paste0(names(tb), " (", format(as.integer(tb), big.mark = ",", trim = TRUE), ")")
    graphics::par(mar = c(0, 1, 0, 1)); graphics::plot.new()
    graphics::legend("top", legend = items, col = col[names(tb)], pch = 16, bty = "n", pt.cex = 2.8,
                     cex = 1.6, ncol = if (length(items) > 24L) 3L else if (length(items) > 8L) 2L else 1L,
                     title = "라벨 (세포 수)", title.adj = 0)
    for (v in items) note(names(labels)[j], "legend", v)
  }
  title <- if (subset) "HSC부터 호중구·단핵구까지 · 같은 R UMAP에서 재라벨링" else "전체 골수 · 같은 R UMAP에서 재라벨링"
  graphics::mtext(title, 3, outer = TRUE, line = 3.3, cex = 2.8, font = 2)
  graphics::mtext(sprintf("GSE253355 · %s 세포 · 기존 26군집 유지 · 군집 평균의 최고 점수 후보", format(nrow(U), big.mark = ",")),
                  3, outer = TRUE, line = 0.9, cex = 1.6)
  foot <- c("단계 정점 마커 76개로 HSC/MPP/CMP/GMP/Neutrophil 패널만 교체. 다른 25유형 유지. 두 비교에 동일한 배경 유전자 사용.",
            "점수는 확률이 아니며 잠정 후보이다. 벌크 단계 정점은 전체 골수 유형 특이성을 보장하지 않는다. 저자 주석은 계산 후 비교용.")
  for (i in seq_along(foot)) graphics::mtext(foot[i], 1, outer = TRUE, line = 1 + (i - 1L) * 1.8, cex = 1.4)
  list(coords = U, text = nf$table(), numbered = numbered)
}

g77_heatmap <- function(path, scores, conditional, clusters, fam) {
  z <- scores[scores$arm == "meta_stage" & scores$type %in% G77_STAGES, ]
  M <- matrix(NA_real_, length(G77_STAGES), length(clusters), dimnames = list(G77_STAGES, clusters))
  for (i in seq_len(nrow(z))) if (z$cluster[i] %in% clusters) M[z$type[i], z$cluster[i]] <- z$score[i]
  grDevices::png(path, width = 1800, height = 1500, res = 130, pointsize = 12)
  on.exit(grDevices::dev.off()); graphics::par(mar = c(5, 5, 6, 2), oma = c(3, 0, 0, 0), family = fam)
  lim <- max(abs(M), na.rm = TRUE); pal <- grDevices::hcl.colors(101, "Blue-Red 3")
  graphics::plot(NA, xlim = c(0.5, ncol(M) + 0.5), ylim = c(0.5, nrow(M) + 0.5), xaxs = "i", yaxs = "i",
                 axes = FALSE, xlab = "군집 번호", ylab = "")
  for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M))) {
    v <- M[i, j]; ci <- round((v + lim) / (2 * lim) * 100) + 1L; y <- nrow(M) + 1L - i
    graphics::rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = if (is.na(v)) "grey90" else pal[ci], border = "white")
    win <- conditional$top1[match(clusters[j], conditional$cluster)] == G77_STAGES[i]
    graphics::text(j, y, paste0(sprintf("%.2f", v), if (isTRUE(win)) " *" else ""), cex = 1.15,
                   col = if (!is.na(v) && abs(v) > 0.65 * lim) "white" else "black")
  }
  graphics::axis(1, at = seq_along(clusters), labels = clusters, cex.axis = 1.2)
  graphics::axis(2, at = rev(seq_along(G77_STAGES)), labels = G77_STAGES, las = 2, cex.axis = 1.2)
  graphics::title(main = "10단계 안에서의 점수 비교", cex.main = 1.7, line = 3)
  graphics::mtext("* = 10단계 내 최고 후보 · 전체 35유형 중 최고 후보와 다를 수 있음", 3, line = 1, cex = 1)
  graphics::mtext("붉음 = 높은 점수 · 파랑 = 낮은 점수 · 단계별 재표준화 없음 · 계통 판정이나 정확도가 아님", 1, outer = TRUE, line = 1, cex = 0.95)
  M
}

g77_build <- function(root, out, public, snap, log, argv) {
  stamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ck <- .mbvq_checks(); files <- character(); roundtrip <- logical()
  wc <- function(x, name) {
    path <- file.path(out, name); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    roundtrip[[name]] <<- atlas_write_csv(x, path); files <<- unique(c(files, name)); invisible(path)
  }
  rdc <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  rdt <- function(p) utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
  src <- file.path(G77_ROOT, "04_results/10_meta_panel")
  reference <- c(report = file.path(G77_ROOT, "01_docs/05_neutrophil_stage_mrna_panel_report.md"),
                 peak = file.path(src, "10_stage_peak_markers.tsv"), transitions = file.path(src, "10_stage_transition_panel.tsv"),
                 curated = file.path(src, "10_curated_panel_grades.tsv"), consensus = file.path(src, "10_consensus_relative_level.tsv.gz"),
                 datasets = file.path(src, "10_dataset_stage_counts.tsv"))
  baseline <- file.path(root, G77_BRIDGE, "codex_meta_stage_baseline.json")
  pins <- jsonlite::fromJSON(baseline)
  baseline_hash <- function() vapply(names(pins), function(p) .nrv_sha(if (startsWith(p, "/")) p else file.path(root, p)), "")
  stopifnot(identical(unname(baseline_hash()), unname(unlist(pins))))
  input <- c(reference, markers = GSH_MARKERS, umap = file.path(root, GSH_UMAP),
             lossless = file.path(root, GSH_SCORES_DIR, "lossless.rds"), object = file.path(root, GSH_OBJECT),
             calc_manifest = file.path(root, GSH_SCORES_DIR, "calc_manifest.json"))
  h0 <- vapply(input, .nrv_sha, "")
  stopifnot(all(h0[names(GSH_PINNED)] == unlist(GSH_PINNED)))
  ck$add("pinned_inputs", TRUE, sprintf("%d protected files and raw counts / UMAP / clusters / markers match their hashes", length(pins)))
  peak <- rdt(reference[["peak"]]); trans <- rdt(reference[["transitions"]]); ds <- rdt(reference[["datasets"]])
  stopifnot(nrow(peak) == 76L, !anyDuplicated(paste(peak$stage, peak$gene)), setequal(peak$stage, G77_STAGES), all(peak$grade %in% c("A", "B")), nrow(trans) == 177L)
  stopifnot(!any(grepl("GSE253355", ds$dataset, fixed = TRUE)))
  cur <- gsh_panels(gsh_read_markers(GSH_MARKERS))$used
  cur <- cur[cur$model == "slide2", c("node_id", "display_label", "gene", "sign")]
  removed <- cur[cur$node_id %in% G77_REMOVE, ]; retained <- cur[!cur$node_id %in% G77_REMOVE, ]
  new <- data.frame(node_id = peak$stage, display_label = peak$stage, gene = peak$gene, sign = 1L)
  arms <- list(baseline_current = cur, meta_stage = rbind(retained, new))
  stopifnot(length(unique(cur$node_id)) == 30L, length(unique(retained$node_id)) == 25L, length(unique(arms$meta_stage$node_id)) == 35L)
  ck$add("panel_replacement_exact", TRUE, "Old five positive AND negative panels removed; 76 source rows inserted unchanged and positive-only; other 25 types retained exactly; target dataset not in reference dataset list")
  for (p in reference[names(reference) != "consensus"]) {
    name <- file.path("reference_input", basename(p)); q <- file.path(out, name); dir.create(dirname(q), showWarnings = FALSE)
    stopifnot(file.copy(p, q, overwrite = FALSE), identical(.nrv_sha(p), .nrv_sha(q))); files <- c(files, name)
  }
  cp <- file.path(out, "reference_input/current_markers.csv"); stopifnot(file.copy(GSH_MARKERS, cp, overwrite = FALSE)); files <- c(files, "reference_input/current_markers.csv")
  con <- gzfile(reference[["consensus"]], "rt"); reflevel <- rdt(con); close(con)
  stopifnot(all(peak$gene %in% reflevel$gene))
  wc(reflevel[match(peak$gene, reflevel$gene), c("gene", G77_STAGES)], "reference_peak_relative_levels.csv")
  rm(reflevel)
  um <- readRDS(input[["umap"]]); ls <- readRDS(input[["lossless"]]); cells <- ls$cells
  cluster <- as.character(ls$clusters); donor <- as.character(ls$donor); clusters <- as.character(0:25)
  obj <- gsx_open_rds(input[["object"]])$object; X <- gsx_counts(obj, "RNA")$counts; rm(obj); invisible(gc())
  stopifnot(identical(um$cells, cells), identical(rownames(um$U), cells), identical(colnames(X), cells),
            identical(as.character(um$new_cluster), cluster), setequal(unique(cluster), clusters), length(cells) == 82742L, sum(um$subset_mask) == 18892L)
  lib <- Matrix::colSums(X); stopifnot(identical(unname(as.numeric(lib)), ls$library))
  ck$add("same_cells_clusters_umap", TRUE, "82,742 cell IDs / order, 26 clusters, 18,892 plotting subset and all-gene libraries identical; no UMAP refit")
  log("Raw counts ready; aggregating log1p(CP10k)")
  agg <- gsh_aggregate(X, lib, cluster, donor, clusters); bins <- gsh_bins(agg$gene_mean)
  active <- unique(c(cur$gene, arms$meta_stage$gene)); measured <- intersect(active, names(bins))
  controls <- gsh_controls(bins, measured, active)
  scores <- do.call(rbind, lapply(G77_ARMS, function(m) data.frame(arm = m, gsh_type_scores(agg, arms[[m]], controls$controls))))
  ranks <- lapply(setNames(G77_ARMS, G77_ARMS), function(m) gsh_rank(scores[scores$arm == m, ]))
  donor_ranks <- lapply(setNames(G77_ARMS, G77_ARMS), function(m) gsh_rank(scores[scores$arm == m, ], "donor_equal_score"))
  conditional <- gsh_rank(scores[scores$arm == "meta_stage" & scores$type %in% G77_STAGES, ])
  peak$measured <- peak$gene %in% measured; peak$sign_used <- 1L
  wc(peak, "stage_marker_panel.csv"); wc(cur, "baseline_panel_used.csv"); wc(arms$meta_stage, "meta_panel_used.csv"); wc(removed, "old_panel_rows_removed.csv")
  wc(scores, "cluster_type_scores.csv"); wc(controls$table, "control_genes.csv")
  ss <- scores; ss$model <- ss$arm; wc(gsh_type_peaks(ss), "type_top_clusters.csv")
  for (m in G77_ARMS) { wc(ranks[[m]], paste0(m, "_ranks.csv")); wc(donor_ranks[[m]], paste0(m, "_donor_equal_ranks.csv")) }
  wc(conditional, "conditional_stage_ranks.csv")
  genes <- intersect(peak$gene, rownames(agg$C))
  expr <- do.call(rbind, lapply(genes, function(g) data.frame(gene = g, stage = peak$stage[match(g, peak$gene)], cluster = clusters,
    mean_log1p_cp10k = unname(agg$C[g, ]), detection_fraction = unname(agg$D[g, ]), background_mean = colMeans(agg$C[controls$controls[[g]], , drop = FALSE]))))
  wc(expr, "stage_marker_cluster_expression.csv")
  ver <- do.call(rbind, lapply(list(c("meta_stage", "SN"), c("meta_stage", "PM"), c("baseline_current", "HSC")), function(v) {
    p <- arms[[v[1]]]; p <- p[p$node_id == v[2], ]; pg <- intersect(p$gene[p$sign == 1], measured); ng <- intersect(p$gene[p$sign == -1], measured)
    x <- gsh_cell_level(X, lib, cluster, clusters, pg, ng, controls$controls)
    s <- scores[scores$arm == v[1] & scores$type == v[2], ]; s <- s[match(clusters, s$cluster), ]
    data.frame(arm = v[1], type = v[2], cluster = clusters, pos = s$pos_enrichment, pos_cell = x$pos_cell_mean, neg = s$neg_enrichment, neg_cell = x$neg_cell_mean)
  }))
  err <- max(abs(c(ver$pos - ver$pos_cell, ver$neg - ver$neg_cell)), na.rm = TRUE)
  wc(ver, "check_cell_level_recompute.csv"); ck$add("cell_level_recalculation", err < 1e-10, sprintf("SN, PM, HSC, all 26 clusters; max abs component error %.3g", err))
  rm(X); invisible(gc())
  ck$add("shared_controls", !any(unlist(controls$controls) %in% active) && all(vapply(names(controls$controls), function(g) all(bins[controls$controls[[g]]] == bins[g]), TRUE)),
         sprintf("%d measured marker genes; both arms share same-bin controls excluding union of %d active genes", length(measured), length(active)))
  ck$add("score_universe", nrow(scores) == 1690L && !anyDuplicated(paste(scores$arm, scores$type, scores$cluster)) && all(is.finite(scores$score[scores$n_pos_measured > 0])), "26 x (30 + 35), including all low scores; no forced stage assignment")
  calc_files <- files
  calc_hash <- setNames(vapply(file.path(out, calc_files), .nrv_sha, ""), calc_files)
  jsonlite::write_json(list(written_at = stamp(), phase = "Scores and ranks fixed before author-label comparison", files = as.list(calc_hash)), file.path(out, "calc_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  files <- c(files, "calc_manifest.json")
  saveRDS(list(scores = scores, ranks = ranks, donor_ranks = donor_ranks, conditional = conditional, agg = agg, controls = controls, arms = arms), file.path(out, "calculation.rds"))
  log("All 1690 scores fixed; post hoc comparisons and figures")
  # Author labels and earlier output labels are used only after the calculation.
  au <- as.character(um$author_label); comp <- gsm_author_composition(cluster, au, clusters)
  side <- data.frame(cluster = clusters, n_cells = as.integer(agg$n_cells[clusters]), n_subset_cells = as.integer(table(factor(cluster[um$subset_mask], levels = clusters))))
  side$author_top_label <- vapply(clusters, function(k) comp$author_label[comp$cluster == k][1], "")
  side$author_top_fraction <- vapply(clusters, function(k) comp$frac_of_cluster[comp$cluster == k][1], 0)
  for (m in G77_ARMS) {
    r <- ranks[[m]][match(clusters, ranks[[m]]$cluster), ]; d <- donor_ranks[[m]][match(clusters, donor_ranks[[m]]$cluster), ]
    for (k in c("status", "top1", "top1_label", "top1_score", "top2_label", "top2_score", "gap")) side[[paste0(m, "_", k)]] <- r[[k]]
    side[[paste0(m, "_donor_equal_top1")]] <- d$top1_label
  }
  ci <- conditional[match(clusters, conditional$cluster), ]
  side$conditional_stage_candidate <- ci$top1_label; side$conditional_stage_gap <- ci$gap
  side$conditional_disagrees_overall <- ci$top1 != side$meta_stage_top1
  side$note <- "Provisional candidates; conditional stage does not establish lineage; author labels post hoc"
  wc(side, "cluster_side_by_side.csv"); wc(comp, "posthoc_author_composition.csv")
  wc(side[, c("cluster", "conditional_stage_candidate", "conditional_stage_gap", "meta_stage_top1_label", "conditional_disagrees_overall", "note")], "cluster_stage_candidates.csv")
  r75 <- rdc(file.path(root, GSP_R75, "cluster_provisional_annotation.csv"))
  r76 <- rdc(file.path(root, "results/neutrophil_atlas_rna_2026-10-19/gse253355_manual_panel_01/cluster_side_by_side.csv"))
  shift <- data.frame(cluster = clusters, r75 = r75$slide2_top1[match(clusters, r75$cluster)], r76_current = r76$current_top1[match(clusters, r76$cluster)], this_baseline = side$baseline_current_top1_label)
  shift$changed_from_r75 <- shift$r75 != shift$this_baseline; shift$changed_from_r76_current <- shift$r76_current != shift$this_baseline
  shift$note <- "Same current panel and scoring; changed common background genes can alter baseline winners"
  wc(shift, "baseline_background_shifts.csv")
  ck$add("calculation_unchanged_posthoc", identical(unname(vapply(file.path(out, calc_files), .nrv_sha, "")), unname(calc_hash)), "All fixed score/rank/input evidence files unchanged after author labels read")
  fam <- if (mbvr_hangul_font()) MBVR_TITLE_FAMILY else ""
  keys <- unique(unlist(lapply(ranks, function(r) r$top1_label)))
  palette <- setNames(grDevices::hcl.colors(length(keys), "Dark 3"), keys)
  if ("undetermined" %in% keys) palette["undetermined"] <- GSH_UNDET_COL
  dir.create(file.path(out, "figures"), showWarnings = FALSE)
  f1 <- "figures/fig1_umap_full_author_baseline_meta.png"; f2 <- "figures/fig2_umap_lineage_author_baseline_meta.png"; f3 <- "figures/fig3_stage_scores.png"
  r1 <- g77_draw(file.path(out, f1), um, cluster, ranks, palette, FALSE, fam)
  r2 <- g77_draw(file.path(out, f2), um, cluster, ranks, palette, TRUE, fam)
  major <- c("8", "17", "16", "7", "12", "15", "6")
  heat <- g77_heatmap(file.path(out, f3), scores, conditional, major, fam)
  files <- c(files, f1, f2, f3); wc(rbind(r1$text, r2$text), "figure_text.csv")
  wc(data.frame(label = names(palette), colour = unname(palette)), "figure_colours.csv")
  same_xy <- function(x, y) identical(as.vector(x), as.vector(y)) && identical(dimnames(x), dimnames(y)) && identical(dim(x), dim(y))
  ck$add("figures_use_original_coordinates", same_xy(r1$coords, um$U) && same_xy(r2$coords, um$U[um$subset_mask, , drop = FALSE]), "Original R UMAP for all cells and exact author-defined plotting subset; no recomputation / trajectory arrows")
  wc(data.frame(input = names(input), path = unname(input), sha256 = unname(h0)), "input_sha256.csv")
  wc(data.frame(item = c("R_executable", "R_version", "libPaths", "direct_packages", "seed", "command", "written_at", "implementation"), value = c(file.path(R.home("bin"), "R"), R.version.string,
    paste(.libPaths(), collapse = ";"), paste(vapply(c("Matrix", "digest", "jsonlite"), function(p) paste0(p, "=", utils::packageVersion(p)), ""), collapse = ";"),
    GSH_SEED, paste(c("run-r dcta R/branch_biology/77_gse253355_meta_stage.R", argv), collapse = " "), stamp(), "Codex; user authorized one-task R implementation exception")), "environment.csv")
  wc(snap, "source_snapshot.csv")
  ck$add("protected_sources_unchanged", identical(unname(baseline_hash()), unname(unlist(pins))) && identical(unname(vapply(input, .nrv_sha, "")), unname(h0)), sprintf("All %d prior/source files and %d inputs unchanged", length(pins), length(input)))
  ck$add("csv_roundtrip", all(roundtrip), sprintf("%d aggregate tables re-read successfully", length(roundtrip)))
  checks <- ck$table(); wc(checks, "checks.csv"); stopifnot(all(checks$pass))
  selected <- side[match(major, side$cluster), ]
  sm <- c("# 보고서 mRNA 정점 패널 재라벨링", "", "- 기존 82,742세포·26군집·R UMAP 유지. HSC→SN 76개 정점 마커만 교체, 다른 25유형 유지.",
    "- 점수 계산은 R75와 동일. 30유형 기존판과 35유형 새판에 공통 배경을 사용하며 기본 패널은 바꾸지 않음.", "", "|군집|저자 최다 주석|기존 후보|새 전체 후보|10단계 내 후보|", "|---|---|---|---|---|")
  sm <- c(sm, vapply(seq_len(nrow(selected)), function(i) paste0("|", paste(selected[i, c("cluster", "author_top_label", "baseline_current_top1_label", "meta_stage_top1_label", "conditional_stage_candidate")], collapse = "|"), "|"), ""),
    sprintf("- 패널 측정: %d/76 유전자. 검사 %d/%d 통과. 배경만 달라 기존판 후보가 R75 대비 %d/26, R76 대비 %d/26 변경.", sum(peak$measured), sum(checks$pass), nrow(checks), sum(shift$changed_from_r75), sum(shift$changed_from_r76_current)),
    "- 10단계 내 후보는 해당 계통이라는 판정이 아님. 다른 계통과의 경쟁 결과를 함께 볼 것.",
    "- 한계: 벌크 분획 정점은 전체 골수 특이성과 다름; 단계 경계 ±1, GMP→MB 추론. 참조 구축 자료와의 ρ=0.91은 독립 검증 아님.",
    "- 이행 증감 177행·기능69개는 근거 사본만 보존. 이번 점수에는 섞지 않음. 형태학 정답이 없어 정확도는 확인 못함.")
  stopifnot(length(sm) <= 20L); writeLines(sm, file.path(out, "summary_ko.md")); files <- c(files, "summary_ko.md")
  # Only aggregate results and small source-panel evidence are published.
  pubfiles <- setdiff(files, "source_snapshot.csv")
  stopifnot(!any(grepl("\\.rds$|\\.gz$", pubfiles)))
  for (d in unique(dirname(file.path(public, pubfiles)))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  stopifnot(all(file.copy(file.path(out, pubfiles), file.path(public, pubfiles), overwrite = FALSE)))
  jsonlite::write_json(list(source_output = normalizePath(out), files = as.list(setNames(vapply(file.path(public, pubfiles), .nrv_sha, ""), pubfiles))), file.path(public, "public_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  log("Complete:", nrow(checks), "checks passed; public files", length(pubfiles))
  invisible(side)
}

if (!interactive() && sys.nframe() == 0L) {
  bb <- file.path(getwd(), "R/branch_biology")
  for (f in c("17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R", "76_gse253355_manual_panel.R")) source(file.path(bb, f))
  tryCatch(g77_main(commandArgs(TRUE)), error = function(e) { message(conditionMessage(e)); quit(status = 4L) })
}
