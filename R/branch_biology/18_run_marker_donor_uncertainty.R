# =============================================================================
# 18_run_marker_donor_uncertainty.R — donor 기반 marker 불확실성 실행기 (스냅샷 실행, provenance 고정)
# =============================================================================
#
# 실행 (worktree 에서, dcta 환경만)
#   /Users/edemhw/Documents/Codex/2026-10-07/x20/outputs/agent-workflow/scripts/run-r dcta \
#     R/branch_biology/18_run_marker_donor_uncertainty.R \
#     --data-root /Users/edemhw/Honglab/CS/04_Claude/14_dcta \
#     --outdir results/marker_donor_uncertainty_2026-10-11/attempt_01 [--replicates 2000] [--seed 20261011]
#
# 순서
#   0. 출력 폴더가 비어 있지 않으면 거절(exit 6).
#   1. [bootstrap = mdu_main] 05 / 17 / 18 / 19 를 outdir/source_snapshot 에 복사·hash 고정(하나라도 없으면 실패).
#      사본 hash 를 다시 확인한 뒤 05 / 17 / 18 사본을 새 환경에 적재하고, 그 환경의 mdu_run_body 를 실행한다.
#   2. 적합 단계 (저자 라벨 없음): 원 RNA 를 chunk 로 읽어 (군집, donor) 합 → pooled 를 저장 기준과 대조
#      → donor 요약, 재표본, 부호 검정, 하나 빼기, 유전자 근거 → 표 저장 → fit_manifest.json (hash 고정)
#   3. 사후 단계: 저자 라벨·평가 대응표·기존 판정을 읽어 대조표를 만든다. 적합 표 hash 가 그대로인지 다시 확인한다.
# exit code: 0 정상, 2 인자, 4 입력 불일치/검사 실패, 6 출력 폴더 비어 있지 않음
# =============================================================================

mdu_find_root <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) {
    p <- normalizePath(file.path(dirname(sub("^--file=", "", a[1])), "..", ".."), mustWork = FALSE)
    if (dir.exists(file.path(p, "R", "branch_biology"))) return(p)
  }
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) { if (dir.exists(file.path(p, "R", "branch_biology"))) return(p); p <- dirname(p) }
  stop("프로젝트 루트를 찾지 못했다.")
}

MDU_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R",
                 "19_view_marker_donor_uncertainty.R")                                                 # 모두 필수
MDU_RUN_SOURCES <- c("05_cell_trajectory_io.R", "17_marker_donor_uncertainty.R", "18_run_marker_donor_uncertainty.R")
MDU_H5AD <- "data/external/hay_hca_source.h5ad"
MDU_MEMBERSHIP <- "results/mgp_context_2026-09-27/attempt_01/BoneMarrowMap/structure/membership.tsv.gz"
MDU_MODULE_DIR <- "results/mgp_context_2026-09-27/attempt_01/BoneMarrowMap/module_features/descriptive"
MDU_REFERENCE <- "results/partial_branch_annotation_2026-09-30/attempt_02/BoneMarrowMap/fit/res4.0_seed0/descriptive/unit_state_evidence.tsv.gz"
MDU_READOUT <- "results/simple_joint_annotation_2026-10-01/attempt_03/BoneMarrowMap/unit_readouts.tsv.gz"
MDU_MAPPING <- "results/path_context_2026-09-26/attempt_01/BoneMarrowMap/evaluation_mapping_used.tsv"
MDU_PANEL <- "results/projection_compact_2026-09-25/GSE181989/v0_graph/panel_adapted.tsv"   # 빠진 유전자 이름 기록용 (읽기 전용)
MDU_PARTITION <- "unit_res4.0_seed0"


# ── 스냅샷 (15 와 같은 방식) ─────────────────────────────────────────────────
mdu_sha <- function(p) digest::digest(p, algo = "sha256", file = TRUE)

mdu_snapshot_sources <- function(files, snap_dir) {
  miss <- files[!file.exists(files)]
  if (length(miss)) stop("필수 소스 파일 없음: ", paste(basename(miss), collapse = ", "))
  dir.create(snap_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(snap_dir, basename(files))
  if (!all(file.copy(files, dest, overwrite = FALSE))) stop("소스 스냅샷 복사 실패")
  sha <- vapply(files, mdu_sha, character(1)); sha2 <- vapply(dest, mdu_sha, character(1))
  if (!identical(unname(sha), unname(sha2))) stop("스냅샷 hash 가 원본과 다르다")
  data.frame(original_path = files, snapshot_path = dest, sha256 = unname(sha), stringsAsFactors = FALSE)
}

mdu_load_snapshot <- function(snap, files = MDU_RUN_SOURCES) {
  ix <- match(files, basename(snap$snapshot_path))
  if (anyNA(ix)) stop("스냅샷에 없는 소스: ", paste(files[is.na(ix)], collapse = ", "))
  loaded <- snap[ix, c("snapshot_path", "sha256")]
  now <- vapply(loaded$snapshot_path, mdu_sha, character(1))
  if (!identical(unname(now), loaded$sha256)) stop("스냅샷 사본 hash 가 고정값과 다르다")
  env <- new.env(parent = globalenv())
  for (p in loaded$snapshot_path) sys.source(p, envir = env)
  rownames(loaded) <- NULL
  list(env = env, loaded = loaded)
}

mdu_compare_sources <- function(snap) {
  now <- vapply(snap$original_path, function(f) if (file.exists(f)) mdu_sha(f) else NA_character_, character(1))
  data.frame(original_path = snap$original_path, sha256_at_start = snap$sha256, sha256_at_end = unname(now),
             same = !is.na(now) & unname(now) == snap$sha256, stringsAsFactors = FALSE)
}

#' R 에서 다시 계산할 때 쓰는 무손실 사본(RDS, 이진 double 그대로).
#' R 의 십진 파서(as.numeric/read.delim)는 17 자리 문자열을 1 ulp 틀리게 읽을 수 있어서 TSV 만으로는 R 왕복이 정확하지 않다.
mdu_save_lossless <- function(tables, path) saveRDS(tables, path, version = 3)

#' 표 저장. double 은 17 유효숫자(%.17g)로 쓴다 (Codex 13: 기본 15 자리 반올림 금지).
#' 정확히 반올림하는 파서에는 무손실이다. R 에서는 1 ulp 이내이므로, R 재검산은 mdu_save_lossless 사본을 쓴다.
mdu_write <- function(df, path) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_along(df)) if (is.double(df[[j]])) df[[j]] <- ifelse(is.na(df[[j]]), NA_character_, sprintf("%.17g", df[[j]]))
  if (grepl("\\.gz$", path)) { con <- gzfile(path, "w"); on.exit(close(con)) } else con <- path
  utils::write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}


#' bootstrap: 인자 해석 → 스냅샷 → 사본 적재 → 사본의 mdu_run_body 실행. root / after_snapshot 은 테스트용.
mdu_main <- function(argv, root = NULL, after_snapshot = NULL) {
  a <- list(outdir = NULL, data_root = NULL, replicates = 2000L, seed = 20261011L, chunk = 5000L)
  i <- 1L
  while (i <= length(argv)) {
    v <- if (i < length(argv)) argv[i + 1L] else NA_character_
    switch(argv[i], "--outdir" = { a$outdir <- v }, "--data-root" = { a$data_root <- v },
           "--replicates" = { a$replicates <- as.integer(v) }, "--seed" = { a$seed <- as.integer(v) }, "--chunk" = { a$chunk <- as.integer(v) },
           { cat("알 수 없는 인자:", argv[i], "\n"); return(2L) })
    i <- i + 2L
  }
  if (is.null(a$outdir)) { cat("--outdir 가 필요하다.\n"); return(2L) }
  if (dir.exists(a$outdir) && length(list.files(a$outdir, all.files = TRUE, no.. = TRUE))) { cat("출력 폴더가 비어 있지 않다:", a$outdir, "\n"); return(6L) }
  if (is.null(a$data_root)) { cat("--data-root 가 필요하다.\n"); return(2L) }
  if (is.na(a$replicates) || a$replicates < 1L) { cat("--replicates 는 양의 정수\n"); return(2L) }
  t0 <- Sys.time(); out <- a$outdir; dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), ...); cat(line, "\n"); cat(line, "\n", file = file.path(out, "run.log"), append = TRUE) }
  if (is.null(root)) root <- mdu_find_root()
  tryCatch({
    snap <- mdu_snapshot_sources(file.path(root, "R", "branch_biology", MDU_SOURCES), file.path(out, "source_snapshot"))
    log("source snapshot fixed:", paste(basename(snap$snapshot_path), substr(snap$sha256, 1, 12), collapse = "; "))
    if (is.function(after_snapshot)) after_snapshot()
    ld <- mdu_load_snapshot(snap)
    log("loaded from snapshot:", paste(basename(ld$loaded$snapshot_path), collapse = ", "))
    body <- get("mdu_run_body", envir = ld$env, inherits = FALSE)
    body(a = a, argv = argv, out = out, root = root, snap = snap, loaded = ld$loaded, t0 = t0, log = log)
    0L
  }, error = function(e) { log("ERROR", conditionMessage(e)); 4L })
}


#' 계산 본체 (스냅샷 사본에서 적재되어 실행된다)
mdu_run_body <- function(a, argv, out, root, snap, loaded, t0, log) {
  P <- function(x) file.path(a$data_root, x)
  inputs <- c(h5ad = P(MDU_H5AD), membership = P(MDU_MEMBERSHIP), state_table = P(file.path(MDU_MODULE_DIR, "state_table.tsv")),
              controls = P(file.path(MDU_MODULE_DIR, "controls.tsv.gz")), selection = P(file.path(MDU_MODULE_DIR, "selection.json")),
              reference = P(MDU_REFERENCE), readout = P(MDU_READOUT), mapping = P(MDU_MAPPING), panel = P(MDU_PANEL))
  miss <- inputs[!file.exists(inputs)]
  if (length(miss)) stop("입력 파일 없음: ", paste(names(miss), collapse = ", "))
  in_sha <- vapply(inputs, mdu_sha, character(1)); log("inputs hashed (", length(in_sha), ")")

  # ---- 1. 세포 메타 (저자 라벨 없음)
  f <- hdf5r::H5File$new(inputs[["h5ad"]], "r")
  idx_key <- hdf5r::h5attr(f[["obs"]], "_index")
  cell_id <- as.character(f[[paste0("obs/", idx_key)]]$read())
  donor <- ct_h5_obs(f, "donor_id"); study <- ct_h5_obs(f, "Study")
  gs <- f[["raw/var/gene_symbols"]]
  symbols <- if (inherits(gs, "H5Group")) { cats <- as.character(gs[["categories"]]$read()); cats[as.integer(gs[["codes"]]$read()) + 1L] } else as.character(gs$read())
  f$close_all()
  if (anyNA(symbols)) stop("raw/var/gene_symbols 에 결측이 있다.")
  memb <- utils::read.delim(gzfile(inputs[["membership"]]), stringsAsFactors = FALSE, colClasses = c(cell = "character"))
  mdu_check_ids(memb$cell, cell_id)
  cl <- as.character(memb[[MDU_PARTITION]])
  fitcells <- data.frame(cell = cell_id, cluster = cl, donor = donor, study = study, stringsAsFactors = FALSE)
  mdu_check_fit_columns(names(fitcells))
  if (anyNA(fitcells$cluster)) stop("군집이 없는 세포가 있다.")
  n_missing_donor <- sum(is.na(donor) | donor == "")
  log(sprintf("cells %d, clusters %d, donors %d (missing donor %d), studies %d, genes %d", nrow(fitcells), length(unique(cl)),
              length(unique(donor[!is.na(donor)])), n_missing_donor, length(unique(study)), length(symbols)))
  donor_key <- ifelse(is.na(donor) | donor == "", "__missing_donor__", donor)
  study_of <- tapply(study, donor_key, function(x) paste(sort(unique(x)), collapse = "|"))

  # ---- 2. module 정의와 스트리밍 집계
  st <- utils::read.delim(inputs[["state_table"]], stringsAsFactors = FALSE)
  ctrl <- utils::read.delim(gzfile(inputs[["controls"]]), stringsAsFactors = FALSE)
  def <- mdu_module_definition(st, ctrl, symbols)
  avail <- mdu_panel_availability(st, utils::read.delim(inputs[["panel"]], stringsAsFactors = FALSE, na.strings = character(0)))
  nct <- table(ctrl$module)
  chk <- c(setNames(st$pos_n_controls, paste0(st$state, "+")), setNames(st$neg_n_controls, paste0(st$state, "-")))
  chk <- chk[chk > 0]
  if (!identical(as.integer(nct[names(chk)]), as.integer(chk))) stop("controls.tsv 의 대조군 수가 state_table 과 다르다.")
  log(sprintf("module definition: %d states, %d panel genes (duplicate-symbol genes %d), %d control columns",
              length(def$states), length(def$genes), sum(Matrix::colSums(def$P) > 1), length(def$ctrl_cols)))
  gkey <- paste(cl, donor_key, sep = "\r"); ug <- unique(gkey); gid <- match(gkey, ug)
  acc <- mdu_stream_h5(inputs[["h5ad"]], def, gid, length(ug), chunk = a$chunk, log = log)
  grp <- data.frame(cluster = sub("\r.*$", "", ug), donor = sub("^.*\r", "", ug), n_cells = acc$n, stringsAsFactors = FALSE)
  grp$study <- unname(study_of[grp$donor])
  log(sprintf("streamed: %d (cluster, donor) groups, cells %d", nrow(grp), sum(grp$n_cells)))

  # ---- 3. pooled 를 저장 기준과 대조 (theta 는 저장값)
  ref <- utils::read.delim(gzfile(inputs[["reference"]]), stringsAsFactors = FALSE)
  units <- as.character(sort(as.integer(unique(cl))))
  pool <- function(m) { s <- rowsum(m, grp$cluster); s[units, , drop = FALSE] / as.numeric(tapply(grp$n_cells, grp$cluster, sum)[units]) }
  Ap <- pool(acc$pos); Bp <- pool(acc$neg)
  val <- mdu_validate_reference(Ap, Bp, ref, tol = 1e-9)
  theta <- val$theta
  log(sprintf("reference check: max abs diff A %.3g, B %.3g, M %.3g (tol %.0e)", val$max_abs_diff[["A"]], val$max_abs_diff[["B"]], val$max_abs_diff[["M"]], val$tol))
  meas_ref <- ref[ref$state %in% def$states, ]
  measured <- as.vector(tapply(toupper(as.character(meas_ref$measured)) == "TRUE", meas_ref$state, all)[def$states])
  if (!identical(measured, as.vector(def$measured))) stop("저장 기준의 measured 와 패널 상태 표가 다르다.")

  # ---- 4. donor 층
  Ad <- acc$pos / grp$n_cells; Bd <- acc$neg / grp$n_cells
  k <- def$states[def$measured]
  dl <- do.call(rbind, lapply(k, function(s) data.frame(cluster = grp$cluster, donor = grp$donor, study = grp$study, n_cells = grp$n_cells,
                                                         state = s, A = Ad[, s], B = Bd[, s], theta = unname(theta[s]), stringsAsFactors = FALSE)))
  dl$M <- dl$A - pmax(0, dl$B - dl$theta)
  real <- dl$donor != "__missing_donor__"
  missing_donor_cells <- grp[grp$donor == "__missing_donor__", ]
  dsum <- mdu_donor_summary(dl[real, setdiff(names(dl), "M")])
  boot_type <- list(); boot_pair <- list(); boot_win <- list(); boot_meta <- list(); sa <- list(); sp <- list(); lo <- list()
  for (u in units) {
    r <- grp$cluster == u & grp$donor != "__missing_donor__"
    Au <- Ad[r, k, drop = FALSE]; Mu <- mdu_marker_M(Au, Bd[r, k, drop = FALSE], theta); rownames(Au) <- rownames(Mu) <- grp$donor[r]
    b <- mdu_bootstrap_cluster(Au, Mu, R = a$replicates, seed = a$seed + as.integer(u))
    boot_type[[u]] <- cbind(cluster = u, b$type); boot_pair[[u]] <- cbind(cluster = u, b$pairs); boot_win[[u]] <- cbind(cluster = u, b$winner)
    boot_meta[[u]] <- data.frame(cluster = u, status = b$status, n_donors = b$n_donors, n_replicates = b$n_replicates, seed = b$seed,
                                 tie_draw_frac = b$tie_draw_frac, n_draws_available = b$n_draws_available)
    t2 <- mdu_sign_tests_cluster(Au, Mu, u); sa[[u]] <- t2$a_support; sp[[u]] <- t2$pairwise
    lo[[u]] <- cbind(cluster = u, mdu_leave_one_out(Mu, rownames(Mu), "donor"))
    if (length(unique(grp$study[r])) > 1) lo[[u]] <- rbind(lo[[u]], cbind(cluster = u, mdu_leave_one_out(Mu, grp$study[r], "study")))
  }
  bind <- function(x) do.call(rbind, unname(x))
  tests <- mdu_adjust(rbind(cbind(bind(sa)[, c("cluster", "family", "state")], state_t = NA_character_, bind(sa)[, -(1:3)]),
                            setNames(bind(sp)[, c("cluster", "family", "state_s", "state_t", names(bind(sp))[-(1:4)])],
                                     c("cluster", "family", "state", "state_t", names(bind(sp))[-(1:4)]))))
  names(tests)[names(tests) == "state"] <- "state_s"
  log(sprintf("bootstrap %d replicates per cluster; sign tests: A_support %d, pairwise %d", a$replicates, sum(tests$family == "A_support"),
              sum(tests$family == "pairwise")))

  # 기존 판정과 비교할 0 점 unassigned 재현 (안정성 판독이며 새 주석이 아니다)
  win_pool <- t(sapply(units, function(u) { m <- mdu_winner_mass(val$M[u, k]); c(cluster = u, pooled_winner = paste(names(m)[m > 0], collapse = "|")) }))
  win_donor <- sapply(units, function(u) { x <- dsum[dsum$cluster == u, ]; m <- mdu_winner_mass(setNames(x$M_donor_mean, x$state)); paste(names(m)[m > 0], collapse = "|") })
  winners <- data.frame(win_pool, equal_donor_winner = unname(win_donor), stringsAsFactors = FALSE, row.names = NULL)

  # study 별 기술표
  study_tab <- stats::aggregate(cbind(A = dl$A[real], M = dl$M[real]), list(cluster = dl$cluster[real], study = dl$study[real], state = dl$state[real]), mean)
  names(study_tab)[4:5] <- c("A_donor_mean", "M_donor_mean")
  study_n <- stats::aggregate(list(n_donors = rep(1L, sum(real)), n_cells = dl$n_cells[real]),
                              list(cluster = dl$cluster[real], study = dl$study[real], state = dl$state[real]), sum)
  study_tab <- merge(study_n, study_tab, by = c("cluster", "study", "state"))

  # 유전자 근거 (donor 별)
  gm <- acc$panel_log / grp$n_cells; dm <- acc$det / grp$n_cells; pcm <- acc$pos_ctrl / grp$n_cells; ncm <- acc$neg_ctrl / grp$n_cells
  ge <- do.call(rbind, lapply(k, function(s) {
    rows <- list()
    for (side in c("pos", "neg")) {
      gs <- if (side == "pos") def$pos[[s]] else if (def$neg_mode[[s]] == "ok") def$neg[[s]] else character(0)
      cm <- if (side == "pos") pcm[, s] else ncm[, s]
      for (g in gs) rows[[length(rows) + 1]] <- data.frame(cluster = grp$cluster, donor = grp$donor, state = s, side = side, gene = g,
                                                         gene_mean = gm[, g], ctrl_mean = cm, detect_frac = dm[, g], stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }))
  ge <- ge[ge$donor != "__missing_donor__", ]
  gsum <- mdu_gene_summary(ge)

  # ---- 5. 적합 표 저장과 hash 고정
  fit_files <- list("donor_scores.tsv.gz" = dl, "cluster_type_donor_summary.tsv.gz" = dsum, "bootstrap_type_ci.tsv.gz" = bind(boot_type),
                    "bootstrap_pair_contrast_ci.tsv.gz" = bind(boot_pair), "bootstrap_winner_stability.tsv.gz" = bind(boot_win),
                    "bootstrap_cluster_status.tsv" = bind(boot_meta), "sign_tests.tsv.gz" = tests, "leave_one_out.tsv.gz" = bind(lo),
                    "study_descriptive.tsv.gz" = study_tab, "gene_evidence_donor.tsv.gz" = ge, "gene_evidence_summary.tsv.gz" = gsum,
                    "winner_reproduction_zero_unassigned.tsv" = winners, "panel_availability.tsv" = avail, "missing_donor_cells.tsv" = missing_donor_cells,
                    "reference_check.tsv" = data.frame(quantity = names(val$max_abs_diff), max_abs_diff = unname(val$max_abs_diff), tol = val$tol))
  for (nm in names(fit_files)) mdu_write(fit_files[[nm]], file.path(out, nm))
  mdu_save_lossless(fit_files, file.path(out, "fit_tables_lossless.rds"))
  fit_sha <- vapply(c(names(fit_files), "fit_tables_lossless.rds"), function(nm) mdu_sha(file.path(out, nm)), character(1))
  fit_man <- list(status = "fit_outputs_fixed_before_author_labels", written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  file_sha256 = as.list(fit_sha), replicates = a$replicates, seed_base = a$seed, seed_rule = "seed_base + cluster id",
                  reference_max_abs_diff = as.list(val$max_abs_diff), tolerance = val$tol, numerical_zero_tol = MDU_NUM_TOL, winner_tie_tol = MDU_TIE_TOL, n_groups = nrow(grp),
                  author_labels_read = FALSE)
  writeLines(jsonlite::toJSON(fit_man, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "fit_manifest.json"))
  log("fit outputs fixed:", length(fit_sha), "files")

  # ---- 6. 사후 저자 대조 (여기서 처음 저자 라벨을 읽는다)
  f <- hdf5r::H5File$new(inputs[["h5ad"]], "r")
  id2 <- as.character(f[[paste0("obs/", idx_key)]]$read()); author <- ct_h5_obs(f, "AuthorCellType"); f$close_all()
  mdu_check_ids(id2, cell_id)
  ro <- utils::read.delim(gzfile(inputs[["readout"]]), stringsAsFactors = FALSE, colClasses = c(unit = "character", call = "character"))
  ro <- ro[ro$fit == "descriptive" & ro$readout == "M_new_simple", c("unit", "call", "status")]
  if (!setequal(ro$unit, units) || anyDuplicated(ro$unit)) stop("기존 판정 표의 unit 이 51 군집과 다르다.")
  map <- utils::read.delim(inputs[["mapping"]], stringsAsFactors = FALSE, na.strings = character(0))
  au <- mdu_author_audit(data.frame(cell = cell_id, cluster = cl, donor = donor_key, author = author, stringsAsFactors = FALSE), ro, map)
  if (au$overall$n_total != nrow(fitcells)) stop("저자 대조 분모가 전체 세포 수와 다르다.")
  hm <- au$cells[au$cells$author %in% c("HSC", "MPP-MkEry", "MPP-MyLy") | au$cells$call %in% c("HSC", "MPP"), ]
  hsc_mpp <- stats::aggregate(list(n_cells = rep(1L, nrow(hm))), list(cluster = hm$cluster, call = hm$call, call_status = hm$call_status,
                                                                         author_label = ifelse(is.na(hm$author), "<NA>", hm$author), donor = hm$donor), sum)
  stab <- mdu_stability_table(ro, au$per_cluster, bind(boot_meta), bind(boot_win))
  stab <- merge(stab, winners, by = "cluster", all.x = TRUE, sort = FALSE)
  stab <- stab[order(as.integer(stab$cluster)), ]
  stab$note <- "boot_freq_existing_call_equal_donor: equal-donor resampling frequency of the existing (cell-weighted) call; reflects donor stability AND donor weighting AND the aggregation order of the negative penalty max(0, B - theta) (pooled B once vs per donor); not a probability correct"
  aud_files <- list("author_cells_outcome.tsv.gz" = au$cells, "author_overall.tsv" = au$overall, "author_per_donor.tsv" = au$per_donor,
                    "author_per_cluster.tsv" = au$per_cluster, "author_per_type.tsv" = au$per_type, "author_per_predicted_ppv.tsv" = au$per_predicted,
                    "author_cluster_donor_label_counts.tsv.gz" = au$cluster_donor_author, "author_hsc_mpp_detail.tsv" = hsc_mpp,
                    "stability_vs_author_agreement.tsv" = stab)
  for (nm in names(aud_files)) mdu_write(aud_files[[nm]], file.path(out, nm))
  log(sprintf("author audit: total %d, evaluable %d, correct %d, wrong %d, held %d, unavailable %d, unmapped %d",
              au$overall$n_total, au$overall$n_evaluable, au$overall$n_correct, au$overall$n_wrong, au$overall$n_held,
              au$overall$n_unavailable, au$overall$n_unmapped_reference))
  now <- vapply(c(names(fit_files), "fit_tables_lossless.rds"), function(nm) mdu_sha(file.path(out, nm)), character(1))
  if (!identical(now, fit_sha)) stop("저자 대조 뒤 적합 표 hash 가 바뀌었다.")

  # ---- 7. manifest
  srcs <- mdu_compare_sources(snap); mdu_write(srcs, file.path(out, "source_hash_start_end.tsv"))
  man <- list(status = if (all(srcs$same)) "complete" else "complete_but_source_changed_during_run",
              command = paste(c("run-r dcta R/branch_biology/18_run_marker_donor_uncertainty.R", argv), collapse = " "),
              execution = list(bootstrap = "mdu_main from the original 18 (argument parsing, snapshot, loading only)",
                               loaded_from_snapshot = setNames(as.list(loaded$sha256), basename(loaded$snapshot_path)),
                               body = "mdu_run_body and helpers taken from the snapshot environment",
                               snapshotted_not_run = "19_view_marker_donor_uncertainty.R (figures drawn later)"),
              source_sha256_at_start = setNames(as.list(snap$sha256), basename(snap$original_path)),
              source_same_at_end = setNames(as.list(srcs$same), basename(srcs$original_path)),
              input_sha256 = as.list(in_sha), fit_file_sha256 = as.list(fit_sha),
              author_file_sha256 = lapply(setNames(names(aud_files), names(aud_files)), function(nm) mdu_sha(file.path(out, nm))),
              n_cells = nrow(fitcells), n_clusters = length(units), n_donors = length(unique(donor_key[donor_key != "__missing_donor__"])),
              n_missing_donor_cells = n_missing_donor, replicates = a$replicates, seed_base = a$seed,
              reference_max_abs_diff = as.list(val$max_abs_diff),
              r_version = R.version.string, rscript = file.path(R.home("bin"), "Rscript"), lib_paths = .libPaths(),
              packages = lapply(c(Matrix = "Matrix", hdf5r = "hdf5r", jsonlite = "jsonlite", digest = "digest"), function(p) as.character(utils::packageVersion(p))),
              seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
              interpretation = c("donor-resampling intervals are descriptive, conditional on frozen panel/controls/theta/clusters; not a posterior or calibrated probability",
                                 "sign tests are exploratory: clustering and controls used the same data; not post-selection valid",
                                 "no output is a probability that the marker call is correct; stable wrong panels remain possible"))
  writeLines(jsonlite::toJSON(man, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null"), file.path(out, "manifest.json"))
  writeLines(capture.output(sessionInfo()), file.path(out, "environment_sessionInfo.txt"))
  log("done; sources unchanged during run:", all(srcs$same))
  invisible(TRUE)
}


if (!interactive() && sys.nframe() == 0L) quit(status = mdu_main(commandArgs(trailingOnly = TRUE)))
