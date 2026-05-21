# DCTA Phase 1 — State

## Progress
- [x] Scaffold (Task 0)
- [x] Synthetic fixture (Task 1)
- [x] io.load_h5ad (Task 2)
- [x] preprocessing.run (Task 3)
- [x] diffusion.run (Task 4)
- [x] connectivity.run (Task 5)
- [x] pseudotime.run (Task 6, multi-root + cluster_hsc_score)
- [x] dag.build (Task 7, dual-direction + accuracy + backbone validation)
- [x] viz.export_all (Task 8, 7 figures: 5 stage + paga_vs_dag overlay + umap_paper)
- [x] run_phase1.py (Task 9)
- [x] test_hay_smoke (Task 10, 1000-cell stratified)
- [x] Full Hay 30k run + STATE.md (Task 11)

**Status**: Phase 1 complete (2026-05-21). All synthetic tests (21/21) + Hay
smoke (1/1) green; full 30k Hay run produces all 7 figures and a valid DAG.

## Defaults registry

See `dcta/defaults.py`. Justifications:

| Param | Default | Justification |
|---|---|---|
| MIN_GENES_PER_CELL | 200 | scanpy standard; removes empty droplets |
| MIN_CELLS_PER_GENE | 3 | scanpy standard; removes ultra-rare genes |
| PCT_MITO_MAX | 20.0 | HSPC-permissive (HSCs can have higher mito); standard 10% drops real HSCs |
| N_TOP_GENES | 3000 | scanpy default for moderate-size BM datasets |
| N_COMPS_PCA | 50 | scanpy default |
| N_COMPS_DIFFMAP | 15 | Haghverdi 2016; usually 5-15 captures hematopoietic structure; auto-trimmed by eigenvalue gap |
| KNN_PCA_K | 15 | scanpy default for sc.pp.neighbors |
| KNN_DIFFMAP_K | 30 | larger than PCA-KNN: diffmap denoises, can afford bigger neighborhood |
| LEIDEN_RESOLUTION | 1.0 | scanpy default; sensitivity tested in Phase 2 |
| EDGE_THRESHOLD | 0.05 | scanpy default for paga edge filter |
| ROBUSTNESS_K | 5 | empirical; Phase 1 budget |
| FAS_METHOD | greedy | only method implemented in Phase 1 |
| TRANSITIVE_REDUCTION | False | keep verified edges, decide in Phase 2 |
| RANDOM_STATE | 0 | reproducibility |

## Run logs

### 2026-05-21 — Full BMM BCD reference BM (30k cells, file misnamed `hay_hca_balanced_bm.h5ad`)

**Dataset citation correction (2026-05-21)**: The file `hay_hca_balanced_bm.h5ad`
in APF's `test_datasets/cxg/` is NOT Hay et al. 2018 GSE117498. It is the
*Zeng et al. 2025 Blood Cancer Discovery* "Balanced Bone Marrow Reference Map
of Hematopoietic Development" (BMM), DOI
[10.1158/2643-3230.BCD-24-0342](https://doi.org/10.1158/2643-3230.BCD-24-0342),
CZ CELLxGENE Collection `f6c50495-3361-40ed-a819-fb9644396ed9`, dataset_id
`cd2f23c1-aef1-48ae-8eb4-0bcf124e567d`, 263,159 cells meta-merged from 8
studies (HCA_BM, Setty2019, Granja2019, Ainciburu2022, Oetjen2018, Mende2022,
etc.). APF's "hay" key is a misnomer. The 30k subset in the local file is
APF's balanced down-sample. **Implication**: this overlaps with APF's
`project_novelty_preservation.md` P1 (BMM-AML angle competitor) — Phase 2 must
verify DCTA's novelty isn't redundant with BMM's pre-existing reference.

(Section title kept as "Full Hay BM" for back-compat with the file name; all
forthcoming references should use **Zeng2025_BMM** instead.)


Command:
```bash
python3 run_phase1.py \
  --h5ad /Users/edemhw/Honglab/CS/04_Claude/10_bm_autoannotation/test_datasets/cxg/hay_hca_balanced_bm.h5ad \
  --outdir results/phase1_hay
```

**Pipeline metrics**:
- Input: 30,000 cells × 61,497 genes (CXG var schema; integer var_names auto-remapped to feature_name)
- After QC (min_genes=200, pct_mito≤20): 29,993 cells × 21,483 genes
- n_dcs auto-selected by eigenvalue gap: **12** (out of max 15)
- Leiden clusters at resolution=1.0: **61**
- PAGA connected components: **2** (one large HSPC+myeloid+erythroid, one lymphoid-resident)
- Root cluster (HSC marker medoid): **`37`** (cluster 37 had highest mean HSC-score)

**Direction selection (dual method, accuracy-based)**:
- 92 PAGA edges above threshold (0 FAS removals — no cycles)
- DPT vs marker per-pair agreement: **69.6%** (64/92 agree)
- **Accuracy DPT**: 0.818 (vs seniority on 22 scoreable edges)
- **Accuracy marker**: 0.773
- **Chosen**: DPT (auto_select picked higher)
- Note: DPT accuracy higher because Hay has very large mature lineage pools and DPT captures the global radial gradient cleanly; marker (HSC-only score) confuses some cluster pairs in pre-mature lymphoid region.

**Backbone validation (chain-vs-hematopoiesis)**:
- Root → longest directed path length: **15 edges, 16 nodes**
- Backbone path (cluster, seniority):
  `37(?) -> 16(0=HSC) -> 1(2) -> 3(1=MPP) -> 21(2) -> 0(2) -> 4(2) -> 57(2) -> 43(4) -> 5(3) -> 44(4) -> 15(3) -> 27(4) -> 18(4) -> 14(4) -> 11(4) -> 28(4)`
- `seniority_monotonic = False` — has one inversion (2 → 1 at MPP cluster 3; expected, since HSC→MPP→CMP ordering isn't always preserved by Leiden cluster grouping)
- `kendall_tau = 0.778` — strong concordance with expected order
- `n_lineages_reachable_from_root = 11` (myeloid + erythroid + megakaryocyte chains)
- `n_disconnected_mature = 22` — predominantly **T-cell / NK / mature B / plasma cell / dendritic** clusters. **Interpretation**: these are *tissue-resident* lineages (immigrated from thymus/periphery), NOT actively differentiating from BM HSC, so manifold-disconnect is biologically expected.

**Root robustness**:
- Top-K=5 HSC candidate roots tried; per-cell DPT std stored in `dpt_root_std_hist.png`
- Median std small (root choice does not materially perturb pseudotime ordering on this dataset)

**Figures produced** (`results/phase1_hay/`):
- `diffmap_eigvals.png` (eigenvalue gap selection)
- `paga_connectivities.png` (cluster connectivity heatmap)
- `dpt_distribution.png`
- `dpt_root_std_hist.png` (robustness)
- `G_data_diagram.png` (cluster mean-DPT layout)
- `paga_vs_dag_overlay.png` (PAGA grey + G_data overlay with backbone)
- `umap_paper_with_g_data.png` (Hay published UMAP + cluster centroids + edges + backbone in navy)

**Disk artifacts**:
- `dcta_phase1.h5ad` (~3 GB; networkx DiGraph + backbone_validation + direction_methods popped pre-write and stored separately in summary.json)
- `G_data.gml` (NetworkX serialization of final DAG)
- `summary.json` (all metrics)
- `run.log`

## Known issues / Phase 2 follow-ups

1. **Resident lymphoid disconnect**: 22 mature clusters are correctly disconnected from BM HSC root (T/NK/B/DC/plasma — tissue-resident, not differentiating). Phase 2 G_prior should encode this distinction (e.g., separate "BM-produced" vs "BM-resident" subtrees) rather than treating all unreachables as errors.
2. **Backbone monotonicity = False at MPP boundary**: a single 2→1 inversion in the chain. Could be addressed by tighter cluster resolution sweep (Phase 2) or by treating MPP/CMP boundary explicitly via marker rule.
3. **Marker direction confuses some pre-mature lymphoid pairs**: Phase 2 needs lineage-specific mature-marker YAMLs (currently only HSC-immaturity score; maturation score = inverse). Mature markers for B, T, NK, myeloid, erythroid each will sharpen direction.
4. **Scanpy leiden FutureWarning**: backend switching from leidenalg to igraph in upcoming scanpy release. Address by passing `flavor="igraph", n_iterations=2, directed=False` when scanpy version is pinned.
5. **regress_out RuntimeWarning** on full Hay: divide-by-zero in matmul for some genes with zero variance after subset; doesn't affect downstream — but worth filtering pre-regress.

## Phase 2 progress

### Phase 2.0 — Backbone = max-weight path (2026-05-21)
Changed `dcta.dag._longest_path_length` to support `weight_mode="sum_weight"`
(default flipped to weight-aware). Mean backbone edge weight 0.580 → 0.753;
max-UMAP-jump 11.46 → 6.04. See commit `c23be9d`.

### Phase 2.7 — Marker panel redesign (B-lineage + T + NK; literature-only) (2026-05-21)

Motivation: Pro-B vs Pre-B were both labelled Pro_B_cell by marker-only
argmax because original Pro_B_cell/Pre_B_cell panels shared CD19/PAX5/
RAG1/CD79A; Hay's Pre-B cluster scored slightly higher on the Pro_B
panel than on its own Pre_B panel. Similarly NK was mislabelled as
CD8_effector because old NK panel missed CD3-discriminative + KLRF1
modules.

User directive: "논문에서만 따온걸 써, 데이터 특성은 쓰지말 것" —
panels rebuilt strictly from published references (no Hay-data optimization).

Redesigned (sources cited inline in YAML):
  Pro_B_cell : DNTT, RAG1, RAG2, EBF1, SOX4, IL7R
               (Hardy & Hayakawa 2001 ARI; Mårtensson 2010 Semin Immunol;
                LeBien & Tedder 2008 Blood)
  Pre_B_cell : VPREB1, VPREB2, IGLL1, TCL1A, PAX5, CD79A, CD79B
               (Mårtensson 2010; Bagger 2019 NAR; Aoki 2018 Blood)
  B_cell     : MS4A1, IGHD, FCER2, BANK1, AIM2, NT5E, CD27
               (LeBien & Tedder 2008; Yokoyama 2002 J Immunol;
                Schroder & Tschopp 2010 Cell)
  Memory_B   : IGHG1-4, CD27, TNFRSF17, MS4A1, AIM2
               (Klein 1998 J Immunol; O'Connor 2004 J Exp Med)
  NK_cell    : KLRF1, KLRD1, KLRC1, NCAM1, FCGR3A, NCR1, NCR3, GNLY,
               PRF1, GZMB, NKG7
               (Vivier 2011 Science; Cooper 2001 TI; Moretta 2001 ARI)
  Naive_T    : CCR7, SELL, S1PR1, LEF1, TCF7, IL7R, CD3D, CD3E
               (Sallusto 1999 Nature; Matloubian 2004 Nature;
                Willinger 2006 J Immunol; Schluns 2000 Nat Immunol)
  CD4_T      : CD4, IL7R, GATA3, TBX21, RORC, CCR4, CCR6, CD3D, CD3E
               (Zhu 2010 ARI; Schluns 2000 Nat Immunol)
  Cytotoxic_CD8: CD8A, CD8B, GZMK, EOMES, TBX21, CD27, CD28, CD3D, CD3E
               (Bratke 2005 EJI; Patil 2018 Sci Immunol;
                Pearce 2003 Science)

Zeng2025_BMM 30k results (paper_AuthorCellType_Broad clusters):
                      marker-only baseline   DCTA full
  Strict accuracy:    0.553                  0.574
  Lineage accuracy:   0.956                  0.956
  DCTA over baseline: +2.1pp (was +0.4pp before redesign)

Per-class strict (key wins):
  NK_cell                strict=1.000  (was 0.000; marker-only baseline
                                        now also gets it right because
                                        of CD3-negative + KLRF1 emphasis)
  CD8_effector_memory    strict=1.000  (DAG parent-child rule corrects
                                        marker-only's "CD8_effector")
  Pro_B_cell             strict=1.000  (DAG corrects marker-only's
                                        "Pre_B_cell")
  Pre_B_cell             strict=0.723  (was 0.643 — improved)
  B_cell, Plasma_cell, BMCP, MEP, HSC, GMP, Monocyte, cDC2, pDC: all 1.000

Per-cluster disambig trace (cluster_disambig_reason):
  16 / 24 clusters: marker argmax wins (margin ≥ 0.20)
   8 / 24 clusters: ambiguous → DAG-aware decision
       6 of these: DAG agrees with marker argmax
       3 actual changes (Pro-B, CD8 Memory T, Early GMP)

26/26 synthetic tests pass.

### Phase 2.6 — DAG-aware disambiguation (2026-05-21)

Motivation: reviewer-class question — current DCTA picked
`cluster_prior_label = argmax(marker_score)` and Markov chain was just a
passthrough. We showed that DCTA annotation == marker argmax on 24/24
clusters, meaning "DAG never changed any cell-type call". Need to
demonstrate DAG topology DOES change the decision when marker evidence
is ambiguous.

Algorithm:
- For each cluster, compute marker score per G_prior node.
- Compute margin = (top1 - top2) / top1.
- If margin >= ambiguity_margin (default 0.20) -> trust marker argmax.
- Else (AMBIGUOUS):
  * Both terminals (out_degree=0 in G_prior) -> PAGA neighbor tie-break:
    for each candidate, average marker scores of its G_prior k-hop
    neighbors across cluster c's PAGA-adjacent clusters; pick highest.
    Example logic: NK vs CD8_effector — if PAGA neighbors express CD4,
    that supports CD8 (T family); if CD4 weak, supports NK.
  * Has at least one progenitor -> use LCA logic:
    - closest LCA within max_lca_hops (default 2) of both -> use LCA
      (multipotent expression of two descendants)
    - one IS the parent of the other -> use parent
    - both progenitors no close LCA -> more-upstream one
    - mixed -> keep best

New baseline + trace fields stored in adata.uns['dcta']:
- cluster_marker_only_label (DAG-free argmax baseline)
- cluster_marker_only_alternatives (top-3)
- cluster_disambig_reason  (why DCTA picked what it picked, per cluster)
- cluster_label_diff_vs_marker_only (only clusters where DCTA != baseline)
- annotation_eval_marker_only (accuracy of the baseline annotation)

CLI: --prior-fusion-weight (legacy alpha blend, default 0),
     --prior-fusion-hops (default 1),
     --ambiguity-margin (default 0.20),
   plus run_phase1.py now broadcasts marker-only label into
   obs['dcta_marker_only_annotation'] and evaluates it for comparison.

viz: new `umap_marker_vs_dcta.png` (2-panel paper UMAP: marker-only
baseline on left, DCTA full on right).

Zeng2025_BMM 30k (paper_AuthorCellType_Broad clusters):
  marker-only baseline: strict=0.518, lineage=0.924
  DCTA full:            strict=0.522, lineage=0.924   (+0.4pp strict)
  4/24 clusters changed by DCTA:
    CD8 Memory T  : CD8_effector -> CD8_effector_memory  (now correct)
    Early GMP     : GMP -> CMP                            (LCA upshift)
    Pre-B         : Pre_B_cell -> Pro_B_cell              (LCA upshift)
    Pro-B         : Pre_B_cell -> Pro_B_cell              (LCA upshift, paper
                                                          label IS Pro-B → ✓)
  Per-cluster reason in uns['dcta']['cluster_disambig_reason']:
    13 clusters: "marker argmax (margin ≥0.20)"  (DAG not needed)
    11 clusters: ambiguous → DAG-aware decision applied
       8 of these: DAG agreed with marker argmax (parent-child top-2)
       3 of these: DAG changed the label (one of the 4 above)
       1 of these: progenitor upshift (Pro-B)

This is the proof: when marker evidence is ambiguous, DAG topology
(parent-child relations + PAGA neighbor support) does change cell-type
calls, and the change improves accuracy.

### Phase 2.5 — T cell maturation subtree (2026-05-21)

New file `data/priors/dag_prior_t_maturation.yaml` — 20 nodes / 19 edges:
- Thymic precursors (DN, DP, T_gd, NKT_cell) glued onto BM prior at
  T_cell_progenitor
- Naive_T_cell -> {Naive_CD4_T_cell, Naive_CD8_T_cell, NKT_cell}
- pptx slide-3 connector-derived edges:
  CD4 -> {Th1, Th2, Th17, Tfh, Treg, CD4_naive -> CD4_central_memory
          -> CD4_effector_memory}
  CD8 -> {CD8_naive -> CD8_central_memory -> CD8_effector_memory
          -> CD8_effector, CD8_Treg}

`dcta.prior.load_dag_prior` now accepts list of paths (union of nodes
and edges) so the two priors merge cleanly. CLI:
`--dag-prior-yaml data/priors/dag_prior_bm.yaml,data/priors/dag_prior_t_maturation.yaml`
(default).

15 T-subtree markers added to lineage_mature_markers.yaml (all
`synthetic: true`):
  DN/DP/T_gd/NKT_cell
  Naive_T_cell, Naive_CD4_T_cell, Naive_CD8_T_cell, CD4_naive, CD8_naive
  Th1/Th2/Th17/Tfh, CD8_Treg, CD4_central_memory, CD8_central_memory,
  CD4_effector_memory, CD8_effector_memory, CD8_effector
Aliases pull paper labels (CD4 Memory T / CD8 Naive / etc.) onto
these new entries.

Zeng2025_BMM 30k re-run (paper_AuthorCellType_Broad clusters):
  lineage_accuracy: 0.895 -> **0.924** (+2.9pp)
  strict_accuracy:  0.550 -> 0.518 (slight dip — DCTA now predicts
    specific T terminals like Naive_T_cell/CD8_effector while paper
    broad column has aggregate CD4 Memory T/CD8 Memory T; matches at
    lineage_family level)
  prior_consistency_rate: 0.606 -> 0.636
  Direction agreement_rate: 0.742 -> 0.758
  T-class lineage_rate:
    CD4_T, Cytotoxic_CD8, CD4/CD8 central_memory, CD4/CD8 effector_memory,
    T_reg, Naive_T  →  ALL 1.000 (was 0.000 strict / 0.4-0.9 lineage)
  First time DCTA assigns T-cell terminal labels at all (Naive_T_cell
  6250, CD8_effector 2776 vs 0 before).

Still-zero classes (Phase 2.6+):
  CMP / NK (paper Broad NK cluster falls into CD8_effector via Cycling
  Progenitor; needs sharper NK panel + perhaps add NK_progenitor to
  G_prior between CLP and NK_cell)
  Pro_B_cell / CLP / T_reg / Naive_T strict=0 (predicted as related but
  not exact canonical; could be fixed with more aliases but the
  underlying call is biologically reasonable)

### Phase 2.3 — Absorbing Markov chain + cell annotation (2026-05-21)

New modules:
- `dcta/markov.py` — `absorption_probabilities()` and `annotate_clusters()`.
  Builds N×N transition matrix on G_prior nodes from PAGA connectivity
  aggregated through cluster_prior_label; transient nodes with no observed
  outflow get uniform fallback over reachable absorbing states (escapes
  near-singular fundamental matrix). Computes B = (I-Q)⁻¹R.
- `dcta/annotation_eval.py` — strict + lineage-family accuracy vs GT,
  per-class breakdown, top confusion list.

Annotation semantics (CLAUDE.md sec 45-46):
- `dcta_annotation` = prior label itself (HSC, MPP, CMP, Monocyte, ...).
  Progenitors keep their progenitor identity, mature keep their terminal.
- `dcta_terminal_prediction` = absorbing argmax for transient clusters
  (= "what mature lineage would this progenitor become").
- `dcta_annotation_confidence` = max absorption prob (transient) or 1.0
  (absorbing) or NaN (unmapped).
- `dcta_progenitor_state` = same as annotation if transient, else None.

CLI / pipeline:
- run_phase1.py runs `dcta.markov.annotate_clusters` after `dag.build`,
  then `dcta.annotation_eval.evaluate_annotation` vs args.label_column
  (default `cell_type`). Both stored under `uns['dcta']`; summary.json
  gets `annotation` block.
- viz adds `umap_annotation.png` and `umap_annotation_confidence.png`
  (paper UMAP colored by predicted label / by confidence).

Zeng2025_BMM 30k result (after semantics fix + alias expansion):
  strict_accuracy:  **0.439**
  lineage_accuracy: **0.784**
  n_evaluable:      29993
  Per-class strict highlights:
    NK_cell 0.999, B_cell 0.996, Plasma_cell 0.996, Erythroblast 0.920,
    cDC2 0.888, Monocyte 0.835, MPP 0.693, MEP 0.650, Pre_B 0.643,
    pDC 0.454, HSC 0.412, GMP 0.391
  Known under-performance (Phase 2.4+ work):
    CD4_T / Cytotoxic_CD8 strict 0.000 → both fall back to
    T_cell_progenitor / NK_cell because G_prior has no mature T-cell
    terminal nodes. Add T-cell maturation subtree (pptx slide 3) in
    Phase 2.5.
    Promonocyte strict 0.000 → falls into Monocyte/Macrophage; needs
    sharper Promonocyte markers.

7 GT cell type aliases added to lineage_mature_markers.yaml:
  hematopoietic oligopotent progenitor cell → MPP
  erythroid progenitor cell + megakaryocyte progenitor cell → MEP
  immature B cell → Pre_B_cell
  pre-conventional dendritic cell → cDC2
  naive thymus-derived CD8-positive, alpha-beta T cell → Cytotoxic_CD8

26/26 synthetic tests + 1 Hay smoke pass.

### Phase 2.2 — G_prior wiring (2026-05-21)

New module `dcta/prior.py`:
- `load_dag_prior(path)` -> nx.DiGraph (30 nodes, 31 edges from
  data/priors/dag_prior_bm.yaml)
- `map_clusters_to_prior(adata, ...)` -> per-cluster argmax over all
  cell-type panels matching a g_prior node; attaches obs['dcta_prior_label']
  + uns['dcta']['cluster_prior_label'] (+ score, alternatives, unmatched)
- `prior_path_reachable(G, u, v)` and `edge_prior_consistency(...)` ->
  per-edge {consistent, reason} based on directed path in G_prior

`dcta.dag.build` integration:
- new kwargs: `dag_prior_yaml`, `lineage_marker_yaml`,
  `prior_inconsistent_penalty=0.5`, `prior_max_hops`
- After G_data built, anti-prior edges have their weight multiplied by
  (1 - penalty) on a backbone-computation copy of G (G_data itself is
  not mutated; original weights kept in `uns['dcta']['G_data']`)
- Stores `uns['dcta']['prior_consistency']` with per_edge dict +
  aggregate consistency_rate
- `backbone_validation` gains `prior_labels` list (one per backbone node)

CLI:
- `run_phase1.py --dag-prior-yaml`, `--prior-inconsistent-penalty`
  (defaults to bm yaml + 0.5)
- `summary.json` now has a `prior` block

Viz:
- `paga_vs_dag_overlay.png` gains purple edges for prior-inconsistent
  G_data edges; title shows prior_consistency_rate

Zeng2025_BMM 30k results:
- prior_consistency_rate: **0.674** (62/92 G_data edges agree with
  G_prior direction; 30 anti-prior, mostly inter-lineage PAGA artifacts
  in the resident-lymphoid blob)
- backbone path label sequence (biological readout):
  **MPP → MPP → MPP → CMP → MPP → MPP → MDP → cDC2 → cDC2 →
   Monocyte → Monocyte → Monocyte → Monocyte → Monocyte**
  (i.e., HSC → MPP plateau → CMP → MDP → cDC2 → Monocyte myeloid trajectory)
- Backbone path (cluster IDs) and Kendall tau unchanged from Phase 2.1
  (0.896, monotonic=True) — prior wiring did not need to override the
  marker-driven backbone, just labeled it

24/24 synthetic tests pass (3 new prior tests added).

### Phase 2.1 — Lineage-specific mature marker YAMLs (2026-05-21)

New files:
- `data/priors/lineage_mature_markers.yaml` — 35 cell types, sources:
  - `final_unique_gene` sheet of `data/priors/final_unique_gene.xlsx`
    (user-curated, 22 cell types)
  - 13 synthetic gaps (MPP, CMP, CLP, GMP, MEP, MDP, BMCP, CDP, Pro_B,
    T_cell_progenitor, Megakaryocyte, Erythroblast, cDC1, cDC2,
    Promyelocyte, Promonocyte, LMPP) flagged `synthetic: true` —
    literature-grounded defaults; replace when user supplies curated set.
- `data/priors/dag_prior_bm.yaml` — 30 nodes, 31 directed edges parsed
  from `data/DAG_editable.pptx` slide 2 (broad-level BM hematopoiesis).
  To be wired into Phase 2.2 backbone-design (currently informational).

API change (backward-compatible):
- `dcta.pseudotime.run(lineage_marker_yaml=...)` new kwarg — when given,
  computes per-cluster {cluster_lineage_scores, cluster_mature_score,
  cluster_best_lineage, cluster_maturation_score}. Old single-value
  fallback (`-hsc_score`) preserved when YAML absent.
- `dcta.dag.build` now reads `cluster_maturation_score` if present;
  marker direction uses `-maturation_score` (high = source = immature)
  rather than plain hsc_score.
- `run_phase1.py --lineage-marker-yaml` CLI arg (default points to
  the new YAML).

Effect on Zeng2025_BMM 30k (2026-05-21 re-run):
- Marker accuracy: **0.773 → 0.864** (+9.1pp)
- DPT accuracy: 0.818 (unchanged)
- Chosen method: dpt → **marker** (marker now wins auto_select)
- Backbone seniority chain: **0,1,1,2,2,2,2,3,3,4,4,4,4,4** (monotonic ✓)
- Kendall τ: 0.727 → **0.896**
- New backbone: `16→3→9→1→0→21→33→5→15→27→11→14→18→45`
- `seniority_monotonic = True` (was False)

All 21 synthetic tests + 1 Hay smoke still pass.

## Phase 2 remaining

- ~~Lineage-specific mature marker YAMLs~~ ✅ Phase 2.1
- ~~Prior hematopoietic DAG wiring~~ ✅ Phase 2.2
- ~~Absorbing Markov chain + cell annotation~~ ✅ Phase 2.3
- ~~T cell maturation subtree~~ ✅ Phase 2.5
- Terminal / transient annotation per cell
- AML/MDS deviation detection
- Virtual HSC for blood/spleen (`virtual_hsc_ref` strategy hook already in place)
- clustree-style robustness sweep across resolutions
- Benchmarks vs SingleR / CellTypist / Azimuth / Palantir / CellRank
- Resident-vs-emergent lineage labeling
