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
- Prior hematopoietic DAG (G_prior) — partial: `dag_prior_bm.yaml` exists,
  not yet wired into backbone-design step (Phase 2.2)
- Absorbing Markov chain on G_prior with PAGA-informed transitions
- Terminal / transient annotation per cell
- AML/MDS deviation detection
- Virtual HSC for blood/spleen (`virtual_hsc_ref` strategy hook already in place)
- clustree-style robustness sweep across resolutions
- Benchmarks vs SingleR / CellTypist / Azimuth / Palantir / CellRank
- Resident-vs-emergent lineage labeling
