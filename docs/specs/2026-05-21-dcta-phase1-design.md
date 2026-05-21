# DCTA Phase 1 Design — Data-driven DAG (G_data) on BM

**Date**: 2026-05-21
**Project**: DCTA (DAG-Consistent Trajectory Annotation)
**Phase**: 1 — Proof of concept on a single BM dataset, myeloid lineage focus
**Status**: Design approved by user, ready for implementation planning
**Source spec**: `~/Desktop/CLAUDE.md` (project root vision)

---

## 1. Scope

### In scope (Phase 1)

Build the **data-driven half** of DCTA — the pipeline that produces `G_data`,
a directed acyclic graph over Leiden clusters derived purely from the manifold,
on a single BM dataset.

Pipeline endpoints:
```
hay_hca_balanced_bm.h5ad
  → preprocessing (QC + cell cycle regress-out)
  → PCA → KNN → diffusion map (eigenvalue-gap selected n_dcs)
  → diffusion-based KNN → Leiden → PAGA
  → root selection → DPT (multi-root aware)
  → edge threshold + DPT direction → FAS cycle removal
  → G_data (NetworkX DiGraph, optionally transitively reduced)
```

### Out of scope (Phase 2+)

- Prior hematopoietic DAG (`G_prior`) construction
- YAML marker rule schema beyond the HSC root rule
- Per-cluster marker scoring against `G_prior` nodes
- Structural consistency check between `G_data` and `G_prior`
- Absorbing Markov chain on `G_prior`
- Terminal / transient annotation
- AML/MDS deviation detection
- Virtual HSC for blood/spleen (API hook only)
- clustree-style robustness sweep across resolutions
- Benchmarks vs SingleR / CellTypist / Azimuth / Palantir / CellRank

---

## 2. Project layout

```
14_dcta/
├── CLAUDE.md                       # copied from ~/Desktop/CLAUDE.md (project-local)
├── README.md
├── STATE.md                        # progress + magic-number registry
├── docs/
│   └── specs/
│       └── 2026-05-21-dcta-phase1-design.md   # this file
├── dcta/                           # flat library, no packaging
│   ├── __init__.py
│   ├── defaults.py                 # all default parameters in one place
│   ├── io.py                       # AnnData load, label normalization
│   ├── preprocessing.py            # QC, normalize, HVG, cell-cycle regress-out
│   ├── diffusion.py                # PCA → KNN → diffusion map, eigvalue-gap
│   ├── connectivity.py             # diffusion-KNN → Leiden → PAGA
│   ├── pseudotime.py               # root selection + DPT (multi-root aware)
│   ├── dag.py                      # edge threshold + direction → FAS → optional TR
│   └── viz.py                      # diffmap UMAP, PAGA plot, G_data diagram
├── data/
│   └── priors/
│       └── hsc_markers.yaml        # Phase 1: HSC root identification only
├── results/
│   └── phase1_hay/                 # outputs
├── tests/
│   ├── conftest.py
│   ├── test_synthetic.py           # synthetic AnnData unit tests per stage
│   └── test_hay_smoke.py           # 1000-cell subset end-to-end smoke
└── run_phase1.py                   # orchestrator CLI

(NO pyproject.toml / setup.cfg / src/ layout — flat module by design)
```

---

## 3. Data contract

All public stage functions follow the scanpy idiom:

```python
def run(adata: AnnData, *, copy: bool = False, **params) -> AnnData: ...
```

- In-place mutation when `copy=False` (default), copy when `copy=True`
- Intermediate state is namespaced under `adata.uns['dcta']` to avoid collision
  with scanpy / other tools
- Per-cell stage outputs go into `adata.obs` with `dcta_` prefix
- Per-cluster outputs go into `adata.uns['dcta']`

### 3.1 Field map

| Stage | adata mutations |
|---|---|
| `preprocessing` | `obs['n_counts']`, `obs['pct_mito']`, `obs['S_score']`, `obs['G2M_score']`, `var['highly_variable']`, raw counts preserved in `adata.raw` |
| `diffusion` | `obsm['X_pca']`, `obsm['X_diffmap']`, `uns['dcta']['diffmap_eigvals']`, `uns['dcta']['n_dcs']` |
| `connectivity` | `obsp['dcta_knn']`, `obs['dcta_leiden']`, `uns['dcta']['paga']` (dict: connectivities, cluster_sizes, cluster_centroids) |
| `pseudotime` | `obs['dcta_dpt']`, `obs['dcta_component']`, `uns['dcta']['root_strategy']`, `uns['dcta']['root_cells']`, `uns['dcta']['root_cluster']`, `uns['dcta']['component_count']`, `uns['dcta']['dpt_root_std']`, **`uns['dcta']['cluster_hsc_score']`**, **`uns['dcta']['cluster_maturation_score']`** (per-cluster mean marker scores) |
| `dag` | `uns['dcta']['G_data']` (edge list with weight, direction, retained flag), `uns['dcta']['cycle_removal_log']`, **`uns['dcta']['direction_methods']`** (dict: per-method direction map + accuracy), **`uns['dcta']['direction_method_chosen']`** (string: "marker" \| "dpt" \| "marker (no GT)"), **`uns['dcta']['biological_seniority_map']`** (label→rank, if cell_type column present) |
| `viz` (Sec 8) | writes `paga_vs_dag_overlay.png` comparing PAGA layout with G_data directed edges |

---

## 4. Pipeline (`run_phase1.py`)

```python
import dcta

adata = dcta.io.load_hay(path="../test_datasets/cxg/hay_hca_balanced_bm.h5ad")
adata = dcta.preprocessing.run(adata, regress_cell_cycle=True)
adata = dcta.diffusion.run(adata, n_comps=15)
adata = dcta.connectivity.run(adata, resolution=1.0)
adata = dcta.pseudotime.run(adata, root_strategy="hsc_marker_medoid",
                            hsc_marker_yaml="data/priors/hsc_markers.yaml",
                            robustness_k=5)
adata = dcta.dag.build(adata, edge_threshold=0.05, fas_method="greedy",
                       transitive_reduction=False,
                       direction_strategy="auto_select",
                       seniority_yaml="data/priors/lineage_seniority.yaml")

adata.write_h5ad("results/phase1_hay/dcta_phase1.h5ad")
dcta.viz.export_all(adata, outdir="results/phase1_hay/")
```

All knobs above are CLI flags on `run_phase1.py` with defaults pulled from
`dcta/defaults.py`.

### Why `direction_strategy="auto_select"`?

DPT is a manifold-distance pseudotime, not a biological clock — it can flip
direction on "branch ends" or in datasets where the manifold is not strictly
monotone with maturation. To avoid trusting DPT blindly, the DAG stage
computes **two** candidate directions per edge:

1. **DPT direction**: low cluster-mean DPT → high cluster-mean DPT.
2. **Marker direction**: high cluster-mean HSC-marker score → low HSC-marker
   score (i.e., immature → mature).

When a biological seniority label is available (`adata.obs['cell_type']` mapped
through `lineage_seniority.yaml`), the stage measures each method's per-edge
agreement with the seniority-derived expected direction, and selects the
method with higher agreement to build `G_data`. Both methods' per-edge
decisions are stored for inspection. If no ground-truth column is available,
the marker method is the default (it's grounded in cell biology, not in the
embedding's idiosyncrasies).

---

## 5. Root strategy (multi-root aware)

A single iroot for DPT is fragile. Three concerns drive this section:

1. The HSC cluster contains many cells equally entitled to be "the root"
2. The PAGA graph can be disconnected (rare in healthy BM, plausible in MDS/AML)
3. Blood/spleen require a virtual HSC (Phase 2)

### 5.1 Default: HSC marker medoid

- Score every cell with `hsc_markers.yaml` (Phase 1: a minimal positive marker
  set: CRHBP, AVP, HLF, MEIS1, MLLT3, HOPX — to be finalized in the YAML file)
- The HSC cluster = the Leiden cluster with the highest mean HSC score
- iroot = medoid of the HSC cluster in diffusion-map space (most central HSC,
  reproducible, not random top-1)

### 5.2 Disconnected components

- If PAGA produces ≥2 connected components, pick one root per component
  (highest-HSC-score cell within that component)
- DPT is computed per component, then concatenated into `obs['dcta_dpt']`
- `obs['dcta_component']` records component ID for each cell

### 5.3 Root robustness

- Compute DPT `K=5` times using the top-K HSC candidates
- Store cell-wise std in `uns['dcta']['dpt_root_std']`
- High std signals a brittle root choice — surfaced in viz and STATE.md

### 5.4 Virtual HSC hook (deferred)

- `pseudotime.run(root_strategy="virtual_hsc_ref", ref_adata=...)` raises
  `NotImplementedError("Phase 2")`. API shape is fixed now so Phase 2 work
  doesn't break callers.

---

## 6. Magic-number / hardcoding policy (rule.md P2)

- **No absolute paths in code.** Dataset path is a CLI arg with optional
  environment variable fallback (`DCTA_DATA_DIR`).
- **No magic numbers buried in stage code.** All defaults live in
  `dcta/defaults.py` and are exposed as CLI flags on `run_phase1.py`.
- **STATE.md registry**: every default value is listed with its justification
  (paper citation, empirical observation, or "scanpy default — needs revisit").
- **Data-driven over fixed**: e.g., `n_dcs` is auto-selected from the
  eigenvalue gap rather than hard-set to a number.

Initial defaults (subject to STATE.md justification):
- `n_top_genes` (HVG) = 3000
- `n_comps_pca` = 50
- `n_comps_diffmap` = 15 (auto-trim by eigvalue gap)
- `knn_diffmap_k` = 30
- `leiden_resolution` = 1.0
- `edge_threshold` (PAGA connectivity) = 0.05
- `robustness_k` (DPT root resampling) = 5

---

## 7. Error handling

Per CLAUDE.md line 79 ("No silent failures"):

- `preprocessing`: empty cells, negative X, double log-transform detection
  → raise `DCTADataError` with diagnostic
- `diffusion`: ambiguous eigenvalue gap → warn + use default, write the gap
  profile to STATE.md (do not silently pick)
- `pseudotime`: no candidate root (HSC marker score all zero, or marker genes
  missing in `var_names`) → raise (no silent fallback)
- `dag`: cycle removal fails to converge → raise + dump residual cycles for
  debug

Error types live in `dcta/__init__.py`: `DCTADataError`, `DCTARootError`,
`DCTAGraphError`.

---

## 8. Testing

### 8.1 Synthetic unit tests (`test_synthetic.py`)

- 3-lineage star-graph synthetic AnnData (HSC → {Mye, Ery, Lym}, 500 cells)
- Each stage tested in isolation: input AnnData → expected mutations
- Includes failure cases for each `DCTA*Error`
- Runs in <30s

### 8.2 Hay smoke test (`test_hay_smoke.py`)

- 1000-cell stratified subset of `hay_hca_balanced_bm.h5ad`
- End-to-end run of the full pipeline
- Asserts: `uns['dcta']['G_data']` exists, is a DAG, has ≥3 nodes, root cluster
  is reachable to ≥1 leaf
- Runs in <2 min

### 8.3 Not in Phase 1

- Benchmark vs SingleR/CellTypist/Azimuth/Palantir/CellRank → Phase 2/3
- Cluster-resolution sweep → Phase 2

---

## 9. Definition of done (Phase 1)

Phase 1 is complete when `results/phase1_hay/` contains:

1. `dcta_phase1.h5ad` — AnnData with all stage outputs
2. `G_data.gml` — NetworkX serialization of the data-driven DAG
3. `diffmap_eigvals.png`, `paga_connectivities.png`, `dpt_distribution.png`,
   `dpt_root_std_hist.png`, `G_data_diagram.png`
4. `STATE.md` updated with: defaults registry, eigvalue-gap report, root
   robustness summary
5. `pytest -v` green for both `test_synthetic.py` and `test_hay_smoke.py`

---

## 10. Honest limitations (carried from CLAUDE.md)

These are NOT to be fixed in Phase 1, only acknowledged:

- Cluster resolution dependency
- Root choice affects all downstream pseudotime (robustness log surfaces this
  but doesn't resolve it)
- Hematopoiesis is not a strict tree; G_data may not be either, hence the FAS
  step
- FAS greedy is one heuristic; alternatives (ILP, randomized) deferred
- Phase 1 is BM-only; blood/spleen require Phase 2 virtual HSC

---

## 11. Addendum (2026-05-21): Edge direction strategy + PAGA↔G_data validation

This addendum revises Sections 4 and 7 in response to the design feedback:
*"DPT 는 biological time point 가 아니므로 CD34/HLF 같은 marker 를 direction
결정에 써야 함; PAGA+DPT tree 와 G_data DAG 의 일치를 검증해야 함."*

### 11.1 Two direction methods computed in parallel (T7 / `dcta.dag`)

For every PAGA edge `(c_i, c_j)` above `edge_threshold`, compute:

- **DPT direction** — original method. Compute cluster-mean DPT `μ_dpt(c)`.
  Orient low→high. (`obs['dcta_dpt']` from T6.)
- **Marker direction** — new. Compute cluster-mean HSC-marker score
  `μ_hsc(c)` (and optionally a maturation score = `μ_mature − μ_hsc` if a
  mature-marker YAML is provided). Orient high→low immaturity (= immature→
  mature). For Phase 1, HSC markers only (CRHBP, AVP, HLF, MEIS1, HOPX, ...
  plus CD34 added per user request). Mature markers are Phase 2 — direction
  is just inverse-immaturity.

Both decisions are stored per edge in `uns['dcta']['direction_methods']`:

```python
{
  "dpt":    {(c_i, c_j): "ci_to_cj" | "cj_to_ci", ...},
  "marker": {(c_i, c_j): "ci_to_cj" | "cj_to_ci", ...},
  "agreement_rate": float,   # fraction of edges where both methods agree
}
```

### 11.2 Accuracy-based selection (when ground-truth labels exist)

If `adata.obs` has a `cell_type` column (Hay dataset) or `true_lineage`
(synthetic fixture), DCTA can score each method against biological seniority:

1. Load `data/priors/lineage_seniority.yaml` — a label→rank dict
   (HSC=0, MPP=1, lineage progenitor=2, mature subset=3, …).
2. For each cluster, take its dominant label (mode of `cell_type` over cells
   in that cluster) and look up seniority. If the label is not in the YAML,
   the cluster contributes NaN and is excluded from accuracy scoring.
3. For each PAGA edge, the **expected** direction = lower-seniority cluster
   → higher-seniority cluster. Ties skipped.
4. **Accuracy** of a method = (# edges where method's direction == expected) /
   (# scoreable edges). Stored as `uns['dcta']['direction_methods']['accuracy']
   = {"dpt": 0.78, "marker": 0.85}`.
5. **Chosen method** = `argmax(accuracy)`. Stored as
   `uns['dcta']['direction_method_chosen']` ∈ {"dpt", "marker"}.
   The chosen method's directions are what build `G_data`.

If no ground-truth column is present, `direction_method_chosen = "marker"` by
default (biology-grounded, embedding-agnostic) and accuracy is not computed.

### 11.3 `data/priors/lineage_seniority.yaml` (new file, Phase 1)

Minimal Phase 1 schema:
```yaml
# Higher rank = more mature.
# Labels not listed are excluded from accuracy scoring (not penalized).
seniority:
  HSC: 0
  MPP: 1
  LMPP: 1
  CMP: 2
  GMP: 2
  MEP: 2
  CLP: 2
  Pro-B: 3
  Pre-B: 3
  Erythroblast: 3
  Granulocyte: 4
  Monocyte: 4
  Mature_B: 4
  Plasma: 5
  # synthetic-fixture aliases:
  Mye: 4
  Ery: 4
  Lym: 4
```

The Phase 1 file is illustrative; T11 STATE.md will record which Hay
`cell_type` labels mapped and which were skipped.

### 11.4 PAGA ↔ G_data validation overlay (merged into T8 `dcta.viz`)

New figure: `paga_vs_dag_overlay.png`. Layout:

- Compute scanpy's PAGA force-directed positions for the cluster nodes
  (`scanpy.tl.paga` already produces a `pos` array when followed by
  `sc.pl.paga(plot=False, layout='fa')`).
- Plot all PAGA edges in light grey (undirected).
- Overlay `G_data` directed edges:
  - Solid green: both `marker` and `dpt` methods agree on direction.
  - Dashed orange: methods disagree (chosen method's direction is drawn).
  - Red ×: edge removed by FAS cycle removal.

New summary fields written to `results/phase1_hay/summary.json`:

```json
"direction": {
  "method_chosen": "marker",
  "accuracy": {"dpt": 0.78, "marker": 0.85},
  "agreement_rate": 0.83,
  "n_edges_scored": 17
}
```

### 11.5 Out-of-scope additions deferred

- A full mature-marker YAML per lineage (Phase 2). Phase 1 uses HSC markers
  only (immaturity → maturity by inverse-score).
- Per-edge confidence scoring beyond agreement_rate (Phase 2).
- Automatic seniority learning from data (e.g., HVG monotonicity); Phase 1
  uses a hand-curated YAML.

---

## 12. Notes for spec reviewer

- **Dataset citation correction (2026-05-21, looked up via cellxgene_census)**:
  the file `hay_hca_balanced_bm.h5ad` in APF's `test_datasets/cxg/` is NOT
  Hay et al. 2018 GSE117498 and NOT Pellin et al. 2019. It is:
    - **Title**: A (Balanced) Bone Marrow Reference Map of Hematopoietic Development
    - **Citation**: Zeng et al. 2025 *Blood Cancer Discovery*
    - **DOI**: 10.1158/2643-3230.BCD-24-0342
    - **CZ Collection**: f6c50495-3361-40ed-a819-fb9644396ed9
    - **dataset_id**: cd2f23c1-aef1-48ae-8eb4-0bcf124e567d
    - **Cell count**: 263,159 (multi-study meta-merge incl. HCA_BM, Setty2019,
      Granja2019, Ainciburu2022, Oetjen2018, Mende2022, ...)
    - **30k local subset**: APF's stratified down-sample.
  Hereafter use the key **Zeng2025_BMM** in code/docs. APF's "hay" key is a
  historical misnomer kept only for file-path back-compat.
- **Phase 2 novelty check (per APF feedback_novelty_preservation)**: Zeng 2025
  BMM is the very paper APF avoids competing with on AML angle. DCTA must
  verify its trajectory-DAG contribution is orthogonal to BMM's reference
  map before publication.
- ~~The CLAUDE.md reference to "Hay et al. 2018 GSE117498 (BM reference)" is
  inconsistent — GSE117498 is Pellin et al. 2019.~~ (Superseded by the
  correction above — neither Hay nor Pellin; the actual source is Zeng 2025.)
- This is a fresh project. APF (`v8_apf_dcbmi/`) is NOT a dependency — DCTA
  may borrow data and may eventually be benchmarked against APF, but the two
  share no code in Phase 1.
