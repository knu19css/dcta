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
| `pseudotime` | `obs['dcta_dpt']`, `obs['dcta_component']`, `uns['dcta']['root_strategy']`, `uns['dcta']['root_cells']`, `uns['dcta']['root_cluster']`, `uns['dcta']['component_count']`, `uns['dcta']['dpt_root_std']` |
| `dag` | `uns['dcta']['G_data']` (edge list with weight, direction, retained flag), `uns['dcta']['cycle_removal_log']` (list of edges removed and why) |

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
                       transitive_reduction=False)

adata.write_h5ad("results/phase1_hay/dcta_phase1.h5ad")
dcta.viz.export_all(adata, outdir="results/phase1_hay/")
```

All knobs above are CLI flags on `run_phase1.py` with defaults pulled from
`dcta/defaults.py`.

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

## 11. Notes for spec reviewer

- The CLAUDE.md reference to "Hay et al. 2018 GSE117498 (BM reference)" is
  inconsistent — GSE117498 is Pellin et al. 2019 (FACS-sorted HSPCs). APF
  stores it under the dataset key `hay_hca_balanced_bm`. The plan uses the
  APF naming for the file but should clarify the citation in Phase 1
  README. (Resolve before Phase 2 paper writing.)
- This is a fresh project. APF (`v8_apf_dcbmi/`) is NOT a dependency — DCTA
  may borrow data and may eventually be benchmarked against APF, but the two
  share no code in Phase 1.
