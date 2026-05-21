# DCTA: DAG-Consistent Trajectory Annotation

## Project goal

Develop a novel automated cell type annotation method for scRNA-seq hematopoietic
data (bone marrow, blood, spleen). The framework reconciles a data-driven
trajectory graph with a prior hematopoietic DAG, anchored at HSCs and terminal
mature cell types, using YAML-encoded marker rules for biological consistency.

Two target outputs:
1. Methods paper: the framework itself
2. Clinical paper: AML/MDS deviation detection

## Core framework (4 building blocks)

1. **Diffusion map** — unsupervised manifold, no lineage-gene weighting
   (avoids circular reasoning; preserves abnormal programs for deviation detection)
2. **PAGA** — cluster-level connectivity from KNN graph (built on diffusion components)
3. **DPT** — root-based pseudotime for directionality
   (root = real HSC in BM, virtual HSC centroid in blood/spleen)
4. **Prior hematopoietic DAG** — biological anchor with YAML marker rules

## Pipeline

```
preprocessing (cell cycle regress out)
  -> PCA -> KNN
  -> diffusion map (n_comps ~ 10-15, check eigenvalue gap)
  -> diffusion-based KNN
  -> Leiden clustering
  -> PAGA (connectivities)
  -> root assignment (real or virtual HSC)
  -> DPT
  -> edge threshold + DPT-based direction
  -> cycle removal (FAS greedy, weakest edge in cycle)
  -> transitive reduction (optional)
  => G_data (data-driven DAG)

G_data + G_prior + YAML rules
  -> per-cluster marker score (initial assignment to prior nodes)
  -> structural consistency check (do directed edges in G_data
     match the direction in G_prior?)
  -> absorbing Markov chain on G_prior with transition probabilities
     informed by G_data (PAGA connectivity * DPT gradient)
  -> absorption probability => terminal mature annotation
  -> transient state probability => progenitor annotation
  -> YAML rule filtering for biological plausibility
  => final annotation + confidence + deviation flag
```

## Key design decisions (DO NOT change without discussion)

- **No lineage-gene weighting on the manifold**: keep diffusion map unsupervised
  so abnormal (AML/MDS) programs are preserved
- **Cell cycle regress out is mandatory** before diffusion map
- **Root is the only mandatory prior for direction** (everything else flows from
  topology + DPT)
- **Virtual HSC for blood/spleen**: use BM reference HSC centroid mapped via
  Symphony/scArches, not synthetic marker-based centroid
- **Cluster resolution sensitivity** is a known issue — every result needs
  a clustree-style robustness check

## Stack

- Python primary (this package), R for some plots (ggalluvial etc., separate scripts)
- Conda env: `sc` (Python 3.10.19)
- Key deps: scanpy, anndata, networkx, pyyaml, scikit-learn, numpy, scipy,
  matplotlib, harmonypy (or scvi-tools), scrublet
- For R interop: anndata2ri only if needed

## Coding conventions

- All public functions take `AnnData` and return `AnnData` (in-place modifications
  with `copy=False` option, scanpy style)
- Add new fields to `adata.uns['dcta']` namespace to avoid collisions
- Type hints on all public functions
- Docstrings: numpy style
- No silent failures — raise informative errors when prior DAG and data
  disagree in unexpected ways

## Testing

- Unit tests with pytest, on small synthetic AnnData
- Integration test on one public BM dataset (e.g., Hay et al. 2018 GSE117498)
- Benchmark vs SingleR, CellTypist, Azimuth (label-level)
- Benchmark vs Palantir, CellRank (trajectory-level)

## Honest limitations (acknowledge in paper)

- Cluster resolution dependency
- Root choice affects all downstream pseudotime
- Hematopoiesis is not a strict tree (lineage priming, plasticity) — prior DAG
  is a simplification
- FAS greedy for cycle removal is one heuristic among many
- Blood/spleen has limited progenitor representation; prior dependence is higher

## References to revisit when implementing

- Haghverdi et al. 2016, Nat Methods (DPT)
- Wolf et al. 2019, Genome Biol (PAGA)
- Setty et al. 2019, Nat Biotechnol (Palantir)
- Lange et al. 2022, Nat Methods (CellRank)
- Pellin et al. 2019, Nat Commun (hematopoietic trajectory, SPRING + diffusion)
- Hay et al. 2018, Exp Hematol (BM reference)

## Current phase

Phase 1: proof of concept on a single BM dataset, myeloid lineage only.
