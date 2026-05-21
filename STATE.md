# DCTA Phase 1 — State

## Progress
- [x] Scaffold (Task 0)
- [ ] Synthetic fixture (Task 1)
- [ ] io.load_h5ad (Task 2)
- [ ] preprocessing.run (Task 3)
- [ ] diffusion.run (Task 4)
- [ ] connectivity.run (Task 5)
- [ ] pseudotime.run (Task 6)
- [ ] dag.build (Task 7)
- [ ] viz.export_all (Task 8)
- [ ] run_phase1.py (Task 9)
- [ ] test_hay_smoke (Task 10)

## Defaults registry

See `dcta/defaults.py`. Justifications:

| Param | Default | Justification |
|---|---|---|
| MIN_GENES_PER_CELL | 200 | scanpy standard; removes empty droplets |
| MIN_CELLS_PER_GENE | 3 | scanpy standard; removes ultra-rare genes |
| PCT_MITO_MAX | 20.0 | HSPC-permissive (HSCs can have higher mito); standard 10% drops real HSCs |
| N_TOP_GENES | 3000 | scanpy default for moderate-size BM datasets |
| N_COMPS_PCA | 50 | scanpy default |
| N_COMPS_DIFFMAP | 15 | Haghverdi 2016; usually 5-15 captures hematopoietic structure |
| KNN_PCA_K | 15 | scanpy default for sc.pp.neighbors |
| KNN_DIFFMAP_K | 30 | larger than PCA-KNN: diffmap denoises, can afford bigger neighborhood |
| LEIDEN_RESOLUTION | 1.0 | scanpy default; sensitivity tested in Phase 2 |
| EDGE_THRESHOLD | 0.05 | scanpy default for paga edge filter |
| ROBUSTNESS_K | 5 | empirical; Phase 1 budget |
| FAS_METHOD | greedy | only method implemented in Phase 1 |
| TRANSITIVE_REDUCTION | False | keep verified edges, decide in Phase 2 |
| RANDOM_STATE | 0 | reproducibility |

## Run logs

(populated as Phase 1 runs are executed)
