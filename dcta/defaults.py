"""Default parameters for DCTA Phase 1. Every default is justified in STATE.md."""

# Preprocessing
MIN_GENES_PER_CELL = 200          # standard scanpy QC floor
MIN_CELLS_PER_GENE = 3            # standard scanpy QC floor
PCT_MITO_MAX = 20.0               # percent
N_TOP_GENES = 3000                # HVG count

# Diffusion / PCA / KNN
N_COMPS_PCA = 50
N_COMPS_DIFFMAP = 15              # max; auto-trim by eigvalue gap
KNN_PCA_K = 15                    # for initial KNN before diffusion
KNN_DIFFMAP_K = 30                # for diffusion-based KNN

# Connectivity
LEIDEN_RESOLUTION = 1.0

# PAGA edge filter
EDGE_THRESHOLD = 0.05             # PAGA connectivity threshold

# Pseudotime
ROBUSTNESS_K = 5                  # top-K HSC candidates for DPT root resampling

# DAG
FAS_METHOD = "greedy"             # feedback arc set heuristic
TRANSITIVE_REDUCTION = False      # off by default; keep all confirmed edges

# Reproducibility
RANDOM_STATE = 0
