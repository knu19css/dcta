"""DCTA: DAG-Consistent Trajectory Annotation — Phase 1."""

__version__ = "0.1.0-phase1"


class DCTAError(Exception):
    """Base exception for DCTA pipeline errors."""


class DCTADataError(DCTAError):
    """Input AnnData violates an expected invariant."""


class DCTARootError(DCTAError):
    """Cannot identify a valid root for pseudotime."""


class DCTAGraphError(DCTAError):
    """G_data construction failed (e.g., cycle removal did not converge)."""


# Submodule re-exports (added in T2-T5; pseudotime/dag/viz added in T6+)
from dcta import io, preprocessing, diffusion, connectivity, pseudotime  # noqa: E402,F401
