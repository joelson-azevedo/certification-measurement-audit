"""Binary-marginal-preserving, tetrachoric-to-tetrachoric parallel analysis."""

import time

import numpy as np
import pandas as pd

from .tetrachoric import tetrachoric_matrix


def permute_columns(x, rng):
    x = np.asarray(x, float)
    permuted = np.empty_like(x)
    for j in range(x.shape[1]):
        valid = np.isfinite(x[:, j])
        permuted[:, j] = x[:, j]
        # Missing positions, if any, stay fixed; only measured binary values move.
        permuted[valid,j] = rng.permutation(x[valid,j])
    return permuted


def parallel_analysis(x, observed_matrix, item_ids, replicates=100, seed=20261021,
                      percentile=95, config=None, max_attempts=None, progress=None):
    if replicates < 1 or not 0 < percentile < 100:
        raise ValueError("Invalid parallel-analysis replication/percentile request")
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    if not np.isfinite(observed_matrix).all():
        raise ValueError("Observed tetrachoric matrix is undefined")
    observed = np.linalg.eigvalsh(observed_matrix)[::-1]
    references, logs = [], []
    limit = max_attempts or max(replicates+20, int(1.2*replicates))
    start = time.perf_counter()
    for attempt in range(1, limit+1):
        null = permute_columns(x, rng)
        exact_margins = bool(np.array_equal(np.nansum(null, axis=0), np.nansum(x, axis=0))
                             and np.array_equal(np.isnan(null), np.isnan(x)))
        if not exact_margins:
            raise RuntimeError("Reference permutation did not preserve margins/missingness")
        matrix, pairs, _ = tetrachoric_matrix(null, item_ids, config)
        valid = bool(pairs.converged.all() and np.isfinite(matrix).all())
        eigen = np.linalg.eigvalsh(matrix)[::-1] if valid else np.full(x.shape[1], np.nan)
        logs.append({"attempt": attempt, "valid": valid, "margins_exact": exact_margins,
                     "continuity_corrected_pairs": int(pairs.continuity_corrected.sum()),
                     "nonconverged_pairs": int((~pairs.converged).sum()),
                     "reference_minimum_eigenvalue": float(eigen[-1]),
                     "elapsed_seconds": time.perf_counter()-start})
        if valid:
            references.append(eigen)
        if progress and (attempt % 10 == 0 or len(references) == replicates):
            progress(len(references), replicates)
        if len(references) >= replicates:
            break
    if len(references) < replicates:
        raise RuntimeError(f"Only {len(references)} valid tetrachoric references in {limit} attempts")
    refs = np.asarray(references)
    cutoff = np.percentile(refs, percentile, axis=0)
    exceeds = observed > cutoff
    suggested = next((j for j, value in enumerate(exceeds) if not value), len(exceeds))
    table = pd.DataFrame({"eigenvalue_rank": np.arange(1, len(observed)+1), "observed_tetrachoric_eigenvalue": observed,
                          "reference_mean": refs.mean(axis=0), "reference_percentile": cutoff,
                          "percentile_level": percentile, "exceeds_reference_percentile": exceeds,
                          "retained_by_leading_sequence_rule": np.arange(len(observed)) < suggested,
                          "suggested_dimensions": suggested, "valid_replications": len(refs), "seed": seed})
    return table, refs, pd.DataFrame(logs)
