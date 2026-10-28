"""Fixed-calibration, local-independent Q3 reference with full EAP re-estimation."""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
import pandas as pd

from .irt_2pl import response_probability
from .local_dependence import pairwise_q3, residuals_at_eap


@dataclass(frozen=True)
class ReferenceConfig:
    valid_replicates: int = 500
    maximum_attempts: int = 550
    nodes: int = 81
    familywise_quantiles: tuple = (.95, .99)
    quantile_method: str = "higher"

    def __post_init__(self):
        if self.valid_replicates < 2 or self.maximum_attempts < self.valid_replicates or self.nodes < 3:
            raise ValueError("Unsupported simulation configuration")
        if self.quantile_method != "higher" or tuple(self.familywise_quantiles) != (.95, .99):
            raise ValueError("Screening quantiles/method are fixed in this implementation")


def generate_responses(intercept, slope, observed_mask, seed, attempt):
    mask = np.asarray(observed_mask, bool)
    if mask.ndim != 2 or mask.shape[1] != len(slope):
        raise ValueError("Observed missing-mask/item dimension mismatch")
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(attempt)]))
    theta = rng.standard_normal(len(mask))
    p = response_probability(theta, intercept, slope)
    x = (rng.random(mask.shape) < p).astype(float)
    x[~mask] = np.nan
    return x


def simulate_reference(intercept, slope, observed_mask, seed, config=ReferenceConfig(), progress=None):
    start = time.perf_counter()
    values, rows = [], []
    n_pairs = len(slope)*(len(slope)-1)//2
    for attempt in range(config.maximum_attempts):
        t = time.perf_counter()
        row = {"attempt": attempt, "seed_entropy": seed, "seed_spawn_component": attempt,
               "valid_replicate": False, "failure": "", "eap_reestimated": False}
        try:
            x = generate_responses(intercept, slope, observed_mask, seed, attempt)
            if not np.array_equal(np.isfinite(x), observed_mask):
                raise ValueError("Simulated missing mask changed")
            residual, mean, sd = residuals_at_eap(x, intercept, slope, config.nodes)
            row["eap_reestimated"] = True
            if not np.isfinite(mean).all() or not np.isfinite(sd).all():
                raise FloatingPointError("Nonfinite simulated posterior estimate")
            pairs, _, diagnostics = pairwise_q3(residual)
            if diagnostics["n_finite_pairs"] != n_pairs:
                raise FloatingPointError("At least one simulated Q3 pair is unavailable")
            adjusted = pairs.Q3_adjusted.to_numpy()
            values.append(adjusted)
            row.update({"valid_replicate": True, "valid_replicate_number": len(values),
                        "mean_off_diagonal_Q3": diagnostics["mean_off_diagonal_Q3"],
                        "maximum_positive_adjusted_Q3": float(max(0, adjusted.max())),
                        "maximum_raw_Q3": diagnostics["raw_max"], "minimum_raw_Q3": diagnostics["raw_min"],
                        "minimum_pairwise_n": diagnostics["minimum_paired_n"],
                        "mean_simulated_eap": float(mean.mean()), "sd_simulated_eap": float(mean.std(ddof=1))})
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            row["failure"] = str(exc)
        row["elapsed_seconds"] = time.perf_counter()-t
        rows.append(row)
        if progress is not None and (len(rows) % 50 == 0 or len(values) == config.valid_replicates):
            progress(f"{len(values)} valid / {len(rows)} attempted, {time.perf_counter()-start:.1f}s")
        if len(values) == config.valid_replicates:
            break
    return np.asarray(values).reshape(-1, n_pairs), pd.DataFrame(rows)


def reference_thresholds(null_values, observed, config=ReferenceConfig()):
    null_values = np.asarray(null_values, float)
    observed = np.asarray(observed, float)
    if null_values.ndim != 2 or null_values.shape[1] != len(observed) or len(null_values) < 2:
        raise ValueError("Reference/observed pair dimension mismatch")
    if not np.isfinite(null_values).all() or not np.isfinite(observed).all():
        raise ValueError("Cannot calibrate thresholds with nonfinite Q3")
    maxima = np.maximum(0, null_values.max(axis=1))
    q95, q99 = np.quantile(maxima, [.95, .99], method=config.quantile_method)
    n = len(maxima)
    maximum_observed = float(max(0, observed.max()))
    exceedances = int((maxima >= maximum_observed).sum())
    # Monte Carlo uncertainty in the reference tail probability, conditional on
    # the fixed calibration/observed statistic. This is not model uncertainty.
    z = 1.959963984540054
    proportion = exceedances/n
    center = (proportion+z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(proportion*(1-proportion)/n+z*z/(4*n*n))/(1+z*z/n)
    thresholds = {"valid_replicates": n, "number_of_pairs_in_family": len(observed),
                  "observed_maximum_positive_adjusted_Q3": maximum_observed,
                  "familywise_95_threshold": float(q95), "familywise_99_threshold": float(q99),
                  "quantile_method": config.quantile_method,
                  "maximum_reference_exceedances": exceedances,
                  "maximum_monte_carlo_tail_area_plus_one": (1+exceedances)/(n+1),
                  "reference_tail_probability_mc_wilson_95_low": float(max(0, center-half)),
                  "reference_tail_probability_mc_wilson_95_high": float(min(1, center+half)),
                  "monte_carlo_tail_area_resolution": 1/(n+1),
                  "expected_reference_draws_above_95_percentile": .05*n,
                  "expected_reference_draws_above_99_percentile": .01*n,
                  "tail_proportion_standard_error_at_05": float(np.sqrt(.05*.95/n)),
                  "tail_proportion_standard_error_at_01": float(np.sqrt(.01*.99/n)),
                  "reference_scope": "fixed fitted item parameters; new N(0,1) theta; fixed observed missing mask; EAP reestimated; no model refitting"}
    pointwise = pd.DataFrame({"pointwise_adjusted_mean": null_values.mean(axis=0),
                              "pointwise_adjusted_sd": null_values.std(axis=0, ddof=1),
                              "pointwise_adjusted_p01": np.quantile(null_values, .01, axis=0, method="higher"),
                              "pointwise_adjusted_p05": np.quantile(null_values, .05, axis=0, method="higher"),
                              "pointwise_adjusted_median": np.quantile(null_values, .5, axis=0, method="higher"),
                              "pointwise_adjusted_p95": np.quantile(null_values, .95, axis=0, method="higher"),
                              "pointwise_adjusted_p99": np.quantile(null_values, .99, axis=0, method="higher"),
                              "pointwise_reference_exceedances": (null_values >= observed).sum(axis=0),
                              "pointwise_mc_tail_area_plus_one": (1+(null_values >= observed).sum(axis=0))/(n+1),
                              "familywise_mc_tail_area_plus_one": (1+(maxima[:, None] >= observed).sum(axis=0))/(n+1),
                              "reference_replicates": n})
    return thresholds, pointwise


def classify_pairs(observed_pairs, thresholds, pointwise):
    if len(observed_pairs) != len(pointwise):
        raise ValueError("Pairwise reference ordering mismatch")
    result = observed_pairs.reset_index(drop=True).copy()
    for col in pointwise:
        result[col] = pointwise[col].to_numpy()
    result["exceeds_familywise_95"] = result.Q3_adjusted.gt(thresholds["familywise_95_threshold"])
    result["exceeds_familywise_99"] = result.Q3_adjusted.gt(thresholds["familywise_99_threshold"])
    result["exceeds_pointwise_95_exploratory"] = result.Q3_adjusted.gt(result.pointwise_adjusted_p95)
    result["evidence_classification"] = np.select(
        [result.exceeds_familywise_99, result.exceeds_familywise_95,
         result.exceeds_pointwise_95_exploratory & result.Q3_adjusted.gt(0)],
        ["exceeds simulated family-wise 99% threshold", "exceeds simulated family-wise 95% threshold",
         "elevated in exploratory pointwise analysis only"], default="no notable positive residual association")
    result["negative_association_scope"] = "negative values retained; no formal lower-tail or two-sided family screen"
    return result
