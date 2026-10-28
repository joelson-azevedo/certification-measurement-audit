"""Pairwise Q3 and explicitly descriptive leave-two-out follow-up checks."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit

from .irt_2pl import eap, marginal_components, response_probability, validate_responses


def residuals_at_eap(x, intercept, slope, nodes=81):
    """Use fixed item calibration; re-estimate each person's full-response EAP."""
    x = validate_responses(x)
    mean, sd, _ = eap(x, intercept, slope, nodes)
    probability = response_probability(mean, intercept, slope)
    residual = np.where(np.isfinite(x), x-probability, np.nan)
    return residual, mean, sd


def pair_indices(n_items):
    if n_items < 2:
        raise ValueError("At least two administered items required")
    return np.triu_indices(n_items, 1)


def pairwise_q3(residual):
    """Pearson correlations with pair-specific centering and observed denominators.

    Each unordered pair has at least two observed values and positive variances.
    Near-zero centered sums (<=1e-12*max(1, uncentered sum of squares)) are
    unavailable, not zero correlations. The mean is unweighted over finite pairs.
    """
    residual = np.asarray(residual, float)
    if residual.ndim != 2 or np.isinf(residual).any():
        raise ValueError("Residuals must be a matrix containing finite values/NaN")
    observed = np.isfinite(residual).astype(float)
    filled = np.nan_to_num(residual, nan=0)
    count = observed.T@observed
    sums = filled.T@observed
    sums2 = (filled*filled).T@observed
    cross = filled.T@filled
    safe_count = np.maximum(count, 1)
    var1 = sums2-sums*sums/safe_count
    var2 = var1.T
    covariance = cross-sums*sums.T/safe_count
    valid = ((count >= 2) & (var1 > 1e-12*np.maximum(1, sums2))
             & (var2 > 1e-12*np.maximum(1, sums2.T)))
    denominator = np.sqrt(np.maximum(var1, 0)*np.maximum(var2, 0))
    matrix = np.divide(covariance, denominator, out=np.full_like(covariance, np.nan), where=valid)
    if np.any(np.abs(matrix[np.isfinite(matrix)]) > 1+1e-9):
        raise FloatingPointError("Residual correlation outside numerical range")
    matrix = np.clip(matrix, -1, 1)
    i, j = pair_indices(residual.shape[1])
    raw = matrix[i, j]
    mean = float(np.nanmean(raw)) if np.isfinite(raw).any() else np.nan
    status = np.where(count[i, j] < 2, "unavailable: fewer than two paired responses",
                      np.where(valid[i, j], "available", "unavailable: constant/negligible residual variance"))
    table = pd.DataFrame({"item_index_1": i, "item_index_2": j,
                          "n_paired": count[i, j].astype(int), "Q3_raw": raw,
                          "Q3_adjusted": raw-mean, "pair_status": status})
    np.fill_diagonal(matrix, np.nan)
    finite = raw[np.isfinite(raw)]
    diagnostics = {"n_candidates": len(residual), "n_items": residual.shape[1],
                   "n_expected_pairs": len(i), "n_finite_pairs": len(finite),
                   "n_unavailable_pairs": int((~np.isfinite(raw)).sum()),
                   "mean_off_diagonal_Q3": mean,
                   "mean_adjusted_Q3": float(np.nanmean(raw-mean)) if len(finite) else np.nan,
                   "minimum_paired_n": int(count[i, j].min()), "maximum_paired_n": int(count[i, j].max())}
    for label, values in [("raw", finite), ("adjusted", finite-mean)]:
        diagnostics.update({f"{label}_{name}": float(value) for name, value in zip(
            ["min", "p01", "p05", "median", "p95", "p99", "max"],
            np.quantile(values, [0, .01, .05, .5, .95, .99, 1]) if len(values) else [np.nan]*7)})
        diagnostics[f"{label}_sd_across_pairs"] = float(np.std(values, ddof=1)) if len(values)>1 else np.nan
    return table, matrix, diagnostics


def annotate_pairs(table, members, metadata):
    if len(set(members)) != len(members):
        raise ValueError("Duplicate administered item identifiers")
    meta = metadata.set_index("item_id")
    if not meta.index.is_unique:
        raise ValueError("Duplicate item metadata")
    result = table.copy()
    for side in [1, 2]:
        ids = np.asarray(members)[result[f"item_index_{side}"].to_numpy()]
        result[f"item_{side}"] = ids
        result[f"domain_{side}"] = meta.loc[ids, "domain"].to_numpy()
        anchors = meta.loc[ids, "anchor_item"].astype(str).str.lower().map({"true": True, "false": False})
        if anchors.isna().any():
            raise ValueError("Unrecognized anchor designation")
        result[f"anchor_{side}"] = anchors.to_numpy(dtype=bool)
    result["same_domain"] = result.domain_1.eq(result.domain_2)
    return result


def rank_pairs(table):
    result = table.sort_values(["Q3_adjusted", "item_1", "item_2"], ascending=[False, True, True], na_position="last").copy()
    result["positive_screen_rank"] = np.arange(1, len(result)+1)
    return result


def leave_two_out(x, intercept, slope, j, k, nodes=81):
    """Integrate joint Bernoulli cells over the posterior from OTHER items only.

    Expected counts are sum_i int P_j(u)^v(1-P_j(u))^(1-v)
    P_k(u)^w(1-P_k(u))^(1-w) p(u|x_i,-j,-k)du, on jointly observed i.
    This is neither a Chen-Thissen statistic nor an independent-ability check.
    """
    x = validate_responses(x)
    if j == k or min(j, k) < 0 or max(j, k) >= x.shape[1]:
        raise ValueError("Two distinct valid focal item indices required")
    rest = x.copy(); rest[:, [j, k]] = np.nan
    theta, prior, _, posterior = marginal_components(rest, intercept, slope, nodes)
    probabilities = response_probability(theta, intercept, slope)
    pj, pk = probabilities[:, j], probabilities[:, k]
    cells = np.column_stack([(1-pj)*(1-pk), (1-pj)*pk, pj*(1-pk), pj*pk])
    valid = np.isfinite(x[:, j]) & np.isfinite(x[:, k])
    observed = np.bincount((2*x[valid, j]+x[valid, k]).astype(int), minlength=4)
    expected = (posterior[valid]@cells).sum(axis=0)
    mean = posterior@theta
    rest_p = posterior@probabilities[:, [j, k]]
    point_p = expit(np.asarray(intercept)[[j, k]]+mean[:, None]*np.asarray(slope)[[j, k]])
    residual_point = np.where(np.isfinite(x[:, [j, k]]), x[:, [j, k]]-point_p, np.nan)
    residual_integrated = np.where(np.isfinite(x[:, [j, k]]), x[:, [j, k]]-rest_p, np.nan)
    point_q3 = pairwise_q3(residual_point)[0].Q3_raw.iloc[0]
    integrated_q3 = pairwise_q3(residual_integrated)[0].Q3_raw.iloc[0]
    return {"n_paired": int(valid.sum()), "observed_cells": observed,
            "expected_cells_given_rest": expected, "prior_expected_cells": valid.sum()*(prior@cells),
            "Q3_at_leave_two_out_eap": float(point_q3),
            "residual_correlation_at_rest_posterior_expected_probability": float(integrated_q3),
            "mean_rest_eap": float(mean[valid].mean()) if valid.any() else np.nan,
            "posterior": posterior, "rest_eap": mean,
            "scope": "descriptive; focal likelihoods excluded; fixed fitted parameters reuse all calibration data; no chi-square/p-value"}


def domain_summary(pairs):
    rows = []
    for same, group in pairs.groupby("same_domain", sort=True):
        rows.append({"grouping": "same_domain_vs_between_domain", "category": "within domain" if same else "between domains",
                     "n_pairs": len(group), "mean_adjusted_q3": group.Q3_adjusted.mean(),
                     "median_adjusted_q3": group.Q3_adjusted.median(), "maximum_adjusted_q3": group.Q3_adjusted.max(),
                     "n_above_familywise_95": int(group.exceeds_familywise_95.sum()),
                     "n_above_familywise_99": int(group.exceeds_familywise_99.sum()),
                     "familywise_95_flag_rate": group.exceeds_familywise_95.mean(),
                     "n_pointwise_only": int(group.evidence_classification.eq("elevated in exploratory pointwise analysis only").sum()),
                     "pointwise_only_rate": group.evidence_classification.eq("elevated in exploratory pointwise analysis only").mean(),
                     "proportion_of_all_pairs": len(group)/len(pairs),
                     "scope": "descriptive proportions; overlapping item pairs are not independent; domain labels are not causes"})
    domain_keys = pairs.apply(lambda row: " / ".join(sorted([row.domain_1, row.domain_2])), axis=1)
    for name, group in pairs.groupby(domain_keys, sort=True):
        rows.append({"grouping": "domain_combination", "category": name, "n_pairs": len(group),
                     "mean_adjusted_q3": group.Q3_adjusted.mean(), "median_adjusted_q3": group.Q3_adjusted.median(),
                     "maximum_adjusted_q3": group.Q3_adjusted.max(),
                     "n_above_familywise_95": int(group.exceeds_familywise_95.sum()),
                     "n_above_familywise_99": int(group.exceeds_familywise_99.sum()),
                     "familywise_95_flag_rate": group.exceeds_familywise_95.mean(), "proportion_of_all_pairs": len(group)/len(pairs),
                     "n_pointwise_only": int(group.evidence_classification.eq("elevated in exploratory pointwise analysis only").sum()),
                     "pointwise_only_rate": group.evidence_classification.eq("elevated in exploratory pointwise analysis only").mean(),
                     "scope": "descriptive; item text and stimulus relationships unavailable"})
    return pd.DataFrame(rows)
