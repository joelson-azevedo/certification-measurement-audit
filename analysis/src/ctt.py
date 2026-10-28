"""Exploratory CTT, with explicit denominators and sample-variance conventions."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CTTConfig:
    very_difficult_below: float = 0.20
    very_easy_above: float = 0.90
    low_discrimination_below: float = 0.20
    attention_discrimination_below: float = 0.30
    bootstrap_replicates: int = 2000
    bootstrap_seed: int = 20261021
    alpha_change_tolerance: float = 1e-12
    material_alpha_change: float = 0.005

    def __post_init__(self):
        if not 0 <= self.very_difficult_below < self.very_easy_above <= 1:
            raise ValueError("Invalid difficulty thresholds")
        if not 0 <= self.low_discrimination_below < self.attention_discrimination_below <= 1:
            raise ValueError("Invalid discrimination thresholds")
        if self.bootstrap_replicates < 0 or self.bootstrap_seed < 0:
            raise ValueError("Bootstrap replicates/seed must be nonnegative")

    def as_dict(self):
        return asdict(self)


def binary_matrix(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all() or not np.isin(x, [0, 1]).all():
        raise ValueError("Expected a complete two-dimensional binary score matrix")
    return x


def wilson_interval(correct, administered, z=1.959963984540054):
    if administered < 0 or correct < 0 or correct > administered:
        raise ValueError("Invalid binomial counts")
    if administered == 0:
        return np.nan, np.nan
    p = correct / administered
    scale = 1 + z*z/administered
    center = (p + z*z/(2*administered)) / scale
    half = z * np.sqrt(p*(1-p)/administered + z*z/(4*administered**2)) / scale
    # Wilson endpoints are exactly 0/1 at the corresponding count extremes.
    return (0.0 if correct == 0 else max(0.0, center-half),
            1.0 if correct == administered else min(1.0, center+half))


def corrected_item_total(x):
    """Pearson correlation of X_j with sum(X)-X_j; NaN is never recoded to zero."""
    x = binary_matrix(x)
    n, k = x.shape
    result = np.full(k, np.nan)
    statuses = []
    rest = x.sum(axis=1, keepdims=True) - x
    for j in range(k):
        if n < 2:
            statuses.append("undefined: fewer than two candidates")
            continue
        a, b = x[:, j] - x[:, j].mean(), rest[:, j] - rest[:, j].mean()
        aa, bb = a @ a, b @ b
        if aa == 0:
            statuses.append("undefined: constant item score")
        elif bb == 0:
            statuses.append("undefined: zero rest-score variance")
        else:
            result[j] = (a @ b) / np.sqrt(aa*bb)
            statuses.append("defined")
    return result, statuses


def cronbach_alpha(x):
    x = binary_matrix(x)
    n, k = x.shape
    if n < 2 or k < 2:
        return np.nan
    total_variance = x.sum(axis=1).var(ddof=1)
    if total_variance == 0:
        return np.nan
    return float(k/(k-1) * (1-x.var(axis=0, ddof=1).sum()/total_variance))


def kr20(x):
    """pq adjusted by n/(n-1) to share alpha's ddof=1 denominator convention."""
    x = binary_matrix(x)
    n, k = x.shape
    if n < 2 or k < 2:
        return np.nan
    total_variance = x.sum(axis=1).var(ddof=1)
    if total_variance == 0:
        return np.nan
    p = x.mean(axis=0)
    sample_item_variance_sum = n/(n-1) * np.sum(p*(1-p))
    return float(k/(k-1) * (1-sample_item_variance_sum/total_variance))


def alpha_if_deleted(x):
    x = binary_matrix(x)
    n, k = x.shape
    if n < 2 or k < 3:
        return np.full(k, np.nan)
    variances = x.var(axis=0, ddof=1)
    rest_variance = (x.sum(axis=1, keepdims=True)-x).var(axis=0, ddof=1)
    result = np.full(k, np.nan)
    valid = rest_variance > 0
    result[valid] = (k-1)/(k-2) * (1-(variances.sum()-variances[valid])/rest_variance[valid])
    return result


def bootstrap_alpha(x, replicates=2000, seed=20261021):
    """Resample independent candidate rows within form; percentile 95% interval."""
    x = binary_matrix(x)
    n, k = x.shape
    if replicates < 0:
        raise ValueError("Replicates must be nonnegative")
    draws = np.full(replicates, np.nan)
    if n >= 2 and k >= 2 and replicates:
        rng = np.random.default_rng(seed)
        totals = x.sum(axis=1)
        for b in range(replicates):
            idx = rng.integers(0, n, size=n)
            variance = totals[idx].var(ddof=1)
            if variance > 0:
                p = x[idx].mean(axis=0)
                draws[b] = k/(k-1) * (1-(n/(n-1)*np.sum(p*(1-p)))/variance)
    valid = draws[np.isfinite(draws)]
    interval = np.quantile(valid, [0.025, 0.975]) if len(valid) else [np.nan, np.nan]
    return {"alpha_ci_low": float(interval[0]), "alpha_ci_high": float(interval[1]),
            "bootstrap_replicates": replicates, "bootstrap_valid_replicates": len(valid),
            "bootstrap_undefined_replicates": int(replicates-len(valid)), "bootstrap_seed": seed}, draws


def select_population(admin, long, first_only=True):
    """Explicit Boolean filter; never use truthiness of strings such as 'False'."""
    if admin.is_first_attempt.dtype != bool or long.is_first_attempt.dtype != bool:
        raise ValueError("First-attempt fields must be parsed as actual Booleans")
    selected = admin.loc[admin.is_first_attempt].copy() if first_only else admin.copy()
    if first_only and selected.candidate_id.duplicated().any():
        raise ValueError("Primary population contains repeated candidates")
    responses = long.loc[long.administration_id.isin(selected.administration_id)].copy()
    return selected, responses


def form_matrices(admin, long, metadata, form):
    members = metadata.loc[metadata.form.isin([form, "AB"]), "item_id"].tolist()
    sittings = admin.loc[admin.form.eq(form)].sort_values("administration_id")
    subset = long.loc[long.form.eq(form)]
    if subset.duplicated(["administration_id", "item_id"]).any():
        raise ValueError("Duplicated administration-item pair")
    if set(subset.item_id) != set(members) or set(subset.administration_id) != set(sittings.administration_id):
        raise ValueError("Inconsistent item-to-form membership or administration coverage")
    scores = subset.pivot(index="administration_id", columns="item_id", values="item_score").reindex(index=sittings.administration_id, columns=members)
    responses = subset.pivot(index="administration_id", columns="item_id", values="response").reindex(index=sittings.administration_id, columns=members)
    if scores.isna().any().any() or responses.isna().any().any():
        raise ValueError("Incomplete administered-item mapping")
    x = binary_matrix(scores.to_numpy())
    return sittings, members, x, responses.to_numpy(dtype=str)


def form_ctt(admin, long, metadata, form, config=None, uncertainty=True, seed=None):
    config = config or CTTConfig()
    sittings, members, x, responses = form_matrices(admin, long, metadata, form)
    if uncertainty and sittings.candidate_id.duplicated().any():
        raise ValueError("Uncertainty estimation requires independent candidates; exclude retakes")
    n, k = x.shape
    if n == 0:
        raise ValueError("No candidates in requested form")
    r, r_status = corrected_item_total(x)
    alpha, kr = cronbach_alpha(x), kr20(x)
    if np.isfinite(alpha) and not np.isclose(alpha, kr, rtol=0, atol=1e-12):
        raise ValueError("KR-20/alpha equivalence check failed")
    deleted = alpha_if_deleted(x)
    item = metadata.set_index("item_id").loc[members, ["domain", "anchor_item"]].reset_index()
    item.insert(0, "form", form)
    item["n_administered"] = n
    item["n_correct"] = x.sum(axis=0).astype(int)
    item["n_omitted"] = (responses == "").sum(axis=0)
    item["proportion_correct"] = x.mean(axis=0)
    ci = [wilson_interval(int(v), n) if uncertainty else (np.nan, np.nan) for v in item.n_correct]
    item["p_wilson_95_low"], item["p_wilson_95_high"] = np.array(ci).T
    item["p_interval_status"] = "Wilson; independent first candidates" if uncertainty else "not estimated: retakes are dependent"
    item["corrected_item_total_r"] = r
    item["discrimination_status"] = r_status
    item["very_easy"] = item.proportion_correct.gt(config.very_easy_above)
    item["very_difficult"] = item.proportion_correct.lt(config.very_difficult_below)
    item["negative_discrimination"] = item.corrected_item_total_r.lt(0)
    item["low_discrimination"] = item.corrected_item_total_r.ge(0) & item.corrected_item_total_r.lt(config.low_discrimination_below)
    item["attention_discrimination"] = item.corrected_item_total_r.ge(config.low_discrimination_below) & item.corrected_item_total_r.lt(config.attention_discrimination_below)
    item["undefined_discrimination"] = item.corrected_item_total_r.isna()
    item["discrimination_screen"] = np.select(
        [item.undefined_discrimination, item.negative_discrimination, item.low_discrimination, item.attention_discrimination],
        ["undefined", "negative: critical investigation", "low", "attention"], default="generally acceptable screening")
    deletion = item[["form", "item_id", "domain", "anchor_item"]].copy()
    deletion["full_test_alpha"] = alpha
    deletion["alpha_if_item_deleted"] = deleted
    deletion["alpha_change_if_deleted"] = deleted-alpha
    deletion["alpha_increases_if_deleted"] = deletion.alpha_change_if_deleted.gt(config.alpha_change_tolerance)
    deletion["deleted_alpha_status"] = np.where(np.isfinite(deleted), "defined", "undefined: insufficient items/candidates or zero rest-score variance")
    totals = x.sum(axis=1)
    stats = {"form": form, "n_administrations": n, "n_unique_candidates": int(sittings.candidate_id.nunique()),
             "n_items": k, "mean_raw_score": float(totals.mean()),
             "sd_raw_score": float(totals.std(ddof=1)) if n>1 else np.nan,
             "score_variance": float(totals.var(ddof=1)) if n>1 else np.nan,
             "score_min": int(totals.min()), "score_max": int(totals.max()),
             "kr20": kr, "cronbach_alpha": alpha, "kr20_minus_alpha": kr-alpha,
             "reliability_status": "defined" if np.isfinite(alpha) else "undefined: insufficient observations/items or zero total-score variance",
             "variance_convention": "ddof=1; KR-20 pq multiplied by n/(n-1)"}
    if uncertainty:
        boot, draws = bootstrap_alpha(x, config.bootstrap_replicates, config.bootstrap_seed if seed is None else seed)
    else:
        boot = {"alpha_ci_low": np.nan, "alpha_ci_high": np.nan, "bootstrap_replicates": 0,
                "bootstrap_valid_replicates": 0, "bootstrap_undefined_replicates": 0, "bootstrap_seed": ""}
        draws = np.array([])
    stats.update(boot)
    stats["uncertainty_scope"] = "candidate bootstrap; first attempts only" if uncertainty else "descriptive only; no uncertainty treating retakes as independent"
    distribution = pd.Series(totals.astype(int)).value_counts().sort_index().rename_axis("raw_score").reset_index(name="n")
    distribution.insert(0, "form", form)
    distribution["proportion"] = distribution.n/n
    return stats, item, deletion, distribution, draws
