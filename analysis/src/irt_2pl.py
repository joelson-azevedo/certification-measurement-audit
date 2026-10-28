"""Logistic 2PL with genuinely signed slopes, D=1, and standard-normal quadrature."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp


def validate_responses(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or np.isinf(x).any() or not np.isin(x[np.isfinite(x)], [0, 1]).all():
        raise ValueError("Expected a binary response matrix with NaN omissions")
    return x


def normal_quadrature(nodes=31):
    if nodes < 3:
        raise ValueError("At least three Gauss-Hermite nodes required")
    points, weights = np.polynomial.hermite.hermgauss(nodes)
    theta, weights = np.sqrt(2)*points, weights/np.sqrt(np.pi)
    weights /= weights.sum()
    return theta, weights


def response_probability(theta, intercept, slope):
    """Last axis is item; no positivity constraint, no 1.7 scaling."""
    theta = np.asarray(theta, float)
    return expit(theta[..., None]*np.asarray(slope, float)+np.asarray(intercept, float))


def conditional_log_likelihood(x, theta, intercept, slope):
    x = validate_responses(x)
    intercept, slope = np.asarray(intercept, float), np.asarray(slope, float)
    if len(intercept) != x.shape[1] or slope.shape != intercept.shape:
        raise ValueError("Item parameter dimension mismatch")
    observed = np.isfinite(x).astype(float)
    correct = np.nan_to_num(x, nan=0.0)
    eta = np.asarray(theta, float)[:, None]*slope+intercept
    log_p, log_q = -np.logaddexp(0, -eta), -np.logaddexp(0, eta)
    return correct@log_p.T+(observed-correct)@log_q.T


def marginal_components(x, intercept, slope, nodes=31):
    theta, weights = normal_quadrature(nodes)
    conditional = conditional_log_likelihood(x, theta, intercept, slope)
    joint = conditional+np.log(weights)[None, :]
    log_marginal = logsumexp(joint, axis=1)
    posterior = np.exp(joint-log_marginal[:, None])
    return theta, weights, log_marginal, posterior


def marginal_log_likelihood(x, intercept, slope, nodes=31):
    return float(marginal_components(x, intercept, slope, nodes)[2].sum())


def expected_score(theta, intercept, slope):
    return response_probability(theta, intercept, slope).sum(axis=-1)


def item_information(theta, intercept, slope):
    p = response_probability(theta, intercept, slope)
    return np.asarray(slope, float)**2*p*(1-p)


def test_information(theta, intercept, slope):
    return item_information(theta, intercept, slope).sum(axis=-1)


def eap(x, intercept, slope, nodes=31):
    theta, _, _, posterior = marginal_components(x, intercept, slope, nodes)
    mean = posterior@theta
    variance = posterior@(theta*theta)-mean*mean
    return mean, np.sqrt(np.maximum(variance, 0)), posterior


def orient_parameters(intercept, slope):
    """Resolve global sign indeterminacy by positive aggregate discrimination.

    This flips the global latent orientation only, never responses or an individual
    item's sign relative to the dominant item pattern.
    """
    intercept, slope = np.asarray(intercept, float).copy(), np.asarray(slope, float).copy()
    flipped = bool(slope.sum() < 0)
    return intercept, -slope if flipped else slope, flipped


def interpretable_difficulty(intercept, slope, minimum_abs_slope=.10, slope_se=None):
    intercept, slope = np.asarray(intercept, float), np.asarray(slope, float)
    b = np.full_like(slope, np.nan)
    reason = np.full(slope.shape, "defined: P(correct)=0.5 location", dtype=object)
    near_zero = np.abs(slope) < minimum_abs_slope
    reason[near_zero] = "undefined/unstable: absolute slope below declared threshold"
    supported = ~near_zero
    if slope_se is not None:
        se = np.asarray(slope_se, float)
        unavailable = ~np.isfinite(se)
        crosses_zero = np.abs(slope) <= 1.959963984540054*se
        reason[~near_zero & unavailable] = "undefined/unstable: slope uncertainty unavailable"
        reason[~near_zero & ~unavailable & crosses_zero] = "undefined/unstable: slope interval includes zero"
        supported &= ~unavailable & ~crosses_zero
    b[supported] = -intercept[supported]/slope[supported]
    return b, reason


def approximate_cut_location(intercept, slope, target=42.0, lower=-4, upper=4, points=321):
    """Expected-score mapping only; allow detection of multiple/nonmonotone roots."""
    grid = np.linspace(lower, upper, points)
    scores = expected_score(grid, intercept, slope)
    delta = scores-target
    roots = []
    for j in range(len(grid)-1):
        if delta[j] == 0:
            roots.append(float(grid[j]))
        elif delta[j]*delta[j+1] < 0:
            roots.append(float(brentq(lambda t: float(expected_score(t,intercept,slope))-target, grid[j],grid[j+1],xtol=1e-12)))
    if delta[-1] == 0:
        roots.append(float(grid[-1]))
    roots = sorted(set(round(r,12) for r in roots))
    unique = len(roots)==1
    theta = roots[0] if unique else np.nan
    information = float(test_information(theta,intercept,slope)) if unique else np.nan
    return {"target_expected_full_form_score": target, "theta_at_expected_score":theta,
            "n_roots_in_prespecified_range":len(roots), "all_roots":roots,
            "tcc_monotone_on_grid":bool(np.all(np.diff(scores)>=-1e-10)),
            "information_at_location":information, "conditional_sem_at_location":1/np.sqrt(information) if information>0 else np.nan,
            "mapping_status":"unique expected-score location; not an operational theta cut" if unique else "no unique mapping within prespecified range"}
