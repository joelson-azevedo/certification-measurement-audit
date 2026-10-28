"""Deterministic marginal-threshold tetrachorics; no presumed library API."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri


@dataclass(frozen=True)
class TetrachoricConfig:
    continuity_correction: float = 0.5
    sparse_cell_below: int = 5
    minimum_pair_n: int = 20
    probability_tolerance: float = 1e-10
    quadrature_tolerance: float = 1e-8
    quadrature_nodes: int = 32
    max_iterations: int = 45


@lru_cache(maxsize=4)
def quadrature(nodes):
    return np.polynomial.legendre.leggauss(nodes)


def bivariate_normal_cdf(a, b, rho, nodes=32):
    """Plackett integral, t=sin(theta), fixed Gauss-Legendre quadrature.

    Phi2(a,b,r)=Phi(a)Phi(b) + integral_0^asin(r)
    exp(-(a²-2ab sin(theta)+b²)/(2 cos²(theta))) dtheta/(2pi).
    No randomized bivariate-normal CDF is used.
    """
    a, b, rho = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float), np.asarray(rho, float))
    if np.any(np.abs(rho) >= 1):
        raise ValueError("Correlation must lie strictly between -1 and 1")
    shape = a.shape
    a, b, rho = a.ravel(), b.ravel(), rho.ravel()
    points, weights = quadrature(nodes)
    angle = np.arcsin(rho)
    theta = angle[:, None]*(points[None, :]+1)/2
    s = np.sin(theta)
    exponent = -(a[:, None]**2-2*a[:, None]*b[:, None]*s+b[:, None]**2)/(2*(1-s*s))
    integrals = angle/2*np.sum(np.exp(exponent)*weights[None, :], axis=1)/(2*np.pi)
    return (ndtr(a)*ndtr(b)+integrals).reshape(shape)


def adaptive_cdf(a, b, rho):
    angle = np.arcsin(rho)
    def integrand(theta):
        s = np.sin(theta)
        return np.exp(-(a*a-2*a*b*s+b*b)/(2*(1-s*s)))/(2*np.pi)
    integral, _ = quad(integrand, 0, angle, epsabs=1e-12, epsrel=1e-12, limit=150)
    return float(ndtr(a)*ndtr(b)+integral)


def estimate_tables(tables, config=None):
    """Cells in order n00,n01,n10,n11; thresholds estimated from table marginals."""
    config = config or TetrachoricConfig()
    cells = np.atleast_2d(np.asarray(tables, float))
    if cells.shape[1] != 4 or not np.isfinite(cells).all() or (cells < 0).any():
        raise ValueError("Expected finite nonnegative 2x2 table counts")
    n = cells.sum(axis=1)
    nonconstant = ((cells[:, :2].sum(axis=1)>0) & (cells[:, 2:].sum(axis=1)>0)
                   & (cells[:, [0, 2]].sum(axis=1)>0) & (cells[:, [1, 3]].sum(axis=1)>0))
    valid = nonconstant & (n >= config.minimum_pair_n)
    corrected = np.any(cells == 0, axis=1) & valid
    adjusted = cells + corrected[:, None]*config.continuity_correction
    rho = np.full(len(cells), np.nan)
    iterations = np.zeros(len(cells), dtype=int)
    errors = np.full(len(cells), np.nan)
    quadrature_errors = np.full(len(cells), np.nan)
    converged = np.zeros(len(cells), dtype=bool)
    methods = np.full(len(cells), "not estimated", dtype=object)
    statuses = np.where(n < config.minimum_pair_n, "insufficient pair N", "constant binary margin").astype(object)
    idx = np.flatnonzero(valid)
    if len(idx):
        t = adjusted[idx]
        denom = t.sum(axis=1)
        pa, pb, target = (t[:,0]+t[:,1])/denom, (t[:,0]+t[:,2])/denom, t[:,0]/denom
        a, b = ndtri(pa), ndtri(pb)
        phi = (target-pa*pb)/np.sqrt(pa*(1-pa)*pb*(1-pb))
        r = np.clip(np.sin(np.pi/2*phi), -.95, .95)
        lower, upper = np.full(len(idx), -.999999), np.full(len(idx), .999999)
        active = np.ones(len(idx), dtype=bool)
        counts = np.zeros(len(idx), dtype=int)
        for iteration in range(config.max_iterations):
            which = np.flatnonzero(active)
            if not len(which):
                break
            rr, aa, bb = r[which], a[which], b[which]
            error = bivariate_normal_cdf(aa, bb, rr, config.quadrature_nodes)-target[which]
            counts[which] += 1
            done = np.abs(error) <= config.probability_tolerance
            active[which[done]] = False
            unfinished = ~done
            q, e = which[unfinished], error[unfinished]
            lower[q] = np.where(e < 0, r[q], lower[q])
            upper[q] = np.where(e > 0, r[q], upper[q])
            density = np.exp(-(a[q]**2-2*a[q]*b[q]*r[q]+b[q]**2)/(2*(1-r[q]**2)))/(2*np.pi*np.sqrt(1-r[q]**2))
            with np.errstate(divide="ignore", invalid="ignore"):
                proposal = r[q]-e/density
            safe = np.isfinite(proposal) & (proposal>lower[q]) & (proposal<upper[q])
            r[q] = np.where(safe, proposal, (lower[q]+upper[q])/2)
        values32 = bivariate_normal_cdf(a, b, r, config.quadrature_nodes)
        values64 = bivariate_normal_cdf(a, b, r, max(64, 2*config.quadrature_nodes))
        err = np.abs(values64-target)
        quadrature_error = np.abs(values64-values32)
        good = (err <= config.probability_tolerance*2) & (quadrature_error <= config.quadrature_tolerance)
        methods[idx] = "safeguarded vector Newton; Gauss-Legendre verified at doubled nodes"
        if config.max_iterations > 0:
            for local in np.flatnonzero(~good):
                try:
                    root = brentq(lambda value: adaptive_cdf(a[local], b[local], value)-target[local], -.999999, .999999,
                                  xtol=1e-11, maxiter=150)
                    r[local] = root
                    err[local] = abs(adaptive_cdf(a[local], b[local], root)-target[local])
                    quadrature_error[local] = abs(adaptive_cdf(a[local], b[local], root)-float(bivariate_normal_cdf(a[local], b[local], root, 64)))
                    good[local] = err[local] <= config.probability_tolerance*2
                    methods[idx[local]] = "adaptive quadrature + Brent fallback"
                except (ValueError, RuntimeError):
                    methods[idx[local]] = "adaptive quadrature + Brent failed"
        rho[idx], errors[idx], quadrature_errors[idx] = r, err, quadrature_error
        iterations[idx], converged[idx] = counts, good
        statuses[idx] = np.where(good, "converged", "numerical nonconvergence")
        rho[idx[~good]] = np.nan
    return pd.DataFrame({"n_pair": n.astype(int), "n00": cells[:,0], "n01": cells[:,1], "n10": cells[:,2], "n11": cells[:,3],
                         "sparse_table": (cells < config.sparse_cell_below).any(axis=1), "continuity_corrected": corrected,
                         "correction_added_to_each_cell": corrected*config.continuity_correction,
                         "tetrachoric_r": rho, "converged": converged, "iterations": iterations,
                         "absolute_probability_error": errors, "quadrature_check_difference": quadrature_errors,
                         "boundary_estimate": np.abs(rho)>.98, "estimation_method": methods, "status": statuses})


def tetrachoric_matrix(x, item_ids=None, config=None):
    x = np.asarray(x, float)
    if x.ndim != 2 or np.isinf(x).any() or not np.isin(x[np.isfinite(x)], [0, 1]).all():
        raise ValueError("Expected 2D binary observations, with NaN allowed for missing")
    p = x.shape[1]
    names = list(item_ids) if item_ids is not None else [f"item_{j}" for j in range(p)]
    if len(names) != p or len(set(names)) != p:
        raise ValueError("Invalid or duplicated item ordering")
    valid, zero = np.isfinite(x).astype(float), (x == 0).astype(float)
    n, zero_margin, n00 = valid.T@valid, zero.T@valid, zero.T@zero
    ii, jj = np.triu_indices(p, 1)
    a = n00[ii,jj]
    b, c = zero_margin[ii,jj]-a, zero_margin[jj,ii]-a
    d = n[ii,jj]-a-b-c
    pairs = estimate_tables(np.column_stack([a,b,c,d]), config)
    pairs.insert(0, "item_j", [names[j] for j in jj])
    pairs.insert(0, "item_i", [names[i] for i in ii])
    matrix = np.eye(p)
    matrix[ii,jj], matrix[jj,ii] = pairs.tetrachoric_r, pairs.tetrachoric_r
    return matrix, pairs, n.astype(int)
