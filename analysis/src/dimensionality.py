"""Missing-aware matrix preparation, matrix diagnostics and tetrachoric MINRES."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .ctt import form_matrices, select_population


def response_specifications(admin, long, metadata, form, excluded_item="ANC_14"):
    selected, responses = select_population(admin, long, True)
    sittings, members, operational, marked = form_matrices(selected, responses, metadata, form)
    if len(members) != 60:
        raise ValueError("Stage 3 requires the complete 60-item form")
    latent = operational.copy()
    latent[marked == ""] = np.nan
    complete = np.isfinite(latent).all(axis=1)
    if excluded_item not in members:
        raise ValueError("Prespecified Stage 2 sensitivity item is absent")
    keep = np.array([item != excluded_item for item in members])
    specs = {
        "primary_complete_cases": (latent[complete], members, sittings.loc[complete].copy()),
        "exclude_ANC_14_matched_cases": (latent[complete][:,keep], list(np.asarray(members)[keep]), sittings.loc[complete].copy()),
        "omissions_zero_all_first": (operational, members, sittings.copy()),
    }
    info = {"form": form, "n_first_attempts": len(sittings), "n_complete_cases": int(complete.sum()),
            "n_excluded_for_omission": int((~complete).sum()), "n_omitted_cells": int(np.isnan(latent).sum()),
            "n_items_full_form": len(members), "retention_proportion": float(complete.mean()),
            "primary_missingness_method": "all-item complete cases; administered omissions missing",
            "sensitivity_exclusion_cases": "same complete candidates as full primary form; no case recovery"}
    return specs, info, latent, sittings


def matrix_diagnostics(matrix, tolerance=1e-8):
    r = np.asarray(matrix, float)
    square = r.ndim == 2 and r.shape[0] == r.shape[1]
    finite = bool(np.isfinite(r).all())
    result = {"square": square, "finite": finite, "n_items": r.shape[0] if square else 0}
    if not square or not finite:
        return {**result, "valid_correlation_structure": False, "positive_definite": False,
                "minimum_eigenvalue": np.nan, "maximum_eigenvalue": np.nan, "condition_number": np.nan,
                "symmetry_max_error": np.nan, "diagonal_max_error": np.nan, "negative_eigenvalues": -1}
    sym = np.max(np.abs(r-r.T))
    diag = np.max(np.abs(np.diag(r)-1))
    eigen = np.linalg.eigvalsh((r+r.T)/2)
    valid = sym <= tolerance and diag <= tolerance and np.max(np.abs(r)) <= 1+tolerance
    return {**result, "valid_correlation_structure": bool(valid), "positive_definite": bool(valid and eigen[0] > tolerance),
            "minimum_eigenvalue": float(eigen[0]), "maximum_eigenvalue": float(eigen[-1]),
            "condition_number": float(np.linalg.cond(r)), "symmetry_max_error": float(sym),
            "diagonal_max_error": float(diag), "negative_eigenvalues": int((eigen < -tolerance).sum())}


def ordered_eigenvalues(r):
    d = matrix_diagnostics(r)
    if not d["valid_correlation_structure"]:
        raise ValueError("Invalid correlation matrix")
    return np.linalg.eigvalsh(r)[::-1]


def condition_for_factor_analysis(r, floor=1e-6, maximum_adjustment=.05, maximum_relative_frobenius=.05):
    """Explicit Higham alternating projections with Dykstra correction, if required.

    Raw matrices are never overwritten. PCA/parallel analysis use raw tetrachorics;
    only EFA uses the conditioned copy, and only within prespecified small-change
    safeguards. These are numerical safeguards, not psychometric fit standards.
    """
    r = np.asarray(r, float)
    diagnostic = matrix_diagnostics(r)
    if not diagnostic["valid_correlation_structure"]:
        raise ValueError("Cannot condition an invalid/nonfinite correlation matrix")
    if diagnostic["minimum_eigenvalue"] >= floor:
        return r.copy(), {"correction_applied": False, "correction_method": "none", "correction_iterations": 0,
                          "correction_converged": True, "maximum_absolute_adjustment": 0.0, "frobenius_adjustment": 0.0,
                          "relative_frobenius_adjustment": 0.0, "pairs_changed_over_0_01": 0, "efa_numerically_supported": True}
    y, delta = r.copy(), np.zeros_like(r)
    converged = False
    for iteration in range(1, 501):
        residual = y-delta
        values, vectors = np.linalg.eigh(residual)
        x = (vectors*np.maximum(values, 0))@vectors.T
        delta = x-residual
        new = x.copy()
        np.fill_diagonal(new, 1)
        change = np.linalg.norm(new-y, "fro")
        y = new
        if change < 1e-10:
            converged = True
            break
    values, vectors = np.linalg.eigh((y+y.T)/2)
    positive = (vectors*np.maximum(values, floor))@vectors.T
    scale = np.sqrt(np.diag(positive))
    adjusted = positive/scale[:,None]/scale[None,:]
    difference = adjusted-r
    maximum, frobenius = np.max(np.abs(difference)), np.linalg.norm(difference, "fro")
    relative = frobenius/np.linalg.norm(r, "fro")
    return adjusted, {"correction_applied": True, "correction_method": "Higham PSD/unit-diagonal projections with Dykstra; eigenvalue floor 1e-6 and diagonal rescaling",
                      "correction_iterations": iteration, "correction_converged": converged,
                      "maximum_absolute_adjustment": float(maximum), "frobenius_adjustment": float(frobenius),
                      "relative_frobenius_adjustment": float(relative),
                      "pairs_changed_over_0_01": int(np.count_nonzero(np.triu(np.abs(difference)>.01, 1))),
                      "efa_numerically_supported": bool(converged and maximum<=maximum_adjustment and relative<=maximum_relative_frobenius)}


def varimax(loadings, max_iterations=1000, tolerance=1e-9):
    p, k = loadings.shape
    scale = np.sqrt((loadings**2).sum(axis=1))
    safe = np.where(scale>1e-12, scale, 1)
    normalized = loadings/safe[:,None]
    rotation = np.eye(k)
    previous = 0.0
    for iteration in range(1, max_iterations+1):
        rotated = normalized@rotation
        u, values, vt = np.linalg.svd(normalized.T@(rotated**3-rotated@np.diag((rotated**2).sum(axis=0))/p))
        rotation = u@vt
        objective = values.sum()
        if abs(objective-previous) <= tolerance*max(1, abs(previous)):
            return loadings@rotation, rotation, True, iteration
        previous = objective
    return loadings@rotation, rotation, False, max_iterations


def promax(loadings, power=4):
    if loadings.shape[1] == 1:
        sign = 1 if loadings[np.argmax(np.abs(loadings[:,0])),0] >= 0 else -1
        return loadings*sign, np.ones((1,1)), True, "none; one factor"
    orthogonal, rotation, converged, _ = varimax(loadings)
    target = np.sign(orthogonal)*np.abs(orthogonal)**power
    regression = np.linalg.lstsq(orthogonal, target, rcond=None)[0]
    transform = rotation@regression
    if np.linalg.cond(transform) > 1e10:
        raise ValueError("Promax transform is ill-conditioned")
    covariance = np.linalg.inv(transform.T@transform)
    normalization = np.diag(np.sqrt(np.diag(covariance)))
    transform = transform@normalization
    pattern = loadings@transform
    phi = np.linalg.inv(transform.T@transform)
    signs = np.array([1 if pattern[np.argmax(np.abs(pattern[:,j])),j]>=0 else -1 for j in range(pattern.shape[1])])
    pattern = pattern*signs
    phi = phi*signs[:,None]*signs[None,:]
    order = np.argsort(-np.sum((pattern@phi)**2,axis=0), kind="stable")
    return pattern[:,order], phi[np.ix_(order,order)], converged, "oblique Promax power=4 after Kaiser-normalized Varimax"


def minres(r, n_factors, max_iterations=1000):
    """ULS/MINRES uniqueness optimization with analytic spectral gradient.

    Minimize squared residual eigenvalues of R-diag(psi) after its best positive
    rank-m approximation. At interior optima uniqueness absorbs the diagonal and
    the objective minimizes off-diagonal common-factor residuals. This is common
    factor extraction, not principal-components loadings from the full R.
    """
    r = np.asarray(r, float)
    diagnostic = matrix_diagnostics(r)
    p = len(r)
    if not diagnostic["valid_correlation_structure"] or not diagnostic["positive_definite"]:
        raise ValueError("EFA requires a documented positive-definite analysis matrix")
    if n_factors < 1 or n_factors >= p or (p-n_factors)**2 < p+n_factors:
        raise ValueError("Requested factor solution is not identified by covariance degrees of freedom")

    def objective(uniqueness):
        reduced = r-np.diag(uniqueness)
        values, vectors = np.linalg.eigh(reduced)
        load = vectors[:,-n_factors:]*np.sqrt(np.maximum(values[-n_factors:], 0))
        residual = reduced-load@load.T
        return float(np.sum(residual**2)), -2*np.diag(residual)

    attempts = []
    for initial in [.5, .8]:
        result = minimize(objective, np.full(p, initial), method="L-BFGS-B", jac=True, bounds=[(1e-5,1)]*p,
                          options={"maxiter": max_iterations, "ftol": 1e-12, "gtol": 1e-7, "maxls": 40})
        attempts.append(result)
    eligible = [fit for fit in attempts if fit.success] or attempts
    result = min(eligible, key=lambda fit: fit.fun)
    values, vectors = np.linalg.eigh(r-np.diag(result.x))
    unrotated = vectors[:,-n_factors:][:,::-1]*np.sqrt(np.maximum(values[-n_factors:][::-1],0))
    try:
        pattern, phi, rotation_converged, rotation_name = promax(unrotated)
        rotation_warning = ""
    except (ValueError, np.linalg.LinAlgError) as exc:
        pattern, phi, rotation_converged, rotation_name = unrotated, np.eye(n_factors), False, "rotation failed; unrotated retained as diagnostic"
        rotation_warning = str(exc)
    implied = pattern@phi@pattern.T
    communalities = np.diag(implied)
    residual = r-implied
    off_diagonal = residual[np.triu_indices(p,1)]
    heywood = (communalities>1+1e-6) | (communalities< -1e-8)
    boundary = (result.x <= 1.0001e-5)
    warnings = []
    if not result.success:
        warnings.append("optimizer did not converge")
    if not rotation_converged:
        warnings.append("rotation did not converge: "+rotation_warning)
    if boundary.any():
        warnings.append("uniqueness lower-bound hits: potential Heywood constraint")
    if heywood.any():
        warnings.append("inadmissible communality")
    diagnostics = {"n_factors": n_factors, "converged": bool(result.success), "optimizer_message": str(result.message),
                   "optimizer_iterations": int(result.nit), "objective": float(result.fun), "optimizer_starts": len(attempts),
                   "rotation": rotation_name, "rotation_converged": rotation_converged,
                   "rmsr_off_diagonal": float(np.sqrt(np.mean(off_diagonal**2))),
                   "mean_absolute_residual": float(np.mean(np.abs(off_diagonal))),
                   "residual_absolute_95th_percentile": float(np.percentile(np.abs(off_diagonal),95)),
                   "maximum_absolute_residual": float(np.max(np.abs(off_diagonal))),
                   "communalities_min": float(communalities.min()), "communalities_max": float(communalities.max()),
                   "n_heywood_cases": int(heywood.sum()), "n_uniqueness_lower_bound_hits": int(boundary.sum()),
                   "n_cross_loading_items_abs_0_30": int(((np.abs(pattern)>=.30).sum(axis=1)>1).sum()),
                   "factor_correlation_max_abs_off_diagonal": float(np.max(np.abs(phi-np.eye(n_factors)))) if n_factors>1 else 0,
                   "admissible": bool(not heywood.any() and np.linalg.eigvalsh(phi).min()>0), "warnings": "; ".join(warnings)}
    return {"pattern": pattern, "structure": pattern@phi, "phi": phi, "communalities": communalities,
            "uniqueness": 1-communalities, "optimized_uniqueness": result.x,
            "unrotated": unrotated, "diagnostics": diagnostics}


def factor_tables(result, metadata, members):
    meta = metadata.set_index("item_id").loc[members]
    loadings, correlations, domains = [], [], []
    pattern, structure, phi = result["pattern"], result["structure"], result["phi"]
    factors = pattern.shape[1]
    for j,item in enumerate(members):
        for factor in range(factors):
            loadings.append({"n_factors": factors, "item_id": item, "domain": meta.at[item,"domain"],
                             "anchor_item": meta.at[item,"anchor_item"], "factor": factor+1,
                             "pattern_loading": pattern[j,factor], "structure_loading": structure[j,factor],
                             "communality": result["communalities"][j], "uniqueness": result["uniqueness"][j],
                             "optimized_uniqueness": result["optimized_uniqueness"][j],
                             "n_pattern_loadings_abs_ge_0_30": int((np.abs(pattern[j])>=.30).sum()),
                             "salient_pattern_loading": bool(abs(pattern[j,factor])>=.30)})
    for i in range(factors):
        for j in range(factors):
            correlations.append({"n_factors": factors, "factor_i": i+1, "factor_j": j+1, "factor_correlation": phi[i,j]})
    for factor in range(factors):
        high = np.abs(pattern[:,factor]) >= .30
        all_high = int(high.sum())
        for domain in metadata.domain.unique():
            belongs = meta.domain.eq(domain).to_numpy()
            count = int((belongs & high).sum())
            domains.append({"n_factors": factors, "factor": factor+1, "domain": domain,
                            "n_domain_items": int(belongs.sum()), "n_abs_pattern_ge_0_30": count,
                            "proportion_of_factor_salient_items": count/all_high if all_high else np.nan,
                            "proportion_of_domain_items_salient": count/belongs.sum() if belongs.sum() else np.nan,
                            "mean_absolute_pattern_loading": float(np.abs(pattern[belongs,factor]).mean()) if belongs.any() else np.nan,
                            "mean_absolute_structure_loading": float(np.abs(structure[belongs,factor]).mean()) if belongs.any() else np.nan})
    return pd.DataFrame(loadings), pd.DataFrame(correlations), pd.DataFrame(domains)
