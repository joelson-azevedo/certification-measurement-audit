"""Exact signed-model number-correct tail probabilities, reusing Stage 8 recursion."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from .score_distribution import conditional_pmf,response_probabilities

DISCLOSURE="MODEL-BASED DECISION-CONSISTENCY RESEARCH — NOT A VALIDATION OF THE PASSING STANDARD OR AN AUTHORIZATION TO CHANGE CERTIFICATION DECISIONS."


def score_tail(pmf,cut):
    pmf=np.asarray(pmf,float)
    if not isinstance(cut,(int,np.integer)) or not 0<=cut<=pmf.shape[-1]:raise ValueError("An integer cut in 0..K+1 is required")
    if pmf.ndim<1 or not np.isfinite(pmf).all() or (pmf<0).any() or not np.allclose(pmf.sum(axis=-1),1,rtol=0,atol=1e-11):
        raise ValueError("Normalized finite nonnegative conditional PMF required")
    # Do not clip or normalize the recursive PMF; retain roundoff and verify tolerance.
    return pmf[...,cut:].sum(axis=-1)


def passing_curves(theta,d,a,cuts):
    probabilities=response_probabilities(theta,d,a);pmf=conditional_pmf(probabilities)
    tails={int(cut):score_tail(pmf,cut) for cut in cuts}
    return tails,pmf,probabilities


def limiting_pass_probability(d,a,cut,positive_infinity):
    d,a=np.asarray(d,float),np.asarray(a,float)
    p=np.where(a==0,1/(1+np.exp(-d)),np.where(a>0,1. if positive_infinity else 0.,0. if positive_infinity else 1.))
    return float(score_tail(conditional_pmf(p),cut))


def passing_review(theta,tail,d,a,cut,tolerance=1e-12):
    theta=np.asarray(theta);tail=np.asarray(tail)
    differences=np.diff(tail);bad=np.flatnonzero(differences < -tolerance)
    return {"cut":cut,"n_theta":len(theta),"theta_lower":theta[0],"theta_upper":theta[-1],
            "minimum_passing_probability":float(tail.min()),"maximum_passing_probability":float(tail.max()),
            "minimum_passing_increment":float(differences.min()),"n_decreasing_intervals_beyond_tolerance":len(bad),
            "monotonic_on_checked_grid":len(bad)==0,"passing_at_grid_lower":float(tail[0]),"passing_at_grid_upper":float(tail[-1]),
            "analytic_limit_theta_negative_infinity":limiting_pass_probability(d,a,cut,False),
            "analytic_limit_theta_positive_infinity":limiting_pass_probability(d,a,cut,True),
            "interpretation":"aggregate numerical check; individual slopes remain signed; no global monotonicity assumption"}


def probability_band_intervals(theta,tail,lower,upper,tail_function):
    """All connected grid-detected components, with exact-tail root refinement."""
    theta,tail=np.asarray(theta,float),np.asarray(tail,float)
    if not 0<=lower<upper<=1:raise ValueError("Ordered probability band in [0,1] required")
    inside=(tail>=lower)&(tail<=upper);indices=np.flatnonzero(inside);rows=[]
    if not len(indices):return pd.DataFrame(columns=["component","probability_lower","probability_upper","theta_lower","theta_upper","truncated_by_grid"])
    blocks=np.split(indices,np.flatnonzero(np.diff(indices)>1)+1)
    for component,block in enumerate(blocks,1):
        i,j=block[0],block[-1];left,right=theta[i],theta[j]
        if i>0:
            level=lower if tail[i-1]<lower else upper
            left=brentq(lambda value:float(tail_function(value))-level,theta[i-1],theta[i],xtol=1e-12)
        if j<len(theta)-1:
            level=lower if tail[j+1]<lower else upper
            right=brentq(lambda value:float(tail_function(value))-level,theta[j],theta[j+1],xtol=1e-12)
        rows.append({"component":component,"probability_lower":lower,"probability_upper":upper,
                     "theta_lower":left,"theta_upper":right,"truncated_by_grid":bool(i==0 or j==len(theta)-1)})
    return pd.DataFrame(rows)
