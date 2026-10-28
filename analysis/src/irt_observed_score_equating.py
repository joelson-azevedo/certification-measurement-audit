"""Explicit integer mid-CDF convention; unsupported tails never clipped."""
from __future__ import annotations
import numpy as np


def mid_cdf(pmf):
    pmf=np.asarray(pmf,float)
    if pmf.ndim!=1 or len(pmf)<2 or not np.isfinite(pmf).all() or (pmf<0).any() or not np.isclose(pmf.sum(),1,atol=1e-11,rtol=0):
        raise ValueError("Normalized nonnegative score PMF required")
    # Sum lower and upper tails separately to avoid loss of precision in extreme tails.
    lower=np.cumsum(pmf)-pmf/2
    upper=np.cumsum(pmf[::-1])[::-1]-pmf/2
    return np.where(lower<=.5,lower,1-upper)


def inverse_mid_cdf(percentiles,pmf):
    p=np.asarray(percentiles,float);knots=mid_cdf(pmf);scores=np.arange(len(knots),dtype=float)
    values=np.full(p.shape,np.nan);status=np.full(p.shape,"unavailable: percentile outside reference integer-score mid-CDF range",dtype=object)
    for index in np.ndindex(p.shape):
        value=p[index]
        if not np.isfinite(value) or value<knots[0] or value>knots[-1]:continue
        equal=np.flatnonzero(knots==value)
        if len(equal)>1:
            status[index]="unavailable: nonunique flat percentile interval";continue
        if len(equal)==1:values[index]=scores[equal[0]];status[index]="available";continue
        right=np.searchsorted(knots,value,side="right");left=right-1
        if right>=len(knots) or left<0 or knots[right]<=knots[left]:
            status[index]="unavailable: nonunique flat percentile interval";continue
        values[index]=left+(value-knots[left])/(knots[right]-knots[left]);status[index]="available"
    return values,status


def mid_rank_at(scores,pmf):
    scores=np.asarray(scores,float);grid=np.arange(len(pmf));r=mid_cdf(pmf)
    valid=np.isfinite(scores)&(scores>=0)&(scores<=len(pmf)-1)
    result=np.full(scores.shape,np.nan);result[valid]=np.interp(scores[valid],grid,r)
    return result


def equate_observed(source_scores,source_pmf,reference_pmf):
    percentiles=mid_rank_at(source_scores,source_pmf)
    values,status=inverse_mid_cdf(percentiles,reference_pmf)
    return {"equated":values,"source_percentile":percentiles,"reference_percentile":mid_rank_at(values,reference_pmf),"status":status}


def integer_interpretations(values):
    x=np.asarray(values,float)
    return {"integer_floor":np.floor(x),"integer_nearest_half_up":np.floor(x+.5),"integer_ceiling":np.ceil(x)}


def support_status(source_percentile,source_location,reference_location,overlap,available=True,central=(.025,.975)):
    if not available:return "unavailable/unstable"
    if (central[0]<=source_percentile<=central[1] and np.isfinite(source_location) and np.isfinite(reference_location)
        and overlap[0]<=source_location<=overlap[1] and overlap[0]<=reference_location<=overlap[1]):
        return "supported central research region"
    return "weakly supported tail/model extrapolation"
