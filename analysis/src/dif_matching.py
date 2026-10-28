"""Fixed-calibration, focal-excluded DIF matching and bounded matching pools."""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .irt_2pl import conditional_log_likelihood, normal_quadrature, validate_responses


def encode_regions(values):
    return np.array([0. if value=="South" else 1. if value=="North" else np.nan for value in values])


def matching_pool(members,removed=()):
    forbidden=set(removed)|{"ANC_14"}
    return [item for item in members if item not in forbidden]


def rest_eap_matrix(x,intercept,slope,members,pool=None,nodes=81):
    x=validate_responses(x);d=np.asarray(intercept,float);a=np.asarray(slope,float)
    if x.shape[1]!=len(members) or len(set(members))!=len(members) or d.shape!=a.shape or d.shape!=(len(members),):
        raise ValueError("Matching parameter/member dimension mismatch")
    pool=matching_pool(members) if pool is None else list(pool)
    if not pool or not set(pool)<=set(members) or len(pool)!=len(set(pool)) or "ANC_14" in pool:
        raise ValueError("Invalid common matching pool; ANC_14 cannot anchor matching")
    keep=[members.index(item) for item in pool]
    theta,weights=normal_quadrature(nodes)
    base=conditional_log_likelihood(x[:,keep],theta,d[keep],a[keep])+np.log(weights)
    total_observed=np.isfinite(x[:,keep]).sum(axis=1)
    means=np.empty_like(x);sd=np.empty_like(x);counts=np.empty(x.shape,int);planned=[]
    for j,item in enumerate(members):
        if item in pool:
            # Subtract exactly the target's observed likelihood contribution.
            focal=conditional_log_likelihood(x[:,[j]],theta,d[[j]],a[[j]])
            joint=base-focal
            counts[:,j]=total_observed-np.isfinite(x[:,j])
            planned.append(len(pool)-1)
        else:
            joint=base
            counts[:,j]=total_observed
            planned.append(len(pool))
        posterior=np.exp(joint-logsumexp(joint,axis=1)[:,None])
        means[:,j]=posterior@theta
        sd[:,j]=np.sqrt(np.maximum(posterior@(theta*theta)-means[:,j]**2,0))
    return means,sd,counts,np.array(planned)


def rest_score_matrix(x,members,pool=None):
    x=validate_responses(x)
    pool=matching_pool(members) if pool is None else list(pool)
    if not pool or not set(pool)<=set(members) or "ANC_14" in pool:raise ValueError("Unsupported rest-score pool")
    keep=[members.index(item) for item in pool]
    totals=np.nansum(x[:,keep],axis=1);n=np.isfinite(x[:,keep]).sum(axis=1)
    result=np.full_like(x,np.nan);counts=np.empty(x.shape,int);planned=[]
    for j,item in enumerate(members):
        included=item in pool
        observed=n-np.isfinite(x[:,j]) if included else n
        correct=totals-np.nan_to_num(x[:,j],nan=0) if included else totals
        size=len(pool)-int(included)
        result[:,j]=np.divide(size*correct,observed,out=np.full(len(x),np.nan),where=observed>0)
        counts[:,j]=observed;planned.append(size)
    return result,counts,np.array(planned)


def center_matching(values,observed_count,group):
    values=np.asarray(values,float);observed_count=np.asarray(observed_count)
    group=np.asarray(group,float)
    valid=(observed_count>0)&np.isfinite(values)&np.isfinite(group)[:,None]
    centers=np.array([values[valid[:,j],j].mean() if valid[:,j].any() else np.nan for j in range(values.shape[1])])
    return values-centers,centers,valid


def purification_pool(members,metadata,removed,screened,minimum_matching_items=30,minimum_per_domain=3):
    removed=set(removed)|set(screened)|{"ANC_14"}
    pool=matching_pool(members,removed)
    meta=metadata.set_index("item_id")
    domain_counts=meta.loc[pool,"domain"].value_counts()
    all_domains=set(meta.loc[members,"domain"])
    too_small=len(pool)<minimum_matching_items
    inadequate=any(domain_counts.get(domain,0)<minimum_per_domain for domain in all_domains)
    if too_small or inadequate:
        return {"pool":pool,"removed":sorted(removed&set(members)),"supported":False,
                "status":"stopped: insufficient matching length/domain coverage","domain_counts":domain_counts.to_dict()}
    return {"pool":pool,"removed":sorted(removed&set(members)),"supported":True,
            "status":"supported bounded sensitivity pool","domain_counts":domain_counts.to_dict()}
