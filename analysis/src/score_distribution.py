"""Lord–Wingersky Bernoulli recursion and common-reference marginalization."""
from __future__ import annotations
import numpy as np
from scipy.special import expit,roots_hermitenorm

DISCLOSURE="PROVISIONAL RESEARCH EQUATING — NOT AUTHORIZED FOR OPERATIONAL SCORING OR CERTIFICATION DECISIONS."


def response_probabilities(theta,d,a):
    d,a=np.asarray(d,float),np.asarray(a,float);theta=np.asarray(theta,float)
    if d.ndim!=1 or d.shape!=a.shape or not len(d) or not all(np.isfinite(v).all() for v in [theta,d,a]):
        raise ValueError("Finite matching intercept/slope vectors required")
    return expit(theta[...,None]*a+d)


def conditional_pmf(probabilities):
    """Each leading row is a theta; no slope-sign assumption or normalization fix."""
    p=np.asarray(probabilities,float)
    if p.ndim<1 or p.shape[-1]==0 or not np.isfinite(p).all() or ((p<0)|(p>1)).any():
        raise ValueError("Finite Bernoulli probabilities in [0,1] required")
    k=p.shape[-1];shape=p.shape[:-1];flat=p.reshape(-1,k)
    pmf=np.zeros((len(flat),k+1));pmf[:,0]=1.
    for j in range(k):
        old=pmf[:,:j+1].copy();success=flat[:,j,None]
        pmf[:,:j+2]=0.
        pmf[:,:j+1]+=old*(1-success)
        pmf[:,1:j+2]+=old*success
    return pmf.reshape(shape+(k+1,))


def conditional_diagnostics(p,pmf):
    p=np.asarray(p,float);pmf=np.asarray(pmf,float);scores=np.arange(pmf.shape[-1])
    mean=pmf@scores;variance=pmf@(scores*scores)-mean*mean
    return {"all_probabilities_finite":bool(np.isfinite(pmf).all()),"minimum_probability":float(pmf.min()),
            "maximum_normalization_error":float(abs(pmf.sum(axis=-1)-1).max()),
            "maximum_mean_TCC_error":float(abs(mean-p.sum(axis=-1)).max()),
            "maximum_variance_Bernoulli_error":float(abs(variance-(p*(1-p)).sum(axis=-1)).max())}


def normal_quadrature(nodes,mean=0.,sd=1.):
    if nodes<3 or not np.isfinite(mean) or not np.isfinite(sd) or sd<=0:raise ValueError("Positive normal scale and at least 3 nodes required")
    z,w=roots_hermitenorm(nodes)
    return mean+sd*z,w/np.sqrt(2*np.pi)


def reference_quadrature(name,nodes,scale=1.,shift=0.):
    if scale<=0 or not np.isfinite(scale) or not np.isfinite(shift):raise ValueError("Invalid linking constants")
    if name=="A_standard_normal":return normal_quadrature(nodes)
    if name=="B_identification_on_A":return normal_quadrature(nodes,shift,scale)
    if name=="equal_identification_mixture":
        a,wa=normal_quadrature(nodes);b,wb=normal_quadrature(nodes,shift,scale)
        order=np.argsort(np.r_[a,b],kind="stable")
        return np.r_[a,b][order],np.r_[wa,wb][order]/2
    raise ValueError("Unknown common reference distribution")


def marginal_distribution(d,a,theta,weights):
    theta,weights=np.asarray(theta,float),np.asarray(weights,float)
    if theta.ndim!=1 or weights.shape!=theta.shape or not np.isfinite(weights).all() or (weights<0).any() or not np.isclose(weights.sum(),1,atol=1e-12,rtol=0):
        raise ValueError("Same finite normalized common-reference nodes/weights required")
    p=response_probabilities(theta,d,a);conditional=conditional_pmf(p);pmf=weights@conditional
    scores=np.arange(len(pmf));mean=float(pmf@scores);variance=float(pmf@(scores*scores)-mean*mean)
    conditional_mean=p.sum(axis=1);conditional_var=(p*(1-p)).sum(axis=1)
    diagnostic={**conditional_diagnostics(p,conditional),"marginal_normalization_error":float(abs(pmf.sum()-1)),
                "marginal_minimum_probability":float(pmf.min()),"marginal_mean":mean,"marginal_variance":variance,
                "marginal_mean_integral_error":float(abs(mean-weights@conditional_mean)),
                "marginal_variance_total_variance_error":float(abs(variance-(weights@(conditional_var+conditional_mean**2)-(weights@conditional_mean)**2))),
                "cdf_minimum_increment":float(np.diff(np.cumsum(pmf)).min()),"n_integration_nodes":len(theta)}
    return {"pmf":pmf,"cdf":np.cumsum(pmf),"conditional":conditional,"probabilities":p,"theta":theta,"weights":weights,"diagnostics":diagnostic}
