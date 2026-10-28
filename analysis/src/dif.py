"""Nested logistic DIF, separation safeguards, BH and supported probability effects."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.special import expit, xlogy
from scipy.stats import chi2

Z95=1.959963984540054


@dataclass(frozen=True)
class LogisticConfig:
    max_iterations: int=100
    gradient_tolerance_per_candidate: float=1e-8
    information_maximum_condition_number: float=1e10


def logistic_loglikelihood(beta,x,y):
    eta=np.asarray(x,float)@np.asarray(beta,float)
    y=np.asarray(y,float)
    return float(np.sum(y*eta-np.logaddexp(0,eta)))


def separation_diagnostic(x,y):
    """Find a nonzero separating signed margin by bounded-direction linear programs.

    Column scaling only aids numerics. If all signed margins are >=0 and their
    average can exceed 1e-7, a complete/quasi separating direction exists.
    A second max-min-margin LP distinguishes strict complete separation.
    No logistic coefficient bound, clipping or penalized likelihood is used.
    """
    x=np.asarray(x,float);y=np.asarray(y,float)
    scale=np.sqrt(np.mean(x*x,axis=0));scale[scale==0]=1
    signed=(2*y-1)[:,None]*(x/scale)
    result=linprog(-signed.mean(axis=0),A_ub=-signed,b_ub=np.zeros(len(y)),
                   bounds=[(-1,1)]*x.shape[1],method="highs")
    if not result.success:
        return {"separation_status":"unavailable: separation LP failed","separation_margin_mean":np.nan}
    mean_margin=float(-result.fun)
    if mean_margin<=1e-7:
        return {"separation_status":"none detected","separation_margin_mean":mean_margin}
    constraints=np.column_stack([-signed,np.ones(len(y))])
    strict=linprog(np.r_[np.zeros(x.shape[1]),-1.],A_ub=constraints,b_ub=np.zeros(len(y)),
                   bounds=[(-1,1)]*x.shape[1]+[(0,None)],method="highs")
    status=("separation detected; strict-margin LP unavailable" if not strict.success else
            "complete separation" if strict.x[-1]>1e-7 else "quasi-complete separation")
    return {"separation_status":status,"separation_margin_mean":mean_margin}


def fit_logistic(x,y,config=LogisticConfig(),initial=None):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if x.ndim!=2 or y.shape!=(len(x),) or not np.isfinite(x).all() or not np.isin(y,[0,1]).all():
        raise ValueError("Finite logistic design and binary target required")
    p=x.shape[1];n=len(y)
    unavailable={"beta":np.full(p,np.nan),"covariance":np.full((p,p),np.nan),"log_likelihood":np.nan,
                 "converged":False,"numerically_supported":False,"iterations":0,
                 "maximum_scaled_gradient_per_candidate":np.nan,"information_condition_number":np.nan,
                 "information_minimum_eigenvalue":np.nan,"n":n,"n_correct":int(y.sum())}
    if n<=p or np.linalg.matrix_rank(x)<p:
        return {**unavailable,"status":"unavailable: insufficient observations or rank-deficient design","separation_status":"not evaluated"}
    separation=separation_diagnostic(x,y)
    if separation["separation_status"]!="none detected":
        return {**unavailable,**separation,"status":"unavailable: "+separation["separation_status"]}
    scale=np.sqrt(np.mean(x*x,axis=0));scale[scale==0]=1
    z=x/scale
    beta=np.zeros(p) if initial is None else np.asarray(initial,float).copy()
    if initial is None and np.allclose(x[:,0],1):
        mean=y.mean();beta[0]=np.log(mean/(1-mean))
    gamma=beta*scale
    converged=False;status="maximum iterations reached";iterations=0
    for iteration in range(config.max_iterations):
        iterations=iteration+1
        eta=z@gamma;prob=expit(eta);gradient=z.T@(prob-y)
        if np.max(abs(gradient))/n<=config.gradient_tolerance_per_candidate:
            converged=True;status="converged: first-order condition";break
        h=z.T@((prob*(1-prob))[:,None]*z)
        try:step=np.linalg.solve(h,gradient)
        except np.linalg.LinAlgError:
            status="unavailable: singular Newton information";break
        objective=float(np.sum(np.logaddexp(0,eta)-y*eta))
        factor=1.;accepted=False
        for _ in range(40):
            proposal=gamma-factor*step
            proposal_objective=-logistic_loglikelihood(proposal,z,y)
            if proposal_objective<=objective-1e-4*factor*(gradient@step)+1e-12:
                gamma=proposal;accepted=True;break
            factor*=.5
        if not accepted:
            status="unavailable: Newton line search failed";break
    beta=gamma/scale
    prob=expit(x@beta);gradient=z.T@(prob-y)
    h=x.T@((prob*(1-prob))[:,None]*x)
    eigen=np.linalg.eigvalsh(h);condition=float(np.linalg.cond(h))
    supported=(converged and eigen.min()>1e-10 and condition<config.information_maximum_condition_number)
    covariance=np.linalg.solve(h,np.eye(p)) if supported else np.full((p,p),np.nan)
    return {"beta":beta if supported else np.full(p,np.nan),"covariance":covariance,
            "log_likelihood":logistic_loglikelihood(beta,x,y) if supported else np.nan,
            "converged":converged,"numerically_supported":supported,"iterations":iterations,
            "maximum_scaled_gradient_per_candidate":float(np.max(abs(gradient))/n),
            "information_condition_number":condition,"information_minimum_eigenvalue":float(eigen.min()),
            "n":n,"n_correct":int(y.sum()),"status":status if supported else "unavailable: "+status,
            **separation}


def designs(t,g):
    t=np.asarray(t,float);g=np.asarray(g,float)
    return [np.column_stack([np.ones(len(t)),t]),np.column_stack([np.ones(len(t)),t,g]),
            np.column_stack([np.ones(len(t)),t,g,t*g])]


def fit_nested(t,g,y,config=LogisticConfig()):
    if not np.isin(g,[0,1]).all():raise ValueError("Group must be South=0 / North=1")
    x=designs(t,g);fits=[]
    for matrix in x:
        initial=np.r_[fits[-1]["beta"],0.] if fits and fits[-1]["numerically_supported"] else None
        fits.append(fit_logistic(matrix,y,config,initial))
    tests={"n":len(y),"all_nested_models_supported":all(f["numerically_supported"] for f in fits)}
    for label,i,j,df in [("uniform",0,1,1),("nonuniform",1,2,1),("omnibus",0,2,2)]:
        valid=fits[i]["numerically_supported"] and fits[j]["numerically_supported"]
        statistic=2*(fits[j]["log_likelihood"]-fits[i]["log_likelihood"]) if valid else np.nan
        if valid and statistic < -1e-6:raise FloatingPointError("Nested likelihood order violated")
        statistic=max(0.,statistic) if valid else np.nan
        tests[f"lr_{label}"]=statistic;tests[f"df_{label}"]=df
        tests[f"p_{label}"]=float(chi2.sf(statistic,df)) if valid else np.nan
    return fits,tests


def bh_adjust(pvalues,family_size=None):
    p=np.asarray(pvalues,float)
    if np.any((p[np.isfinite(p)]<0)|(p[np.isfinite(p)]>1)):raise ValueError("Invalid p-value")
    m=len(p) if family_size is None else family_size
    if m<len(p):raise ValueError("Prespecified family cannot be smaller than supplied tests")
    # Unavailable hypotheses retain their family slots conservatively as p=1.
    filled=np.r_[np.where(np.isfinite(p),p,1),np.ones(m-len(p))]
    order=np.argsort(filled,kind="stable")
    adjusted=np.minimum.accumulate((filled[order]*m/np.arange(1,m+1))[::-1])[::-1]
    result=np.empty(m);result[order]=np.minimum(adjusted,1)
    result=result[:len(p)];result[~np.isfinite(p)]=np.nan
    return result


def common_support(t,g,trim_quantile=.05,bins=10,minimum_per_group_per_bin=20,grid_points=201):
    t=np.asarray(t,float);g=np.asarray(g,int)
    empty={"available":False,"status":"unavailable: insufficient shared support","reference_values":np.array([]),
           "reference_weights":np.array([]),"grid":pd.DataFrame(),"bins":pd.DataFrame(),"quantile_values":np.full(3,np.nan)}
    if any((g==group).sum()<minimum_per_group_per_bin for group in [0,1]):return empty
    lower=max(np.quantile(t[g==0],trim_quantile),np.quantile(t[g==1],trim_quantile))
    upper=min(np.quantile(t[g==0],1-trim_quantile),np.quantile(t[g==1],1-trim_quantile))
    if not upper>lower:return empty
    edges=np.linspace(lower,upper,bins+1)
    labels=np.clip(np.searchsorted(edges,t,side="right")-1,0,bins-1)
    inside=(t>=lower)&(t<=upper)
    counts=np.array([[np.sum(inside&(labels==b)&(g==group)) for group in [0,1]] for b in range(bins)])
    supported=np.all(counts>=minimum_per_group_per_bin,axis=1)
    selected=inside&supported[labels]
    if not selected.any() or not supported.any():return empty
    values=t[selected];groups=g[selected]
    weights=np.where(groups==0,.5/(groups==0).sum(),.5/(groups==1).sum())
    order=np.argsort(values,kind="stable");cum=np.cumsum(weights[order]);cum[-1]=1
    quantiles=values[order[np.searchsorted(cum,[.1,.5,.9],side="left")]]
    grid=np.linspace(lower,upper,grid_points)
    grid_labels=np.clip(np.searchsorted(edges,grid,side="right")-1,0,bins-1)
    return {"available":True,"status":"available","lower":float(lower),"upper":float(upper),
            "reference_values":values,"reference_weights":weights,"quantile_values":quantiles,
            "reference_n_south":int((groups==0).sum()),"reference_n_north":int((groups==1).sum()),
            "fraction_supported_south":float(selected[g==0].mean()),"fraction_supported_north":float(selected[g==1].mean()),
            "grid":pd.DataFrame({"matching_centered":grid,"support_bin":grid_labels+1,"supported":supported[grid_labels]}),
            "bins":pd.DataFrame({"support_bin":np.arange(1,bins+1),"lower":edges[:-1],"upper":edges[1:],
                                 "n_south":counts[:,0],"n_north":counts[:,1],"supported":supported})}


def probability_contrasts(beta,covariance,t):
    t=np.asarray(t,float);beta=np.asarray(beta,float);covariance=np.asarray(covariance,float)
    south=np.column_stack([np.ones(len(t)),t,np.zeros(len(t)),np.zeros(len(t))])
    north=np.column_stack([np.ones(len(t)),t,np.ones(len(t)),t])
    ps,pn=expit(south@beta),expit(north@beta)
    grad=pn[:,None]*(1-pn[:,None])*north-ps[:,None]*(1-ps[:,None])*south
    se=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',grad,covariance,grad),0))
    result={"p_south":ps,"p_north":pn,"difference_north_minus_south":pn-ps,
            "difference_se_delta":se,"difference_95_low_pointwise":np.maximum(-1,pn-ps-Z95*se),
            "difference_95_high_pointwise":np.minimum(1,pn-ps+Z95*se)}
    for label,design in [("south",south),("north",north)]:
        eta=design@beta
        eta_se=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',design,covariance,design),0))
        result[f"p_{label}_95_low_pointwise"]=expit(eta-Z95*eta_se)
        result[f"p_{label}_95_high_pointwise"]=expit(eta+Z95*eta_se)
    return pd.DataFrame(result)


def effect_sizes(fits,y,support):
    result={"effect_status":"unavailable: nested model or common support unsupported"}
    if not all(f["numerically_supported"] for f in fits):return result,pd.DataFrame()
    m0,m1,m2=fits;n=len(y);mean=np.mean(y)
    ll_null=float(n*(xlogy(mean,mean)+xlogy(1-mean,1-mean)))
    denominator=-np.expm1(2*ll_null/n)
    r2=[-np.expm1(2*(ll_null-f["log_likelihood"])/n)/denominator for f in fits]
    result.update({"nagelkerke_r2_M0":r2[0],"nagelkerke_r2_M1":r2[1],"nagelkerke_r2_M2":r2[2],
                   "delta_nagelkerke_omnibus":r2[2]-r2[0],"delta_nagelkerke_uniform":r2[1]-r2[0],
                   "delta_nagelkerke_nonuniform":r2[2]-r2[1]})
    for label,model,index in [("M1_group",m1,2),("M2_group",m2,2),("M2_interaction",m2,3)]:
        value=model["beta"][index];se=np.sqrt(model["covariance"][index,index])
        result.update({label+"_coefficient":value,label+"_se":se,label+"_95_low":value-Z95*se,label+"_95_high":value+Z95*se})
    result.update({"M1_uniform_odds_ratio":float(np.exp(result["M1_group_coefficient"])),
                   "M1_uniform_odds_ratio_95_low":float(np.exp(result["M1_group_95_low"])),
                   "M1_uniform_odds_ratio_95_high":float(np.exp(result["M1_group_95_high"])),
                   "odds_ratio_scope":"M1 additive-model North/South odds; not a constant summary of nonuniform M2"})
    if not support["available"]:return result,pd.DataFrame()
    reference=probability_contrasts(m2["beta"],m2["covariance"],support["reference_values"])
    delta=reference.difference_north_minus_south.to_numpy();w=support["reference_weights"]
    result.update({"effect_status":"available; conditional on fixed matching and common support",
                   "signed_reference_mean_difference":float(w@delta),"absolute_reference_mean_difference":float(w@abs(delta))})
    grid=support["grid"].copy()
    prediction=probability_contrasts(m2["beta"],m2["covariance"],grid.matching_centered.to_numpy())
    grid=pd.concat([grid,prediction],axis=1)
    # Unsupported rows carry no released contrast; retain coordinates as gap markers.
    grid.loc[~grid.supported,prediction.columns]=np.nan
    result["maximum_absolute_difference_on_supported_grid"]=float(grid.difference_north_minus_south.abs().max())
    crossing=False;root=np.nan
    if m2["beta"][3]!=0:
        root=-m2["beta"][2]/m2["beta"][3]
        coordinates=grid.matching_centered.to_numpy();valid=grid.supported.to_numpy(bool)
        crossing=bool(np.any(valid[:-1]&valid[1:]&(coordinates[:-1]<root)&(root<coordinates[1:])))
        for index in np.flatnonzero(np.isclose(coordinates,root,atol=1e-12,rtol=0)):
            if 0<index<len(valid)-1 and valid[index-1:index+2].all():crossing=True
    result["crossing_curves_on_supported_grid"]=crossing
    result["matching_at_supported_curve_crossing"]=float(root) if crossing else np.nan
    contrast=probability_contrasts(m2["beta"],m2["covariance"],support["quantile_values"])
    for index,label in enumerate(["q10","q50","q90"]):
        result[label+"_matching_centered"]=support["quantile_values"][index]
        for field in ["p_south","p_north","difference_north_minus_south","difference_95_low_pointwise","difference_95_high_pointwise"]:
            result[label+"_"+field]=contrast[field].iloc[index]
    return result,grid
