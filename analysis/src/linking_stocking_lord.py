"""Directional characteristic-curve optimization, with explicit B-to-A direction."""
from __future__ import annotations

from dataclasses import dataclass,field
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit


@dataclass(frozen=True)
class LinkingConfig:
    scale_lower: float=.25
    scale_upper: float=4.
    shift_limit: float=3.
    max_iterations: int=500
    gtol: float=1e-10
    ftol: float=1e-15
    verified_gradient_tolerance: float=1e-7
    boundary_tolerance: float=1e-5
    starts: tuple=field(default_factory=lambda:((1.,0.),(.8,-.2),(1.25,.2),(1.,.5),(1.,-.5)))


def normal_grid(lower=-3.,upper=3.,points=601):
    if points<3 or upper<=lower:raise ValueError("Ordered sufficiently resolved grid required")
    theta=np.linspace(lower,upper,points)
    weights=np.exp(-theta*theta/2)/np.sqrt(2*np.pi)
    weights[[0,-1]]*=.5;weights*=theta[1]-theta[0];weights/=weights.sum()
    return theta,weights


def objective_gradient(parameters,d_a,a_a,d_b,a_b,theta,weights,method="Stocking-Lord"):
    """Normalized objective and derivative w.r.t. [log(scale), shift].

    Stocking–Lord divides anchor-TCC squared discrepancy by K² solely for
    optimization; Haebara divides the sum of item squared discrepancies by K.
    The raw objective is recovered in the result table. Both minimize the
    requested directional A-grid discrepancy with normalized reference weights.
    """
    d_a,a_a,d_b,a_b=[np.asarray(value,float) for value in [d_a,a_a,d_b,a_b]]
    theta=np.asarray(theta,float);weights=np.asarray(weights,float)
    if not (d_a.shape==a_a.shape==d_b.shape==a_b.shape) or d_a.ndim!=1 or not len(d_a):
        raise ValueError("Nonempty corresponding anchor parameter arrays required")
    if len(theta)!=len(weights) or not np.isfinite(weights).all() or (weights<0).any() or not np.isclose(weights.sum(),1):
        raise ValueError("Normalized nonnegative grid weights required")
    scale=np.exp(parameters[0]);shift=parameters[1]
    source_theta=(theta-shift)/scale
    pa=expit(d_a+theta[:,None]*a_a)
    pb=expit(d_b+source_theta[:,None]*a_b)
    derivative=np.stack([-a_b*source_theta[:,None],np.broadcast_to(-a_b/scale,pb.shape)],axis=2)
    derivative*=pb[:,:,None]*(1-pb[:,:,None])
    difference=pa-pb;k=len(d_a)
    if method=="Stocking-Lord":
        delta=difference.sum(axis=1)/k;dp=derivative.sum(axis=1)/k
        return float(weights@(delta*delta)),-2*np.sum(weights[:,None]*delta[:,None]*dp,axis=0)
    if method=="Haebara-directional-A":
        return float(np.sum(weights[:,None]*difference*difference)/k),-2*np.sum(weights[:,None,None]*difference[:,:,None]*derivative,axis=(0,1))/k
    raise ValueError("Unknown characteristic-curve method")


def fit_link(d_a,a_a,d_b,a_b,theta,weights,config=None,method="Stocking-Lord",starts=None):
    config=config or LinkingConfig()
    starts=config.starts if starts is None else starts
    lower,upper=np.log(config.scale_lower),np.log(config.scale_upper)
    bounds=[(lower,upper),(-config.shift_limit,config.shift_limit)]
    attempts=[];candidates=[]
    for index,(scale,shift) in enumerate(starts,1):
        if scale<=0:raise ValueError("Starting scales must be positive")
        initial=np.array([np.log(scale),shift])
        result=minimize(objective_gradient,initial,args=(d_a,a_a,d_b,a_b,theta,weights,method),jac=True,
                        method="L-BFGS-B",bounds=bounds,
                        options={"maxiter":config.max_iterations,"gtol":config.gtol,"ftol":config.ftol,"maxls":50})
        value,gradient=objective_gradient(result.x,d_a,a_a,d_b,a_b,theta,weights,method)
        boundary=bool(min(result.x[0]-lower,upper-result.x[0],config.shift_limit-abs(result.x[1]))<=config.boundary_tolerance)
        finite=bool(np.isfinite(value) and np.isfinite(result.x).all())
        first_order=bool(np.max(abs(gradient))<=config.verified_gradient_tolerance)
        valid=bool(result.success and finite and first_order and not boundary)
        row={"start":index,"start_scale":scale,"start_shift":shift,"optimizer_success":bool(result.success),
             "first_order_verified":first_order,"supported":valid,"boundary_hit":boundary,"iterations":int(result.nit),
             "normalized_objective":value,"maximum_absolute_gradient":float(np.max(abs(gradient))),
             "scale_A":float(np.exp(result.x[0])),"shift_B":float(result.x[1]),"message":str(result.message)}
        attempts.append(row);candidates.append((valid,value,result.x,row))
    if not candidates:raise ValueError("At least one deterministic start required")
    chosen=min([c for c in candidates if c[0]] or candidates,key=lambda candidate:candidate[1])
    valid,value,parameters,diagnostic=chosen
    scale,shift=float(np.exp(parameters[0])),float(parameters[1])
    factor=len(d_a)**2 if method=="Stocking-Lord" else len(d_a)
    accepted=[row for row in attempts if row["supported"]]
    return {"scale_A":scale,"shift_B":shift,"method":method,"supported":valid,
            "normalized_objective":value,"raw_objective":value*factor,
            "maximum_absolute_gradient":diagnostic["maximum_absolute_gradient"],"iterations":diagnostic["iterations"],
            "boundary_hit":diagnostic["boundary_hit"],"n_supported_starts":len(accepted),
            "supported_start_scale_range":max(row["scale_A"] for row in accepted)-min(row["scale_A"] for row in accepted) if accepted else np.nan,
            "supported_start_shift_range":max(row["shift_B"] for row in accepted)-min(row["shift_B"] for row in accepted) if accepted else np.nan,
            "attempts":pd.DataFrame(attempts)}


def stocking_lord(d_a,a_a,d_b,a_b,theta,weights,config=None,starts=None):
    return fit_link(d_a,a_a,d_b,a_b,theta,weights,config,"Stocking-Lord",starts)
