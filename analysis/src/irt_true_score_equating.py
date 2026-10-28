"""Expected-score correspondence, with signed derivatives and explicit root checks."""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from .score_distribution import response_probabilities


def tcc(theta,d,a):return response_probabilities(theta,d,a).sum(axis=-1)


def tcc_derivative(theta,d,a):
    p=response_probabilities(theta,d,a)
    return (np.asarray(a)*p*(1-p)).sum(axis=-1)


def tcc_review(d,a,limits=(-6,6),points=12001):
    theta=np.linspace(*limits,points);score=tcc(theta,d,a);derivative=tcc_derivative(theta,d,a)
    return {"theta":theta,"tcc":score,"derivative":derivative,
            "strictly_increasing_on_checked_grid":bool((np.diff(score)>0).all() and (derivative>0).all()),
            "minimum_signed_derivative":float(derivative.min()),"minimum_score_increment":float(np.diff(score).min()),
            "theta_lower":limits[0],"theta_upper":limits[1],"score_lower":float(score[0]),"score_upper":float(score[-1]),
            "scope":"dense finite-grid monotonicity/derivative check, not proof of global monotonicity outside the scan"}


def score_roots(score,d,a,limits=(-6,6),points=2401,xtol=1e-12):
    theta=np.linspace(*limits,points);difference=tcc(theta,d,a)-score;roots=[]
    roots.extend(theta[difference==0].tolist())
    for index in np.flatnonzero(difference[:-1]*difference[1:]<0):
        roots.append(brentq(lambda z:float(tcc(z,d,a)-score),theta[index],theta[index+1],xtol=xtol))
    unique=[]
    for root in sorted(roots):
        if not unique or abs(root-unique[-1])>10*xtol:unique.append(float(root))
    return unique


def true_score_correspondence(source_scores,d_source,a_source,d_reference,a_reference,overlap,limits=(-6,6),scan_points=2401,xtol=1e-12):
    rows=[]
    for score in np.atleast_1d(np.asarray(source_scores,float)):
        roots=score_roots(score,d_source,a_source,limits,scan_points,xtol)
        row={"source_expected_score":float(score),"n_roots_in_scan":len(roots),"theta":np.nan,"reference_expected_score":np.nan,
             "source_signed_TCC_derivative":np.nan,"root_score_error":np.nan,"status":"unavailable: no root in declared theta scan"}
        if len(roots)>1:row["status"]="unavailable: multiple roots in declared theta scan"
        elif len(roots)==1:
            location=roots[0];derivative=float(tcc_derivative(location,d_source,a_source))
            row.update(theta=location,reference_expected_score=float(tcc(location,d_reference,a_reference)),
                       source_signed_TCC_derivative=derivative,root_score_error=float(abs(tcc(location,d_source,a_source)-score)),
                       status="supported central research region" if overlap[0]<=location<=overlap[1] and derivative>1e-6 else "weakly supported tail/model extrapolation")
        rows.append(row)
    return rows
