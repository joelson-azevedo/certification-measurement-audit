"""Hypothetical conditional-independent joint classifications and their integrals."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.special import expit
from .score_distribution import reference_quadrature
from .classification_probabilities import passing_curves

CELLS=["both_pass","first_pass_second_fail","first_fail_second_pass","both_fail"]


def joint_classification(first,second=None):
    p=np.asarray(first,float);q=p if second is None else np.asarray(second,float)
    if p.shape!=q.shape or not np.isfinite(p).all() or not np.isfinite(q).all() or min(p.min(),q.min()) < -1e-11 or max(p.max(),q.max())>1+1e-11:
        raise ValueError("Finite matching conditional pass probabilities in [0,1] required")
    pp=p*q;pf=p*(1-q);fp=(1-p)*q;ff=(1-p)*(1-q)
    return {"both_pass":pp,"first_pass_second_fail":pf,"first_fail_second_pass":fp,"both_fail":ff,
            "agreement":pp+ff,"disagreement":pf+fp,"first_pass":p,"second_pass":q}


def joint_statistics(cells):
    pp,pf,fp,ff=[float(cells[name]) for name in CELLS];first=pp+pf;second=pp+fp
    agreement=pp+ff;chance=first*second+(1-first)*(1-second);denominator=1-chance
    return {**dict(zip(CELLS,[pp,pf,fp,ff])),"agreement":agreement,"disagreement":pf+fp,
            "first_pass":first,"second_pass":second,"chance_agreement_from_marginals":chance,
            "positive_agreement":2*pp/(2*pp+pf+fp) if 2*pp+pf+fp>1e-14 else np.nan,
            "negative_agreement":2*ff/(2*ff+pf+fp) if 2*ff+pf+fp>1e-14 else np.nan,
            "kappa":(agreement-chance)/denominator if denominator>1e-14 else np.nan,
            "kappa_status":"available: model-implied joint distribution" if denominator>1e-14 else "undefined: degenerate marginal agreement",
            "normalization_error":abs(pp+pf+fp+ff-1)}


def integrate_joint(first,second,weights):
    weights=np.asarray(weights,float);joint=joint_classification(first,second)
    if weights.shape!=np.asarray(first).shape or not np.isfinite(weights).all() or (weights<0).any() or not np.isclose(weights.sum(),1,rtol=0,atol=1e-12):
        raise ValueError("Common normalized reference weights required")
    return joint_statistics({name:weights@joint[name] for name in CELLS})


def integrated_models(da,aa,db,ab,scale,shift,reference,scenarios,resolutions,cell_tolerance=1e-8,kappa_tolerance=1e-7):
    previous=None;logs=[];cuts=sorted({int(s["cut_B"]) for s in scenarios})
    for nodes in resolutions:
        theta,w=reference_quadrature(reference,nodes,scale,shift)
        a,pmfA,_=passing_curves(theta,da,aa,[42]);b,pmfB,_=passing_curves(theta,db,ab,cuts)
        rows=[]
        rows.append({"pair_type":"same_form","form":"A","scenario":"A42","cut_first":42,"cut_second":42,**integrate_joint(a[42],a[42],w)})
        for cut in cuts:rows.append({"pair_type":"same_form","form":"B","scenario":f"B{cut}","cut_first":cut,"cut_second":cut,**integrate_joint(b[cut],b[cut],w)})
        for scenario in scenarios:
            cut=int(scenario["cut_B"]);rows.append({"pair_type":"cross_form","form":"A/B","scenario":scenario["scenario"],"cut_first":42,"cut_second":cut,**integrate_joint(a[42],b[cut],w)})
        frame=pd.DataFrame(rows)
        if previous is not None:
            old_nodes,old=previous
            cell_error=float(abs(frame[CELLS].to_numpy()-old[CELLS].to_numpy()).max())
            current_kappa,old_kappa=frame.kappa.to_numpy(),old.kappa.to_numpy()
            matching_status=np.array_equal(np.isfinite(current_kappa),np.isfinite(old_kappa))
            finite=np.isfinite(current_kappa)&np.isfinite(old_kappa)
            kappa_error=float(max(abs(current_kappa[finite]-old_kappa[finite]))) if matching_status and finite.any() else (0. if matching_status else np.inf)
            passed=cell_error<=cell_tolerance and kappa_error<=kappa_tolerance
            logs.append({"lower_nodes":old_nodes,"higher_nodes":nodes,"maximum_joint_cell_difference":cell_error,"maximum_kappa_difference":kappa_error,"converged":passed})
            if passed:return frame,{"theta":theta,"weights":w,"A_tails":a,"B_tails":b,"conditional_pmf_A":pmfA,"conditional_pmf_B":pmfB},pd.DataFrame(logs)
        previous=(nodes,frame)
    raise ValueError("Classification integration did not meet quadrature tolerances: "+str(logs))


def monte_carlo_fixture(draws=500000,seed=20261029,multiplier=6.,floor=1e-4):
    """New controlled parameters; raw response simulation, not a pass-coin shortcut."""
    da=np.array([-.6,.2,1.,-.3,.7,.4]);aa=np.array([1.2,.8,-.7,1.4,1.,.5])
    db=np.array([-1.,-.1,.6,-.7,.4,.2]);ab=np.array([1.,1.1,-.6,1.2,.9,.7])
    theta,w=reference_quadrature("A_standard_normal",321);pa=passing_curves(theta,da,aa,[4])[0][4];pb=passing_curves(theta,db,ab,[3])[0][3]
    rng=np.random.default_rng(np.random.SeedSequence([seed,0]));latent=rng.standard_normal(draws)
    classifications=[]
    for code,d,a,cut in [(1,da,aa,4),(2,da,aa,4),(3,db,ab,3)]:
        generator=np.random.default_rng(np.random.SeedSequence([seed,code]))
        responses=generator.random((draws,len(d)))<expit(d+latent[:,None]*a)
        classifications.append(responses.sum(axis=1)>=cut)
    first,repeat,other=classifications;rows=[]
    for name,second,expected in [("same_form",repeat,integrate_joint(pa,pa,w)),("cross_form",other,integrate_joint(pa,pb,w))]:
        observations={"both_pass":first&second,"first_pass_second_fail":first&~second,"first_fail_second_pass":~first&second,
                      "both_fail":~first&~second,"agreement":first==second,"disagreement":first!=second}
        for statistic,values in observations.items():
            theory=float(expected[statistic]);empirical=float(values.mean());se=float(np.sqrt(theory*(1-theory)/draws));tolerance=max(multiplier*se,floor)
            rows.append({"fixture":name,"statistic":statistic,"theoretical":theory,"Monte_Carlo_estimate":empirical,
                         "Monte_Carlo_SE":se,"absolute_difference":abs(empirical-theory),"tolerance":tolerance,"passed":abs(empirical-theory)<=tolerance,
                         "n_independent_latent_response_sets":draws,"seed":seed,"source":"independent controlled Stage 9 fixture; not certification-data generating truth"})
    return pd.DataFrame(rows)
