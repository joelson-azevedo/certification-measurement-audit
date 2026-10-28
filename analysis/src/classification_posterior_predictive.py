"""Response-vector posterior integration under each source model's declared prior."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from .score_distribution import normal_quadrature
from .classification_probabilities import passing_curves


def response_posterior(x,d,a,theta,weights):
    x=np.asarray(x,float);d,a=np.asarray(d,float),np.asarray(a,float);theta,weights=np.asarray(theta,float),np.asarray(weights,float)
    if x.ndim!=2 or x.shape[1]!=len(d) or a.shape!=d.shape or np.isinf(x).any() or not np.isin(x[np.isfinite(x)],[0,1]).all():raise ValueError("Binary response vectors with NaN omissions required")
    if theta.ndim!=1 or weights.shape!=theta.shape or (weights<0).any() or not np.isclose(weights.sum(),1,rtol=0,atol=1e-12):raise ValueError("Declared normalized prior quadrature required")
    observed=np.isfinite(x).astype(float);correct=np.nan_to_num(x,nan=0.)
    eta=d+theta[:,None]*a;logp=-np.logaddexp(0,-eta);logq=-np.logaddexp(0,eta)
    likelihood=correct@logp.T+(observed-correct)@logq.T
    with np.errstate(divide="ignore"):joint=likelihood+np.log(weights)
    normalizer=logsumexp(joint,axis=1);posterior=np.exp(joint-normalizer[:,None])
    if not np.isfinite(posterior).all() or not np.allclose(posterior.sum(axis=1),1,rtol=0,atol=1e-11):raise ValueError("Posterior normalization failed")
    mean=posterior@theta;variance=posterior@(theta*theta)-mean*mean
    if variance.min() < -1e-11:raise ValueError("Negative posterior variance")
    return {"posterior":posterior,"mean":mean,"sd":np.sqrt(np.maximum(variance,0)),"log_marginal_likelihood":normalizer}


def predictive_moments(posterior,conditional_pass):
    posterior=np.asarray(posterior,float);q=np.asarray(conditional_pass,float)
    if posterior.ndim!=2 or q.ndim!=1 or posterior.shape[1]!=len(q):raise ValueError("Matching posterior and conditional passing curve required")
    mean=posterior@q;second=posterior@(q*q);between=second-mean*mean;within=posterior@(q*(1-q));future=mean*(1-mean)
    if min(between.min(),within.min(),future.min()) < -1e-10:raise ValueError("Predictive variance decomposition invalid")
    return {"future_pass_probability":mean,"posterior_sd_conditional_pass_probability":np.sqrt(np.maximum(between,0)),
            "within_theta_future_Bernoulli_variance":within,"between_theta_pass_probability_variance":between,
            "future_Bernoulli_variance":future,"variance_decomposition_error":abs(within+between-future)}


def predict_candidates(x,source_form,source_d,source_a,dA,aA,dB_linked,aB_linked,scale,shift,cuts_B,resolutions,max_predictive_difference=1e-5,max_eap_difference=1e-5):
    previous=None;logs=[]
    for nodes in resolutions:
        own,w=normal_quadrature(nodes);posterior=response_posterior(x,source_d,source_a,own,w)
        theta_A=own if source_form=="A" else scale*own+shift
        curvesA=passing_curves(theta_A,dA,aA,[42])[0];curvesB=passing_curves(theta_A,dB_linked,aB_linked,cuts_B)[0]
        target={"A42":curvesA[42],**{f"B{cut}":curve for cut,curve in curvesB.items()}}
        mean=np.column_stack([posterior["posterior"]@curve for curve in target.values()])
        if previous is not None:
            low,oldmean,oldposterior=previous
            difference=float(abs(mean-oldmean).max());eap_difference=float(abs(posterior["mean"]-oldposterior["mean"]).max())
            valid=difference<=max_predictive_difference and eap_difference<=max_eap_difference
            logs.append({"source_form":source_form,"lower_nodes":low,"higher_nodes":nodes,"maximum_candidate_future_pass_difference":difference,"maximum_candidate_EAP_difference":eap_difference,"converged":valid})
            if valid:
                moments={name:predictive_moments(posterior["posterior"],curve) for name,curve in target.items()}
                # Plug-in results are diagnostic comparisons, not the prediction estimator.
                eap_A=posterior["mean"] if source_form=="A" else scale*posterior["mean"]+shift
                plugA=passing_curves(eap_A,dA,aA,[42])[0];plugB=passing_curves(eap_A,dB_linked,aB_linked,cuts_B)[0]
                plugins={"A42":plugA[42],**{f"B{cut}":p for cut,p in plugB.items()}}
                return {"nodes":nodes,"theta_source":own,"theta_A":theta_A,"weights":w,"posterior":posterior,
                        "conditional_pass":target,"moments":moments,"EAP_plugins":plugins,"convergence":pd.DataFrame(logs)}
        previous=(nodes,mean,posterior)
    raise ValueError("Posterior-predictive quadrature did not converge: "+str(logs))


def predictive_score_groups(frame):
    rows=[]
    for keys,group in frame.groupby(["source_form","original_raw_score","target_form","future_cut","target_role"],sort=False):
        source,score,target,cut,role=keys;prediction=group.future_pass_probability.to_numpy();agreement=group.same_as_original_recorded_decision_probability.to_numpy()
        row={"source_form":source,"original_raw_score":score,"target_form":target,"future_cut":cut,"target_role":role,
             "n_first_attempts":len(group),"original_recorded_decision":group.original_recorded_decision.iloc[0],"mean_future_pass_probability":float(prediction.mean()),
             "sd_across_candidates_future_pass_probability":float(prediction.std(ddof=1)) if len(group)>1 else np.nan,
             "mean_same_original_decision_probability":float(agreement.mean()),"mean_posterior_sd_conditional_pass_probability":float(group.posterior_sd_conditional_pass_probability.mean()),
             "mean_within_theta_future_Bernoulli_variance":float(group.within_theta_future_Bernoulli_variance.mean()),
             "mean_between_theta_pass_probability_variance":float(group.between_theta_pass_probability_variance.mean()),
             "maximum_posterior_prediction_minus_EAP_plugin_abs_difference":float(group.posterior_prediction_minus_EAP_plugin.abs().max())}
        for values,label in [(prediction,"future_pass"),(agreement,"same_original_decision")]:
            for q,name in [(0,"min"),(.1,"p10"),(.5,"median"),(.9,"p90"),(1,"max")]:row[label+"_"+name]=float(np.quantile(values,q))
        rows.append(row)
    return pd.DataFrame(rows)
