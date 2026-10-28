"""Reuse saved paired parameter/link draws; no calibration or linking refits."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import binom
from .equating_decision_impact import aggregate_impact
from .irt_observed_score_equating import equate_observed,mid_cdf
from .irt_true_score_equating import tcc,tcc_derivative
from .score_distribution import marginal_distribution,reference_quadrature


def converged_pair(da,aa,db,ab,reference,scale,shift,config):
    """Use the higher member of the first passing successive GH comparison."""
    previous=None;log=[];cut_scores=np.arange(config["cut_region_A"][0],config["cut_region_A"][1]+1)
    for nodes in config["quadrature"]["resolutions"]:
        theta,weights=reference_quadrature(reference,nodes,scale,shift)
        pair={"A":marginal_distribution(da,aa,theta,weights),"B":marginal_distribution(db,ab,theta,weights)}
        for form,result in pair.items():
            d=result["diagnostics"];v=config["validation"]
            if not (d["all_probabilities_finite"] and d["minimum_probability"]>=0 and d["maximum_normalization_error"]<=v["probability_tolerance"]
                    and d["marginal_normalization_error"]<=v["probability_tolerance"] and d["maximum_mean_TCC_error"]<=v["conditional_mean_tolerance"]
                    and d["maximum_variance_Bernoulli_error"]<=v["conditional_variance_tolerance"] and d["marginal_mean_integral_error"]<=v["conditional_mean_tolerance"]
                    and d["marginal_variance_total_variance_error"]<=v["conditional_variance_tolerance"] and d["cdf_minimum_increment"]>=0):
                raise ValueError("Conditional/marginal probability validation failed: "+form)
        forward=equate_observed(np.arange(len(da)+1),pair["B"]["pmf"],pair["A"]["pmf"])["equated"]
        inverse=equate_observed(cut_scores,pair["A"]["pmf"],pair["B"]["pmf"])["equated"]
        if previous is not None:
            low,old,oldforward,oldinverse=previous
            percentile=mid_cdf(pair["B"]["pmf"]);central=(percentile>=.025)&(percentile<=.975)
            central&=np.isfinite(forward)&np.isfinite(oldforward)
            difference=float(max(abs(pair[f]["pmf"]-old[f]["pmf"]).max() for f in ["A","B"]))
            mapping=float(max(abs(forward[central]-oldforward[central]))) if central.any() else np.inf
            cut=float(max(abs(inverse-oldinverse))) if np.isfinite(inverse).all() and np.isfinite(oldinverse).all() else np.inf
            passed=difference<=config["quadrature"]["maximum_pmf_difference"] and mapping<=config["quadrature"]["maximum_central_mapping_difference"] and cut<=config["quadrature"]["maximum_cut_location_difference"]
            log.append({"lower_nodes":low,"higher_nodes":nodes,"maximum_pmf_difference":difference,"maximum_central_mapping_difference":mapping,
                        "maximum_inverse_cut_region_mapping_difference":cut,"converged":bool(passed)})
            if passed:
                # No probability clipping, rescaling or empirical smoothing is used.
                return pair,nodes,pd.DataFrame(log)
        previous=(nodes,pair,forward,inverse)
    raise ValueError("Score integration failed declared GH convergence criteria: "+str(log))


def evaluate_replicates(arrays,links,method,config,raw_B,recorded_B,attempts=None,progress=None):
    selected=links.loc[links.method.eq(method)].copy()
    if attempts is not None:selected=selected.loc[selected.attempt.isin(attempts)]
    rows=[];diagnostics=[];pmfs_A=[];pmfs_B=[];map_arrays=[];inverse_arrays=[];true_arrays=[];attempt_numbers=[]
    cut_scores=np.arange(config["cut_region_A"][0],config["cut_region_A"][1]+1)
    for count,row in enumerate(selected.itertuples(),1):
        attempt=row.attempt;failure="";valid=False
        try:
            if not row.valid_paired_draw or not row.link_supported or row.scale_A<=0:raise ValueError("Unsupported stored Stage 7 paired draw/link")
            va=arrays["parameters_A"][attempt];vb=arrays["parameters_B"][attempt]
            k=len(va)//2;da,aa=va[:k],va[k:];ab=vb[k:]/row.scale_A;db=vb[:k]-ab*row.shift_B
            pair,nodes,convergence=converged_pair(da,aa,db,ab,config["primary_reference_distribution"],row.scale_A,row.shift_B,config)
            # Signed derivative checks over declared finite theta scan; never flip individual slopes.
            scan=np.linspace(*config["true_score"]["theta_limits"],config["true_score"]["root_scan_points"])
            if not ((tcc_derivative(scan,da,aa)>0).all() and (tcc_derivative(scan,db,ab)>0).all()):raise ValueError("Replicate TCC monotonicity not supported on declared scan")
            mapping=equate_observed(np.arange(k+1),pair["B"]["pmf"],pair["A"]["pmf"])["equated"]
            inverse=equate_observed(cut_scores,pair["A"]["pmf"],pair["B"]["pmf"])["equated"]
            true=np.array([float(tcc(brentq(lambda z:float(tcc(z,da,aa)-score),*config["true_score"]["theta_limits"],xtol=config["true_score"]["root_xtol"]),db,ab)) for score in cut_scores])
            if not np.isfinite(inverse).all():raise ValueError("Replicate cut-region correspondence unavailable")
            impact,_,_=aggregate_impact(raw_B,recorded_B,pair["B"]["pmf"],pair["A"]["pmf"],config["counterfactual_conventions"])
            for kind,scores,values in [("B_to_A_observed",np.arange(k+1),mapping),("A_to_B_observed_inverse",cut_scores,inverse),("A_to_B_true",cut_scores,true)]:
                rows.extend({"attempt":attempt,"method":method,"kind":kind,"score":float(s),"convention":"not applicable","value":float(value)} for s,value in zip(scores,values))
            for record in impact.to_dict("records"):
                for metric in ["counterfactual_pass","upward_changes","downward_changes","minimum_B_integer_passing_research_rule"]:
                    rows.append({"attempt":attempt,"method":method,"kind":metric,"score":42.,"convention":record["convention"],"value":record[metric]})
            pmfs_A.append(pair["A"]["pmf"]);pmfs_B.append(pair["B"]["pmf"]);map_arrays.append(mapping);inverse_arrays.append(inverse);true_arrays.append(true);attempt_numbers.append(attempt)
            valid=True
            diagnostic={"selected_nodes":nodes,**convergence.iloc[-1].to_dict(),
                        "max_conditional_mean_error":max(pair[f]["diagnostics"]["maximum_mean_TCC_error"] for f in ["A","B"]),
                        "max_conditional_variance_error":max(pair[f]["diagnostics"]["maximum_variance_Bernoulli_error"] for f in ["A","B"])}
        except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:
            failure=str(exc);diagnostic={}
        diagnostics.append({"attempt":attempt,"method":method,"valid_stage8_replicate":valid,"failure":failure,**diagnostic})
        if progress and (count%100==0 or count==len(selected)):progress(f"{method}: {count}/{len(selected)} stored paired draws processed; {sum(d['valid_stage8_replicate'] for d in diagnostics)} valid")
    return {"values":pd.DataFrame(rows,columns=["attempt","method","kind","score","convention","value"]),"diagnostics":pd.DataFrame(diagnostics),"arrays":{
        "attempts":np.asarray(attempt_numbers),"marginal_pmf_A":np.asarray(pmfs_A),"marginal_pmf_B":np.asarray(pmfs_B),
        "B_to_A_observed":np.asarray(map_arrays),"A_cut_scores":cut_scores,"A_to_B_observed_inverse":np.asarray(inverse_arrays),"A_to_B_true":np.asarray(true_arrays)}}


def interval_summary(values,attempted,minimum_valid=100,minimum_fraction=.9):
    if attempted<=0:raise ValueError("A positive attempted replicate count is required")
    rows=[]
    for keys,group in values.groupby(["method","kind","score","convention"],sort=False):
        vector=group.value.to_numpy(float);vector=vector[np.isfinite(vector)];n=len(vector)
        supported=n>=minimum_valid and n/attempted>=minimum_fraction
        row=dict(zip(["method","kind","score","convention"],keys))
        row.update(n_valid=n,n_attempted=attempted,uncertainty_status="available: pointwise fixed-S18/model/reference percentile interval" if supported else "unavailable: insufficient finite supported replicates",
                   mean=float(vector.mean()) if n else np.nan,se=float(vector.std(ddof=1)) if n>1 else np.nan,
                   low_95=float(np.quantile(vector,.025)) if supported else np.nan,high_95=float(np.quantile(vector,.975)) if supported else np.nan,
                   mean_monte_carlo_se=float(vector.std(ddof=1)/np.sqrt(n)) if n>1 else np.nan)
        if n:
            sort=np.sort(vector)
            for q,label in [(.025,"lower_endpoint"),(.975,"upper_endpoint")]:
                lower=max(0,int(binom.ppf(.025,n,q))-1);upper=min(n-1,int(binom.ppf(.975,n,q)))
                row[label+"_mc95_low"]=float(sort[lower]);row[label+"_mc95_high"]=float(sort[upper])
        rows.append(row)
    return pd.DataFrame(rows)
