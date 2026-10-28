"""Paired stored-parameter propagation; no model fitting and no new random draws."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import binom
from .classification_probabilities import passing_curves
from .classification_posterior_predictive import predict_candidates
from .decision_consistency import integrated_models,joint_classification,joint_statistics
from .score_distribution import response_probabilities

STATISTICS=["agreement","disagreement","both_pass","both_fail","first_pass_second_fail","first_fail_second_pass",
            "first_pass","second_pass","positive_agreement","negative_agreement","kappa","chance_agreement_from_marginals",
            "agreement_change_from_O","agreement_change_from_R39","agreement_change_from_R40","second_pass_change_from_O",
            "first_pass_second_fail_change_from_O","first_fail_second_pass_change_from_O","selected_cut_B"]


def replicate_boundary(stage8_values,attempt,convention):
    group=stage8_values.loc[stage8_values.attempt.eq(attempt)&stage8_values.method.eq("Stocking-Lord")&
                            stage8_values.kind.eq("minimum_B_integer_passing_research_rule")&stage8_values.convention.eq(convention)]
    if len(group)!=1 or not np.isfinite(group.value.iloc[0]):raise ValueError("Missing/nonunique attempt-matched Stage 8 research boundary")
    cut=float(group.value.iloc[0])
    if not cut.is_integer() or not 0<=cut<=60:raise ValueError("Unsupported Stage 8 integer boundary")
    return int(cut)


def propagate_replicates(arrays,links,stage8_values,config,fixed_theta,nearcut_candidates=None,attempts=None,progress=None):
    selected=links.loc[links.method.eq("Stocking-Lord")].copy()
    if attempts is not None:selected=selected.loc[selected.attempt.isin(attempts)]
    rows=[];diagnostics=[];convergence=[];posterior_logs=[];band_rows=[];rerooted=[]
    pass_arrays=[];valid_attempts=[];cuts=np.array([42,39,40,38]);root_grid=np.asarray(fixed_theta,float)
    for count,row in enumerate(selected.itertuples(),1):
        attempt=int(row.attempt);failure="";valid=False;band_valid=None
        try:
            if not row.valid_paired_draw or not row.link_supported or row.scale_A<=0:raise ValueError("Unsupported stored paired link")
            va=np.asarray(arrays["parameters_A"][attempt]);vb=np.asarray(arrays["parameters_B"][attempt]);k=len(va)//2
            da,aa=va[:k],va[k:];db=vb[:k]-vb[k:]*row.shift_B/row.scale_A;ab=vb[k:]/row.scale_A
            scenarios=[dict(s) for s in config["scenarios"]];dynamic={}
            for convention in config["uncertainty"]["dynamic_boundaries"]:
                name="D_continuous" if convention=="continuous_ge42" else "D_nearest"
                dynamic[name]=replicate_boundary(stage8_values,attempt,convention)
                scenarios.append({"scenario":name,"cut_A":42,"cut_B":dynamic[name],"role":"attempt-matched Stage 8 boundary; research only"})
            q=config["quadrature"]
            population,details,log=integrated_models(da,aa,db,ab,row.scale_A,row.shift_B,config["primary_reference"],scenarios,q["resolutions"],q["maximum_joint_cell_difference"],q["maximum_kappa_difference"])
            base=population.loc[population.pair_type.eq("cross_form")&population.scenario.eq("O")].iloc[0]
            r39=population.loc[population.pair_type.eq("cross_form")&population.scenario.eq("R39")].iloc[0]
            r40=population.loc[population.pair_type.eq("cross_form")&population.scenario.eq("R40")].iloc[0]
            extra=[]
            for name,cut in dynamic.items():
                same=population.loc[population.pair_type.eq("same_form")&population.form.eq("B")&population.cut_second.eq(cut)].iloc[0].to_dict()
                same["scenario"]="SAME_"+name;extra.append(same)
            if extra:population=pd.concat([population,pd.DataFrame(extra)],ignore_index=True)
            for record in population.to_dict("records"):
                record["attempt"]=attempt;record["selected_cut_B"]=record["cut_second"] if record["form"]!="A" else np.nan
                if record["pair_type"]=="cross_form":
                    record.update(agreement_change_from_O=record["agreement"]-base.agreement,
                                  agreement_change_from_R39=record["agreement"]-r39.agreement,agreement_change_from_R40=record["agreement"]-r40.agreement,
                                  second_pass_change_from_O=record["second_pass"]-base.second_pass,
                                  first_pass_second_fail_change_from_O=record["first_pass_second_fail"]-base.first_pass_second_fail,
                                  first_fail_second_pass_change_from_O=record["first_fail_second_pass"]-base.first_fail_second_pass)
                rows.append(record)
            log.insert(0,"attempt",attempt);convergence.append(log)
            pa=passing_curves(root_grid,da,aa,[42])[0][42];pb=passing_curves(root_grid,db,ab,cuts)[0]
            pass_arrays.append(np.column_stack([pa]+[pb[int(c)] for c in cuts]));valid_attempts.append(attempt)
            root=brentq(lambda theta:float(response_probabilities(theta,da,aa).sum()-42),-6,6,xtol=1e-12)
            ar=passing_curves(np.array([root]),da,aa,[42])[0][42][0]
            br=passing_curves(np.array([root]),db,ab,sorted(set(int(s["cut_B"]) for s in scenarios)))[0]
            for scenario in scenarios:
                b=float(br[int(scenario["cut_B"])][0]);j=joint_statistics({name:float(value) for name,value in joint_classification(ar,b).items() if name in ["both_pass","both_fail","first_pass_second_fail","first_fail_second_pass"]})
                rerooted.append({"attempt":attempt,"scenario":scenario["scenario"],"theta_A":root,**j})
            valid=True
            if nearcut_candidates is not None:
                band_valid=True
                try:
                    for form,data in nearcut_candidates.items():
                        own_d,own_a=(da,aa) if form=="A" else (vb[:k],vb[k:])
                        pc=config["posterior"]
                        result=predict_candidates(data["x"],form,own_d,own_a,da,aa,db,ab,row.scale_A,row.shift_B,cuts,
                                                  pc["resolutions"],pc["maximum_future_pass_difference"],pc["maximum_EAP_difference"])
                        plog=result["convergence"].copy();plog.insert(0,"attempt",attempt);posterior_logs.append(plog)
                        raw=data["raw"];recorded=raw>=42
                        for target,moment in result["moments"].items():
                            passing=moment["future_pass_probability"];agreement=np.where(recorded,passing,1-passing)
                            for score in np.unique(raw):
                                selected_score=raw==score
                                band_rows.append({"attempt":attempt,"source_form":form,"original_raw_score":int(score),"target":target,
                                                  "mean_future_pass_probability":float(passing[selected_score].mean()),
                                                  "mean_same_original_decision_probability":float(agreement[selected_score].mean()),
                                                  "n_original_candidates_fixed":int(selected_score.sum())})
                except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:
                    band_valid=False;failure="posterior-band prediction: "+str(exc)
                    band_rows=[r for r in band_rows if r["attempt"]!=attempt]
            diagnostic={"selected_nodes":int(log.higher_nodes.iloc[-1]),"maximum_joint_cell_difference":float(log.maximum_joint_cell_difference.iloc[-1]),
                        "maximum_kappa_difference":float(log.maximum_kappa_difference.iloc[-1]),"dynamic_continuous_B_cut":dynamic.get("D_continuous"),"dynamic_nearest_B_cut":dynamic.get("D_nearest")}
        except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:
            failure=str(exc);diagnostic={}
        diagnostics.append({"attempt":attempt,"valid_population_replicate":valid,"posterior_band_valid":band_valid,"failure":failure,**diagnostic})
        if progress and (count%50==0 or count==len(selected)):progress(f"{count}/{len(selected)} paired draws; {sum(d['valid_population_replicate'] for d in diagnostics)} valid population calculations")
    valid_ids={d["attempt"] for d in diagnostics if d["valid_population_replicate"]}
    keep=[j for j,attempt in enumerate(valid_attempts) if attempt in valid_ids]
    rows=[r for r in rows if r["attempt"] in valid_ids];rerooted=[r for r in rerooted if r["attempt"] in valid_ids]
    valid_attempts=[valid_attempts[j] for j in keep];pass_arrays=[pass_arrays[j] for j in keep]
    return {"population":pd.DataFrame(rows),"diagnostics":pd.DataFrame(diagnostics),"convergence":pd.concat(convergence,ignore_index=True) if convergence else pd.DataFrame(),
            "bands":pd.DataFrame(band_rows),"posterior_convergence":pd.concat(posterior_logs,ignore_index=True) if posterior_logs else pd.DataFrame(),
            "rerooted":pd.DataFrame(rerooted),"arrays":{"attempts":np.asarray(valid_attempts),"theta_A_fixed":root_grid,"target_labels":np.array(["A42","B42","B39","B40","B38"]),"conditional_pass":np.asarray(pass_arrays)}}


def empirical_intervals(frame,keys,statistics,attempted,minimum_valid=100,minimum_fraction=.9):
    rows=[]
    if attempted<=0:raise ValueError("Positive attempted count required")
    if frame.empty:return pd.DataFrame()
    for group_key,group in frame.groupby(keys,sort=False,dropna=False):
        if not isinstance(group_key,tuple):group_key=(group_key,)
        for statistic in statistics:
            if statistic not in group:continue
            vector=group[statistic].to_numpy(float);vector=vector[np.isfinite(vector)];n=len(vector)
            supported=n>=minimum_valid and n/attempted>=minimum_fraction
            row={**dict(zip(keys,group_key)),"statistic":statistic,"valid_replicates":n,"attempted_replicates":attempted,
                 "interval_status":"available: conditional pointwise percentile interval" if supported else "unavailable: insufficient finite replicates",
                 "mean":float(vector.mean()) if n else np.nan,"se":float(vector.std(ddof=1)) if n>1 else np.nan,
                 "low_95":float(np.quantile(vector,.025)) if supported else np.nan,"high_95":float(np.quantile(vector,.975)) if supported else np.nan,
                 "mean_monte_carlo_se":float(vector.std(ddof=1)/np.sqrt(n)) if n>1 else np.nan}
            if n:
                ordered=np.sort(vector)
                for quantile,label in [(.025,"lower_endpoint"),(.975,"upper_endpoint")]:
                    low=max(0,int(binom.ppf(.025,n,quantile))-1);high=min(n-1,int(binom.ppf(.975,n,quantile)))
                    row[label+"_mc95_low"]=float(ordered[low]);row[label+"_mc95_high"]=float(ordered[high])
            rows.append(row)
    return pd.DataFrame(rows)


def conditional_interval_summary(arrays,attempted,minimum_valid=100,minimum_fraction=.9):
    passing=arrays["conditional_pass"];theta=arrays["theta_A_fixed"];targets=arrays["target_labels"].tolist();rows=[]
    if not len(passing):return pd.DataFrame()
    for j,target in enumerate(targets):
        for index,location in enumerate(theta):
            values=passing[:,index,j]
            for statistic,vector in [("passing",values),("same_form_agreement",values*values+(1-values)**2)]:
                rows.append({"theta_A":location,"target_or_scenario":target,"statistic":statistic,"valid_replicates":len(vector),
                             "low_95":float(np.quantile(vector,.025)) if len(vector)>=minimum_valid and len(vector)/attempted>=minimum_fraction else np.nan,
                             "high_95":float(np.quantile(vector,.975)) if len(vector)>=minimum_valid and len(vector)/attempted>=minimum_fraction else np.nan})
    for scenario,target in [("O","B42"),("R39","B39"),("R40","B40"),("R38","B38")]:
        a=passing[:,:,0];b=passing[:,:,targets.index(target)]
        values={"cross_form_agreement":a*b+(1-a)*(1-b),"A_pass_B_fail":a*(1-b),"A_fail_B_pass":(1-a)*b}
        for statistic,v in values.items():
            for index,location in enumerate(theta):
                rows.append({"theta_A":location,"target_or_scenario":scenario,"statistic":statistic,"valid_replicates":len(v),
                             "low_95":float(np.quantile(v[:,index],.025)) if len(v)>=minimum_valid and len(v)/attempted>=minimum_fraction else np.nan,
                             "high_95":float(np.quantile(v[:,index],.975)) if len(v)>=minimum_valid and len(v)/attempted>=minimum_fraction else np.nan})
    return pd.DataFrame(rows)
