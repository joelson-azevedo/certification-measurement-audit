"""Candidate-level recalibration bootstrap and full joint-covariance propagation."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import binom

from .irt_estimation import fit_mml
from .linking_stocking_lord import fit_link
from .linking_transformations import probabilities,transform_parameters


def covariance_diagnostics(model):
    covariance=np.asarray(model["covariance"],float)
    eigen=np.linalg.eigvalsh(covariance)
    off=covariance-np.diag(np.diag(covariance))
    supported=bool(np.isfinite(covariance).all() and np.allclose(covariance,covariance.T,atol=1e-10)
                   and eigen.min()>0 and np.linalg.cond(covariance)<1e10)
    return {"form":model.get("form","fixture"),"n_parameters":len(covariance),"covariance_supported":supported,
            "minimum_eigenvalue":float(eigen.min()),"maximum_eigenvalue":float(eigen.max()),
            "condition_number":float(np.linalg.cond(covariance)),"off_diagonal_frobenius_norm":float(np.linalg.norm(off)),
            "maximum_information_covariance_identity_error":float(np.max(abs(model["information"]@covariance-np.eye(len(covariance))))) if "information" in model else np.nan,
            "parameter_order":"all intercepts followed by all slopes; full within-form covariance retained",
            "assumption":"asymptotic joint normal approximation; separately calibrated form draws independent"}


def joint_parameter_draw(model,seed,attempt,form_code):
    rng=np.random.default_rng(np.random.SeedSequence([seed,attempt,form_code]))
    vector=np.r_[model["d"],model["a"]]
    return vector+np.linalg.cholesky(model["covariance"])@rng.standard_normal(len(vector))


def experiment(models,members,points,theta,weights,link_config,calibration_config,procedure,target,maximum_attempts,seed,cut_theta,progress=None):
    if procedure not in ["candidate_bootstrap","joint_covariance_normal"]:raise ValueError("Unknown uncertainty procedure")
    if target<2 or maximum_attempts<target:raise ValueError("Unsupported uncertainty stopping settings")
    indices={form:[models[form]["items"].index(item) for item in members] for form in ["A","B"]}
    draws={"A":[],"B":[]};resamples={"A":[],"B":[]};rows=[];calibration_rows=[];optimizer_rows=[];valid_count=0
    started=time.perf_counter()
    for attempt in range(maximum_attempts):
        attempt_started=time.perf_counter();vectors={};failures=[]
        for code,form in enumerate(["A","B"]):
            model=models[form];k=len(model["items"]);vector=np.full(2*k,np.nan)
            if procedure=="candidate_bootstrap":
                rng=np.random.default_rng(np.random.SeedSequence([seed,attempt,code]))
                sample=rng.integers(0,len(model["x"]),len(model["x"]))
                resamples[form].append(sample)
                try:
                    fit=fit_mml(model["x"][sample],calibration_config,81,warm_start=np.r_[model["d"],model["a"]],calculate_uncertainty=False)
                    diagnostic=fit["diagnostics"];vector=fit["parameters"]
                    calibration_rows.append({"attempt":attempt,"form":form,**diagnostic,
                                             "all_optimizer_starts":fit["start_diagnostics"].to_json(orient="records"),
                                             "replicate_information":"not computed; point-estimator bootstrap does not need replicate Hessian SEs"})
                    if not diagnostic["numerically_admissible"]:failures.append(form+": unsupported recalibration")
                except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
                    failures.append(form+": "+str(exc));calibration_rows.append({"attempt":attempt,"form":form,"converged":False,"warnings":str(exc)})
            else:
                try:vector=joint_parameter_draw(model,seed,attempt,code)
                except (ValueError,np.linalg.LinAlgError) as exc:failures.append(form+": "+str(exc))
            if not np.isfinite(vector).all() or np.sum(vector[k:])<=0:
                failures.append(form+": nonfinite parameters or inconsistent global orientation")
            # No clipping of any covariance draw, signed slope or intercept.
            draws[form].append(vector);vectors[form]=vector
        fits={}
        if not failures:
            da=vectors["A"][:len(models["A"]["items"])][indices["A"]]
            aa=vectors["A"][len(models["A"]["items"]):][indices["A"]]
            db=vectors["B"][:len(models["B"]["items"])][indices["B"]]
            ab=vectors["B"][len(models["B"]["items"]):][indices["B"]]
            for method,point in points.items():
                try:
                    fit=fit_link(da,aa,db,ab,theta,weights,link_config,method,
                                 starts=[(point["scale_A"],point["shift_B"]),(1.,0.)])
                    fits[method]=fit
                    log=fit["attempts"].copy();log.insert(0,"method",method);log.insert(0,"attempt",attempt);optimizer_rows.append(log)
                    if not fit["supported"]:failures.append(method+": unsupported linking optimization")
                except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as exc:failures.append(method+": "+str(exc))
        valid=not failures and len(fits)==len(points)
        if valid:valid_count+=1
        for method in points:
            fit=fits.get(method,{})
            row={"procedure":procedure,"attempt":attempt,"seed":seed,"method":method,"valid_paired_draw":valid,
                 "valid_draw_number":valid_count if valid else np.nan,"failure":"; ".join(failures),
                 "scale_A":fit.get("scale_A",np.nan),"shift_B":fit.get("shift_B",np.nan),
                 "link_supported":fit.get("supported",False),"boundary_hit":fit.get("boundary_hit",False),
                 "maximum_absolute_gradient":fit.get("maximum_absolute_gradient",np.nan),
                 "normalized_objective":fit.get("normalized_objective",np.nan),
                 "n_supported_starts":fit.get("n_supported_starts",0),
                 "negative_slopes_A":int((vectors["A"][len(models["A"]["items"]):]<0).sum()),
                 "negative_slopes_B":int((vectors["B"][len(models["B"]["items"]):]<0).sum()),
                 "elapsed_seconds_attempt":time.perf_counter()-attempt_started}
            if valid:
                ka,kb=len(models["A"]["items"]),len(models["B"]["items"])
                dt,at=transform_parameters(vectors["B"][:kb],vectors["B"][kb:],fit["scale_A"],fit["shift_B"])
                row["expected_score_difference_BminusA_at_fixed_A_expected42_location"]=float(
                    probabilities(cut_theta,dt,at).sum()-probabilities(cut_theta,vectors["A"][:ka],vectors["A"][ka:]).sum())
            rows.append(row)
        if progress is not None and ((attempt+1)%20==0 or valid_count==target):
            progress(f"{procedure}: {valid_count} valid paired draws / {attempt+1} attempts; {time.perf_counter()-started:.1f}s")
        if valid_count==target:break
    arrays={"parameters_A":np.asarray(draws["A"]),"parameters_B":np.asarray(draws["B"]),
            "anchor_ids":np.asarray(members),"theta_A":theta,"weights":weights,
            "source_parameters_A":np.r_[models["A"]["d"],models["A"]["a"]],
            "source_parameters_B":np.r_[models["B"]["d"],models["B"]["a"]],
            "source_covariance_A":models["A"]["covariance"],"source_covariance_B":models["B"]["covariance"]}
    if procedure=="candidate_bootstrap":
        arrays.update({"sample_indices_A":np.asarray(resamples["A"]),"sample_indices_B":np.asarray(resamples["B"])})
    return {"draws":pd.DataFrame(rows),"calibrations":pd.DataFrame(calibration_rows),
            "optimizers":pd.concat(optimizer_rows,ignore_index=True) if optimizer_rows else pd.DataFrame(),"arrays":arrays}


def percentile_mc_bounds(values,quantile):
    """Order-statistic Monte Carlo bounds for a reference-distribution quantile."""
    values=np.sort(np.asarray(values,float));n=len(values)
    lower=max(0,int(binom.ppf(.025,n,quantile))-1)
    upper=min(n-1,int(binom.ppf(.975,n,quantile)))
    return float(values[lower]),float(values[upper])


def summarize_experiment(result,target,maximum_failure_fraction=.10):
    rows=[]
    for method,group in result["draws"].groupby("method",sort=False):
        valid=group.loc[group.valid_paired_draw];n=len(valid);attempted=len(group)
        supported=n>=target and (attempted-n)/attempted<=maximum_failure_fraction
        row={"procedure":group.procedure.iloc[0],"method":method,"target_valid_draws":target,"valid_draws":n,
             "attempted_draws":attempted,"failed_draws":attempted-n,"failure_fraction":(attempted-n)/attempted,
             "uncertainty_supported":supported,"interval_method":"empirical percentile, linear interpolation; fixed S18 selection",
             "source":"candidate resampling and full MML parameter recalibration" if group.procedure.iloc[0]=="candidate_bootstrap" else "joint asymptotic within-form covariance draws; independent forms"}
        for constant in ["scale_A","shift_B","expected_score_difference_BminusA_at_fixed_A_expected42_location"]:
            values=valid[constant].dropna().to_numpy() if constant in valid else np.array([])
            if len(values)>=2:
                row.update({constant+"_mean":float(values.mean()),constant+"_se":float(values.std(ddof=1)),
                            constant+"_mean_monte_carlo_se":float(values.std(ddof=1)/np.sqrt(len(values))),
                            constant+"_95_low":float(np.quantile(values,.025,method="linear")) if supported else np.nan,
                            constant+"_95_high":float(np.quantile(values,.975,method="linear")) if supported else np.nan})
                for q,label in [(.025,"lower_endpoint"),(.975,"upper_endpoint")]:
                    low,high=percentile_mc_bounds(values,q)
                    row[constant+"_"+label+"_mc95_low"]=low;row[constant+"_"+label+"_mc95_high"]=high
        row["scale_shift_correlation"]=float(valid.scale_A.corr(valid.shift_B)) if n>1 else np.nan
        row["approximate_relative_mc_se_of_sd_normal_reference"]=float(1/np.sqrt(2*(n-1))) if n>1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def cache_key(root,models,members,points,theta,weights,link_config,calibration_config,procedure,target,attempts,seed,cut_theta):
    root=Path(root);h=hashlib.sha256()
    names=["src/linking_transformations.py","src/linking_stocking_lord.py","src/linking_haebara.py","src/linking_uncertainty.py",
           "src/irt_estimation.py","src/irt_2pl.py"]
    for name in names:h.update((root/name).read_bytes())
    for form in ["A","B"]:
        for key in ["x","d","a","covariance"]:h.update(np.ascontiguousarray(models[form][key]).tobytes())
        h.update(json.dumps(models[form]["items"]).encode())
    h.update(np.ascontiguousarray(theta).tobytes());h.update(np.ascontiguousarray(weights).tobytes())
    fields={"anchors":members,"points":{m:[p["scale_A"],p["shift_B"]] for m,p in points.items()},
            "linking":asdict(link_config),"calibration":asdict(calibration_config),"procedure":procedure,
            "target":target,"attempts":attempts,"seed":seed,"cut_theta":cut_theta}
    h.update(json.dumps(fields,sort_keys=True).encode())
    return h.hexdigest()


def cached_experiment(root,out,models,members,points,theta,weights,link_config,calibration_config,procedure,target,maximum_attempts,seed,cut_theta,progress=None):
    dest=Path(out)/"uncertainty_draws";dest.mkdir(exist_ok=True)
    prefix="bootstrap" if procedure=="candidate_bootstrap" else "covariance"
    key=cache_key(root,models,members,points,theta,weights,link_config,calibration_config,procedure,target,maximum_attempts,seed,cut_theta)
    paths={"draws":dest/f"{prefix}_linking_draws.csv","calibrations":dest/f"{prefix}_calibration_diagnostics.csv",
           "optimizers":dest/f"{prefix}_optimizer_starts.csv","arrays":dest/f"{prefix}_parameters.npz"}
    manifest=dest/f"{prefix}_cache_manifest.json"
    if manifest.is_file():
        old=json.loads(manifest.read_text())
        if old.get("fingerprint")==key and all(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==old["sha256"].get(name) for name,path in paths.items()):
            if progress:progress(procedure+": reusing verified exact-input/source cache")
            result={name:pd.read_csv(path,float_precision="round_trip") for name,path in paths.items() if name not in ["arrays","calibrations"]}
            result["calibrations"]=pd.read_csv(paths["calibrations"],float_precision="round_trip") if paths["calibrations"].read_bytes().strip() else pd.DataFrame()
            with np.load(paths["arrays"]) as arrays:result["arrays"]={name:arrays[name] for name in arrays.files}
            return result
    result=experiment(models,members,points,theta,weights,link_config,calibration_config,procedure,target,maximum_attempts,seed,cut_theta,progress)
    for name in ["draws","calibrations","optimizers"]:result[name].to_csv(paths[name],index=False)
    np.savez_compressed(paths["arrays"],**result["arrays"])
    manifest.write_text(json.dumps({"fingerprint":key,"sha256":{name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()}},indent=2)+"\n")
    return result
