"""Aggregate counterfactual bookkeeping; never replaces any recorded decision."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .irt_observed_score_equating import mid_cdf,mid_rank_at


def counterfactual_by_score(source_pmf,reference_pmf,conventions):
    scores=np.arange(len(source_pmf));percentiles=mid_cdf(source_pmf);rows=[]
    for name,specification in conventions.items():
        boundary=specification["reference_score_boundary"];threshold=float(mid_rank_at(boundary,reference_pmf))
        if specification["comparison"]==">=":passed=percentiles>=threshold
        elif specification["comparison"]==">":passed=percentiles>threshold
        else:raise ValueError("Unknown research classification comparison")
        rows.extend({"convention":name,"original_raw_score":int(score),"reference_score_boundary":boundary,
                     "reference_mid_percentile_boundary":threshold,"source_mid_percentile":percentiles[score],
                     "comparison":specification["comparison"],"counterfactual_pass":bool(passed[score]),
                     "interpretation":"counterfactual threshold ordering; no score clipping or official decision replacement"} for score in scores)
    return pd.DataFrame(rows)


def aggregate_impact(raw_scores,recorded_pass,source_pmf,reference_pmf,conventions,cut=42):
    raw=np.asarray(raw_scores,int);recorded=np.asarray(recorded_pass,bool)
    if raw.shape!=recorded.shape or (raw<0).any() or (raw>=len(source_pmf)).any():raise ValueError("Valid integer first-attempt scores required")
    original=raw>=cut
    if not np.array_equal(original,recorded):raise ValueError("Recorded decisions disagree with unchanged operational rule")
    counts=np.bincount(raw,minlength=len(source_pmf));rules=counterfactual_by_score(source_pmf,reference_pmf,conventions)
    overall=[];strata=[]
    for convention,group in rules.groupby("convention",sort=False):
        predicted=group.counterfactual_pass.to_numpy(bool);up=(np.arange(len(counts))<cut)&predicted;down=(np.arange(len(counts))>=cut)&~predicted
        total=int(counts.sum());newpass=int(counts@predicted)
        overall.append({"convention":convention,"n_first_attempts":total,"original_pass":int(recorded.sum()),"original_fail":int((~recorded).sum()),
                        "counterfactual_pass":newpass,"counterfactual_fail":total-newpass,"upward_changes":int(counts@up),"downward_changes":int(counts@down),
                        "unchanged_pass":int(counts@((np.arange(len(counts))>=cut)&predicted)),
                        "unchanged_fail":int(counts@((np.arange(len(counts))<cut)&~predicted)),
                        "minimum_B_integer_passing_research_rule":int(np.flatnonzero(predicted)[0]) if predicted.any() else np.nan})
        for row in group.itertuples():
            s=row.original_raw_score;n=int(counts[s]);old=s>=cut;new=row.counterfactual_pass
            strata.append({"convention":convention,"original_raw_score":s,"n":n,"original_pass":n if old else 0,"original_fail":n if not old else 0,
                           "counterfactual_pass":n if new else 0,"counterfactual_fail":n if not new else 0,
                           "upward_changes":n if new and not old else 0,"downward_changes":n if old and not new else 0,
                           "near_operational_cut_38_to_46":38<=s<=46})
    return pd.DataFrame(overall),pd.DataFrame(strata),rules
