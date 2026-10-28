"""Rest-score Mantel-Haenszel complement and descriptive cross-form bookkeeping."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .dif import Z95
from .distractors import performance_groups


def mantel_haenszel_tables(tables,minimum_informative_strata=2):
    """North/South correct-response odds, no continuity correction or cell shift.

    Each table is [[North correct, North incorrect],
                   [South correct, South incorrect]].
    Log-OR uncertainty uses the Robins-Breslow-Greenland variance.
    """
    table=np.asarray(tables,float)
    if table.ndim!=3 or table.shape[1:]!=(2,2) or not np.isfinite(table).all() or (table<0).any():
        raise ValueError("Nonnegative K x 2 x 2 count tables required")
    a,b,c,d=table[:,0,0],table[:,0,1],table[:,1,0],table[:,1,1]
    n=a+b+c+d
    informative=(n>1)&((a+b)>0)&((c+d)>0)&((a+c)>0)&((b+d)>0)
    result={"n_strata":len(table),"n_informative_strata":int(informative.sum()),
            "n_strata_with_zero_cells":int((table==0).any(axis=(1,2)).sum()),
            "n_strata_with_observed_cell_below_5":int((table<5).any(axis=(1,2)).sum()),
            "continuity_correction":False,"zero_cell_shift":False,
            "odds_ratio_orientation":"North correct/incorrect odds divided by South correct/incorrect odds",
            "common_odds_ratio":np.nan,"log_odds_ratio_se_rbg":np.nan,"odds_ratio_95_low":np.nan,
            "odds_ratio_95_high":np.nan,"mh_chi_square":np.nan,"mh_p_uniform_only":np.nan}
    if informative.sum()<minimum_informative_strata:
        return {**result,"status":"unavailable: insufficient informative strata"}
    a,b,c,d,n=[v[informative] for v in [a,b,c,d,n]]
    expected_a=(a+b)*(a+c)/n
    variance_a=(a+b)*(c+d)*(a+c)*(b+d)/(n*n*(n-1))
    statistic=(np.sum(a-expected_a)**2)/variance_a.sum()
    r=a*d/n;s=b*c/n;p=(a+d)/n;q=(b+c)/n
    R,S=r.sum(),s.sum()
    result.update({"mh_chi_square":float(statistic),"mh_p_uniform_only":float(chi2.sf(statistic,1)),
                   "n_informative_strata_with_expected_cell_below_5":int(np.any(np.stack([
                       (a+b)*(a+c)/n,(a+b)*(b+d)/n,(c+d)*(a+c)/n,(c+d)*(b+d)/n],axis=1)<5,axis=1).sum())})
    if R<=0 or S<=0:
        return {**result,"status":"test only; common OR/interval unavailable: zero pooled cross-product"}
    variance=(np.sum(p*r)/(2*R*R)+np.sum(p*s+q*r)/(2*R*S)+np.sum(q*s)/(2*S*S))
    if not variance>0:return {**result,"status":"test only; unsupported log-OR variance"}
    se=np.sqrt(variance);log_or=np.log(R/S)
    return {**result,"status":"available; uniform-DIF complement, approximate rest-score strata",
            "common_odds_ratio":float(R/S),"log_odds_ratio_se_rbg":float(se),
            "odds_ratio_95_low":float(np.exp(log_or-Z95*se)),"odds_ratio_95_high":float(np.exp(log_or+Z95*se))}


def mantel_haenszel(y,g,rest_score,requested_strata=10,minimum_informative_strata=2):
    y=np.asarray(y,int);g=np.asarray(g,int);score=np.asarray(rest_score,float)
    if not np.isin(y,[0,1]).all() or not np.isin(g,[0,1]).all() or not np.isfinite(score).all():
        raise ValueError("Valid binary outcome/group and observed rest score required")
    labels,definition=performance_groups(np.round(score,12),requested_strata)
    rows=[];tables=[]
    for label in range(1,definition['n_groups']+1):
        selected=labels==label
        a=int(np.sum(selected&(g==1)&(y==1)));b=int(np.sum(selected&(g==1)&(y==0)))
        c=int(np.sum(selected&(g==0)&(y==1)));d=int(np.sum(selected&(g==0)&(y==0)))
        tables.append([[a,b],[c,d]])
        rows.append({"stratum":label,"rest_score_min":float(score[selected].min()),
                     "rest_score_max":float(score[selected].max()),"n":int(selected.sum()),
                     "north_correct":a,"north_incorrect":b,"south_correct":c,"south_incorrect":d})
    result=mantel_haenszel_tables(tables,minimum_informative_strata)
    result["stratification"]="ten approximate pooled rest-score groups; ties preserved; no group-specific cutpoints"
    return result,pd.DataFrame(rows)


def compare_results(primary,sensitivity):
    keys=['form','item_id']
    columns=keys+['q_omnibus_within_form','screened_within_form','absolute_reference_mean_difference','signed_reference_mean_difference']
    paired=primary[columns].merge(sensitivity[columns],on=keys,suffixes=('_primary','_sensitivity'),validate='one_to_one')
    paired['screen_change']=paired.screened_within_form_primary.ne(paired.screened_within_form_sensitivity)
    paired['absolute_effect_change']=paired.absolute_reference_mean_difference_sensitivity-paired.absolute_reference_mean_difference_primary
    paired['q_change']=paired.q_omnibus_within_form_sensitivity-paired.q_omnibus_within_form_primary
    return paired


def common_item_comparison(primary,common_members):
    common=primary.loc[primary.item_id.isin(common_members)]
    a=common.loc[common.form.eq('A')].set_index('item_id')
    b=common.loc[common.form.eq('B')].set_index('item_id')
    if set(a.index)!=set(common_members) or set(b.index)!=set(common_members):
        raise ValueError("Common items require separate A and B records")
    fields=['screened_within_form','q_omnibus_within_form','p_uniform','p_nonuniform','M2_group_coefficient',
            'M2_group_95_low','M2_group_95_high','M2_interaction_coefficient','M2_interaction_95_low','M2_interaction_95_high',
            'signed_reference_mean_difference','absolute_reference_mean_difference','pattern_interpretation']
    rows=[]
    for item in common_members:
        left,right=a.loc[item],b.loc[item]
        row={'item_id':item,**{field+'_A':left[field] for field in fields},**{field+'_B':right[field] for field in fields}}
        signs=[np.sign(v) if np.isfinite(v) else np.nan for v in [left.signed_reference_mean_difference,right.signed_reference_mean_difference]]
        row['same_integrated_direction']=bool(signs[0]==signs[1]) if np.isfinite(signs).all() else None
        row['screened_in_both_forms']=bool(left.screened_within_form and right.screened_within_form)
        row['scope']='separate conditional evidence; own scales/reference distributions; no pooled estimate, linking or invariance claim'
        rows.append(row)
    return pd.DataFrame(rows)
