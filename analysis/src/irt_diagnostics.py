"""Conditional measurement summaries and explicitly descriptive model checks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit,logsumexp
from scipy.stats import spearmanr

from .distractors import performance_groups
from .irt_2pl import approximate_cut_location,conditional_log_likelihood,eap,expected_score,interpretable_difficulty,item_information,marginal_components,normal_quadrature,test_information,validate_responses


def parameter_tables(fit,members,metadata,minimum_abs_slope=.10):
    d,a=fit['intercept'],fit['slope'];k=len(a)
    covariance=fit['covariance']
    se=np.sqrt(np.maximum(np.diag(covariance),0))
    se_d,se_a=se[:k],se[k:]
    b,status=interpretable_difficulty(d,a,minimum_abs_slope,se_a)
    se_b=np.full(k,np.nan)
    for j in range(k):
        if np.isfinite(b[j]) and np.isfinite(covariance[j,j]):
            derivative=np.array([-1/a[j],d[j]/(a[j]**2)])
            sub=covariance[np.ix_([j,k+j],[j,k+j])]
            variance=derivative@sub@derivative
            if variance>=0:se_b[j]=np.sqrt(variance)
    meta=metadata.set_index('item_id').loc[members]
    grid=np.linspace(-4,4,321)
    info=item_information(grid,d,a)
    parameters=pd.DataFrame({'item_id':members,'domain':meta.domain.to_numpy(),'anchor_item':meta.anchor_item.to_numpy(),
                             'intercept_d':d,'discrimination_a':a,'difficulty_b':b,'difficulty_status':status,
                             'negative_discrimination':a<0,'near_zero_discrimination':np.abs(a)<minimum_abs_slope,
                             'intercept_boundary_hit':fit['boundary_mask'][:k],'slope_boundary_hit':fit['boundary_mask'][k:],
                             'information_theta_zero':item_information(0,d,a),
                             'maximum_item_information_in_grid':info.max(axis=0),
                             'theta_at_maximum_item_information_in_grid':grid[np.argmax(info,axis=0)],
                             'model_numerically_admissible':fit['diagnostics']['numerically_admissible']})
    uncertainty=pd.DataFrame({'item_id':members,'se_intercept_d':se_d,'se_discrimination_a':se_a,'se_difficulty_b_delta':se_b,
                              'intercept_95_low':d-1.959963984540054*se_d,'intercept_95_high':d+1.959963984540054*se_d,
                              'slope_95_low':a-1.959963984540054*se_a,'slope_95_high':a+1.959963984540054*se_a,
                              'difficulty_95_low_delta':b-1.959963984540054*se_b,'difficulty_95_high_delta':b+1.959963984540054*se_b,
                              'difficulty_status':status,'uncertainty_status':fit['diagnostics']['uncertainty_status']})
    return parameters,uncertainty


def candidate_eap(x,fit,sittings):
    mean,sd,_=eap(x,fit['intercept'],fit['slope'],fit['diagnostics']['quadrature_nodes'])
    result=sittings[['administration_id','candidate_id','form']].reset_index(drop=True).copy()
    result['n_observed_items']=np.isfinite(x).sum(axis=1)
    result['n_omitted_items']=np.isnan(x).sum(axis=1)
    result['eap_theta']=mean;result['posterior_sd']=sd
    result['interpretation']='provisional model EAP; not official certification score'
    return result


def observed_expected_checks(x,fit,members):
    x=validate_responses(x);d,a=fit['intercept'],fit['slope']
    theta,prior,_,posterior=marginal_components(x,d,a,fit['diagnostics']['quadrature_nodes'])
    probabilities=expit(theta[:,None]*a+d)
    joint=conditional_log_likelihood(x,theta,d,a)+np.log(prior)[None,:]
    observed=np.isfinite(x)
    eta=theta[:,None]*a+d
    rows,bands=[],[]
    for j,item in enumerate(members):
        focal=np.where(observed[:,j,None],np.nan_to_num(x[:,j],nan=0)[:,None]*(-np.logaddexp(0,-eta[:,j]))[None,:]
                       +(1-np.nan_to_num(x[:,j],nan=0))[:,None]*(-np.logaddexp(0,eta[:,j]))[None,:],0)
        log_rest=joint-focal
        log_rest-=logsumexp(log_rest,axis=1)[:,None]
        rest_posterior=np.exp(log_rest)
        rest_mean=rest_posterior@theta
        predicted_rest=rest_posterior@probabilities[:,j]
        # Mathematically identical rest patterns can differ at machine precision
        # after subtracting the focal likelihood. Preserve those numerical ties.
        groups,definition=performance_groups(np.round(rest_mean,12))
        valid=observed[:,j];n=int(valid.sum())
        group_differences,group_weights=[],[]
        for group in range(1,definition['n_groups']+1):
            mask=(groups==group)&valid
            count=int(mask.sum())
            observed_p=float(x[mask,j].mean()) if count else np.nan
            expected_p=float(predicted_rest[mask].mean()) if count else np.nan
            bands.append({'item_id':item,'rest_eap_group':group,'n_observed':count,
                          'observed_proportion_correct':observed_p,'model_expected_proportion_correct_given_rest':expected_p,
                          'observed_minus_expected':observed_p-expected_p,'mean_rest_eap':float(rest_mean[mask].mean()) if count else np.nan,
                          'grouping':'five approximate groups; posterior means rounded to 12 decimals for numerical ties; fitted parameters reuse same sample'})
            if count:
                group_differences.append(observed_p-expected_p);group_weights.append(count)
        differences=np.asarray(group_differences)
        rows.append({'item_id':item,'n_observed':n,'n_omitted':int((~valid).sum()),'n_correct':int(np.nansum(x[:,j])),
                     'observed_proportion_correct':float(np.nanmean(x[:,j])) if n else np.nan,
                     'prior_marginal_expected_proportion':float(prior@probabilities[:,j]),
                     'in_sample_posterior_fitted_proportion':float((posterior[valid]@probabilities[:,j]).mean()) if n else np.nan,
                     'observed_minus_prior_expected':float(np.nanmean(x[:,j])-prior@probabilities[:,j]) if n else np.nan,
                     'rest_eap_group_weighted_rmse':float(np.sqrt(np.average(differences**2,weights=group_weights))) if len(differences) else np.nan,
                     'rest_eap_group_maximum_absolute_difference':float(np.max(np.abs(differences))) if len(differences) else np.nan,
                     'fit_check_scope':'descriptive in-sample; focal response omitted from grouping posterior, but parameters estimated on same sample'})
    return pd.DataFrame(rows),pd.DataFrame(bands)


def replicated_model_checks(x,fit,members,replicates=200,seed=20261024):
    """Fixed-parameter model replications, not a parameter-posterior predictive test.

    New theta ~ N(0,1) and binary responses are drawn; the observed missing mask
    stays fixed. No refitting, p-values, item-pair diagnostics or classification.
    """
    x=validate_responses(x);n,k=x.shape
    mask=np.isfinite(x);denominator=mask.sum(axis=0)
    rng=np.random.default_rng(seed)
    item_rates=[];score_statistics=[]
    for _ in range(replicates):
        theta=rng.standard_normal(n)
        probability=expit(fit['intercept']+theta[:,None]*fit['slope'])
        replicated=(rng.random((n,k))<probability)*mask
        total=replicated.sum(axis=1)
        item_rates.append(np.divide(replicated.sum(axis=0),denominator,out=np.full(k,np.nan),where=denominator>0))
        score_statistics.append([total.mean(),total.std(ddof=1),*np.quantile(total,[.1,.5,.9])])
    rates=np.asarray(item_rates);stats=np.asarray(score_statistics)
    total=np.nan_to_num(x,nan=0).sum(axis=1)
    observed_stats=[total.mean(),total.std(ddof=1),*np.quantile(total,[.1,.5,.9])]
    names=['mean_observed_correct_total','sd_observed_correct_total','total_p10','total_median','total_p90']
    global_rows=[]
    for j,name in enumerate(names):
        low,high=np.quantile(stats[:,j],[.025,.975])
        global_rows.append({'statistic':name,'observed':float(observed_stats[j]),'replicated_mean':float(stats[:,j].mean()),
                            'replicated_95_low':float(low),'replicated_95_high':float(high),
                            'observed_outside_replicated_95':bool(observed_stats[j]<low or observed_stats[j]>high),
                            'replicates':replicates,'seed':seed,
                            'interpretation':'fixed fitted parameters; latent/model sampling only; descriptive, not independent validation or classification consistency'})
    intervals=np.nanquantile(rates,[.025,.975],axis=0)
    item_rows=pd.DataFrame({'item_id':members,'model_replicated_p_95_low':intervals[0],
                            'model_replicated_p_95_high':intervals[1],'model_replicated_p_mean':np.nanmean(rates,axis=0)})
    return pd.DataFrame(global_rows),item_rows


def information_table(fit,theta_grid):
    information=test_information(theta_grid,fit['intercept'],fit['slope'])
    sem=np.divide(1,np.sqrt(information),out=np.full_like(information,np.nan),where=information>0)
    return pd.DataFrame({'theta':theta_grid,'expected_full_response_correct_total':expected_score(theta_grid,fit['intercept'],fit['slope']),
                         'test_information':information,'conditional_measurement_sem':sem,
                         'scope':'model conditional local independence; own identified form/specification scale'})


def cut_region(fit,full_form=True):
    if not full_form or len(fit['slope'])!=60:
        return {'target_expected_full_form_score':np.nan,'theta_at_expected_score':np.nan,
                'information_at_location':np.nan,'conditional_sem_at_location':np.nan,
                'mapping_status':'not mapped: 59-item sensitivity has no replacement operational threshold',
                'tcc_monotone_on_grid':np.nan,'n_roots_in_prespecified_range':0,'all_roots':[]}
    mapped=approximate_cut_location(fit['intercept'],fit['slope'],42)
    for score in [40,44]:
        mapped[f'theta_at_expected_score_{score}']=approximate_cut_location(fit['intercept'],fit['slope'],score)['theta_at_expected_score']
    mapped['interpretation']='hypothetical fully observed model expected score 42/60; not official theta cut or classification reliability'
    return mapped


def sensitivity_comparison(base,other,base_members,other_members,base_eap,other_eap):
    base_index={item:j for j,item in enumerate(base_members)}
    deltas=[]
    for j,item in enumerate(other_members):
        i=base_index[item]
        deltas.append({'item_id':item,'intercept_base':base['intercept'][i],'intercept_comparison':other['intercept'][j],
                       'slope_base':base['slope'][i],'slope_comparison':other['slope'][j],
                       'intercept_difference':other['intercept'][j]-base['intercept'][i],
                       'slope_difference':other['slope'][j]-base['slope'][i]})
    paired=base_eap[['administration_id','eap_theta']].merge(other_eap[['administration_id','eap_theta']],on='administration_id',suffixes=('_base','_comparison'),validate='one_to_one')
    statistics={'matched_items':len(deltas),'matched_candidates':len(paired),
                'maximum_absolute_intercept_difference':max(abs(row['intercept_difference']) for row in deltas),
                'maximum_absolute_slope_difference':max(abs(row['slope_difference']) for row in deltas),
                'candidate_eap_pearson':float(paired.eap_theta_base.corr(paired.eap_theta_comparison)),
                'candidate_eap_spearman':float(spearmanr(paired.eap_theta_base,paired.eap_theta_comparison).statistic),
                'candidate_eap_mean_absolute_difference':float(np.abs(paired.eap_theta_base-paired.eap_theta_comparison).mean()),
                'candidate_eap_maximum_absolute_difference':float(np.abs(paired.eap_theta_base-paired.eap_theta_comparison).max()),
                'comparison_scope':'within form standardized N(0,1) convention; descriptive sensitivity, not linking; no likelihood/AIC/BIC difference across item sets'}
    return statistics,pd.DataFrame(deltas)
