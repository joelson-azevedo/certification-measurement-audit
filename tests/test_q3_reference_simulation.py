"""Public, controlled independent/dependent fixtures for the Q3 reference."""
import numpy as np
import pandas as pd
import pytest

from src.local_dependence import pairwise_q3, residuals_at_eap
from src.q3_reference_simulation import (ReferenceConfig, classify_pairs, generate_responses,
                                          reference_thresholds, simulate_reference)


@pytest.fixture(scope='module')
def controlled_reference():
    d = np.array([-.8,-.3,0,.4,.8,-.5,.1,.6])
    a = np.array([.8,1.1,.9,1.3,.7,1.2,-.8,.2])
    mask = np.ones((600,8),bool); mask[::17,2]=False; mask[::31,5]=False
    config = ReferenceConfig(60,65,41)
    values, log = simulate_reference(d,a,mask,7841,config)
    return d,a,mask,config,values,log


def test_deterministic_response_generation_and_exact_mask():
    mask = np.ones((100,3),bool);mask[::3,1]=False
    first = generate_responses(np.array([.1,0,-.2]),np.array([1,-1,.2]),mask,88,0)
    same = generate_responses(np.array([.1,0,-.2]),np.array([1,-1,.2]),mask,88,0)
    next_draw = generate_responses(np.array([.1,0,-.2]),np.array([1,-1,.2]),mask,88,1)
    np.testing.assert_array_equal(first,same)
    assert np.array_equal(np.isfinite(first),mask)
    assert not np.array_equal(np.nan_to_num(first),np.nan_to_num(next_draw))


def test_null_reference_reestimates_eap_and_matches_independent_first_replicate(controlled_reference):
    d,a,mask,config,values,log=controlled_reference
    x=generate_responses(d,a,mask,7841,0)
    residual,mean,_=residuals_at_eap(x,d,a,41)
    pairs,_,diag=pairwise_q3(residual)
    np.testing.assert_allclose(values[0],pairs.Q3_adjusted,atol=1e-15)
    assert log.eap_reestimated.all() and log.valid_replicate.all()
    assert log.mean_simulated_eap.iloc[0] == pytest.approx(mean.mean())
    assert log.mean_simulated_eap.nunique() == 60
    again,log2=simulate_reference(d,a,mask,7841,ReferenceConfig(3,3,41))
    np.testing.assert_array_equal(again,values[:3])
    assert log2.minimum_pairwise_n.eq(int((mask[:,2]&mask[:,5]).sum())).all()


def test_familywise_quantiles_are_maximum_based_higher_order_statistics():
    values=np.array([[.01,-.01,.03],[.04,.05,-.03],[-.02,.02,.08],[.06,.02,.01]])
    observed=np.array([.1,.02,.01])
    thresholds,point=reference_thresholds(values,observed)
    assert thresholds['familywise_95_threshold']==.08
    assert thresholds['familywise_99_threshold']==.08
    assert thresholds['maximum_monte_carlo_tail_area_plus_one']==.2
    assert thresholds['monte_carlo_tail_area_resolution']==.2
    assert thresholds['reference_tail_probability_mc_wilson_95_low']==pytest.approx(0,abs=1e-15)
    assert thresholds['reference_tail_probability_mc_wilson_95_high']==pytest.approx(.4898908364545973)
    assert point.pointwise_mc_tail_area_plus_one.tolist()==[.2, .8, .8]
    assert thresholds['familywise_95_threshold'] > np.quantile(values,.95)


def test_classification_distinguishes_familywise_exploratory_negative_and_ties():
    observed=pd.DataFrame({'Q3_adjusted':[.31,.21,.11,.10,-.3]})
    point=pd.DataFrame({'pointwise_adjusted_p95':[.1]*5})
    classified=classify_pairs(observed,{'familywise_95_threshold':.2,'familywise_99_threshold':.3},point)
    assert classified.evidence_classification.tolist()==[
        'exceeds simulated family-wise 99% threshold','exceeds simulated family-wise 95% threshold',
        'elevated in exploratory pointwise analysis only','no notable positive residual association',
        'no notable positive residual association']
    assert classified.exceeds_familywise_95.sum()==2 and classified.exceeds_familywise_99.sum()==1


def test_independent_null_bias_and_dependent_pair_sensitivity(controlled_reference):
    d,a,mask,config,values,log=controlled_reference
    # Independent data's finite-test residual bias is reproduced, not set to zero.
    assert log.mean_off_diagonal_Q3.mean() < -.01
    np.testing.assert_allclose(values.mean(axis=1),0,atol=1e-15)
    # Many independent draws have approximately their nominal maximum exceedance rate.
    cutoff=np.quantile(np.maximum(0,values.max(axis=1)),.95,method='higher')
    assert ((values.max(axis=1)>cutoff).sum()) <= 3
    # Independent held-out null draws use another seed, not the reference quantile itself.
    held_out_maxima=[]
    for attempt in range(20):
        independent=generate_responses(d,a,mask,8110,attempt)
        r,_,_=residuals_at_eap(independent,d,a,41)
        held_out_maxima.append(pairwise_q3(r)[0].Q3_adjusted.max())
    assert (np.asarray(held_out_maxima)>cutoff).sum() <= 5
    x=generate_responses(d,a,mask,45002,0)
    x[:,1]=x[:,0]  # Explicit controlled dependence, unrelated to audit data.
    residual,_,_=residuals_at_eap(x,d,a,41)
    pairs,_,_=pairwise_q3(residual)
    thresholds,point=reference_thresholds(values,pairs.Q3_adjusted.to_numpy(),config)
    flagged=classify_pairs(pairs,thresholds,point)
    target=flagged.query('item_index_1 == 0 and item_index_2 == 1').iloc[0]
    assert target.Q3_adjusted > .9 and target.exceeds_familywise_99


def test_failure_is_retained_not_fabricated():
    mask=np.zeros((20,3),bool)
    values,log=simulate_reference(np.zeros(3),np.ones(3),mask,78,ReferenceConfig(2,3,31))
    assert values.shape==(0,3) and len(log)==3
    assert not log.valid_replicate.any() and log.failure.str.contains('unavailable').all()


@pytest.mark.parametrize('args', [(1,2,31),(3,2,31),(2,2,1)])
def test_unsupported_reference_config_rejected(args):
    with pytest.raises(ValueError): ReferenceConfig(*args)


def test_nonfinite_or_incompatible_reference_rejected():
    with pytest.raises(ValueError): reference_thresholds(np.ones((3,2)),np.ones(3))
    with pytest.raises(ValueError): reference_thresholds(np.array([[1,np.nan],[1,2]]),np.ones(2))
