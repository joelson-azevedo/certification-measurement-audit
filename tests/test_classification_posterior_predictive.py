"""Exact likelihood/posterior arithmetic and prediction integration, not EAP substitution."""
import numpy as np
import pandas as pd
import pytest
from scipy.special import expit
from src.classification_posterior_predictive import response_posterior,predictive_moments,predict_candidates,predictive_score_groups
from src.score_distribution import normal_quadrature


def test_missing_response_is_excluded_from_posterior_likelihood():
    theta=np.array([-1.,1.]);weights=np.array([.5,.5]);d=np.array([0.,1.]);a=np.array([1.,-2.])
    result=response_posterior([[1.,np.nan]],d,a,theta,weights)
    expected=weights*expit(theta);expected/=expected.sum()
    np.testing.assert_allclose(result["posterior"][0],expected,atol=1e-15)
    allmissing=response_posterior([[np.nan,np.nan]],d,a,theta,weights)
    np.testing.assert_allclose(allmissing["posterior"][0],weights,atol=1e-15)
    not_missing=response_posterior([[1.,0.]],d,a,theta,weights)
    assert not np.allclose(not_missing["posterior"],result["posterior"])


def test_posterior_prediction_differs_from_EAP_plugin_and_variance_decomposes():
    theta,w=normal_quadrature(161);result=response_posterior([[np.nan]],[1.],[2.],theta,w)
    q=expit(1+2*theta);moments=predictive_moments(result["posterior"],q)
    predicted=moments["future_pass_probability"][0];plugin=expit(1+2*result["mean"][0])
    assert abs(predicted-plugin)>.05
    assert moments["variance_decomposition_error"][0]<1e-14
    assert predicted==pytest.approx(w@q,abs=1e-14)


def test_B_source_prior_is_mapped_to_A_instead_of_replaced_by_A_prior():
    # Identical 60-item forms simplify an independent binomial-tail comparison.
    d=np.full(60,.3);a=np.ones(60);x=np.zeros((2,60));x[0,:42]=1;x[1,:38]=1
    scale=1.2;shift=-.3;dB=d.copy();aB=a.copy();linked_d=dB-aB*shift/scale;linked_a=aB/scale
    result=predict_candidates(x,"B",dB,aB,d,a,linked_d,linked_a,scale,shift,[42],[81,161,321,641])
    np.testing.assert_allclose(result["theta_A"],scale*result["theta_source"]+shift)
    from scipy.stats import binom
    own=binom.sf(41,60,expit(.3+result["theta_source"]))
    np.testing.assert_allclose(result["moments"]["B42"]["future_pass_probability"],result["posterior"]["posterior"]@own,atol=1e-13)


def test_raw_score_grouping_uses_individual_predictions_and_recorded_rule():
    frame=pd.DataFrame({"source_form":["A","A","A"],"original_raw_score":[42,42,41],"target_form":"A","future_cut":42,"target_role":"same_form_original_rule",
                        "future_pass_probability":[.4,.6,.45],"same_as_original_recorded_decision_probability":[.4,.6,.55],"original_recorded_decision":["PASS","PASS","FAIL"],
                        "posterior_sd_conditional_pass_probability":.2,"within_theta_future_Bernoulli_variance":.15,"between_theta_pass_probability_variance":.03,
                        "posterior_prediction_minus_EAP_plugin":[-.1,.05,.03]})
    original=frame.copy(deep=True);groups=predictive_score_groups(frame)
    g=groups.loc[groups.original_raw_score.eq(42)].iloc[0]
    assert g.n_first_attempts==2 and g.mean_future_pass_probability==pytest.approx(.5)
    assert g.future_pass_p10==pytest.approx(.42) and g.same_original_decision_p90==pytest.approx(.58)
    assert groups.loc[groups.original_raw_score.eq(41),"mean_same_original_decision_probability"].iloc[0]==pytest.approx(.55)
    pd.testing.assert_frame_equal(frame,original)


def test_invalid_response_and_posterior_resolution_failure_are_explicit():
    with pytest.raises(ValueError):response_posterior([[2.]],[0.],[1.],[0.,1.],[.5,.5])
    with pytest.raises(ValueError,match="quadrature"):
        predict_candidates(np.ones((1,60)),"A",np.zeros(60),np.ones(60),np.zeros(60),np.ones(60),np.zeros(60),np.ones(60),1,0,[42],[5,7],1e-20,1e-20)
