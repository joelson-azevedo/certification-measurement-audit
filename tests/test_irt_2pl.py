"""Independent numerical references for logistic signed-slope IRT calculations."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import expit

from src.irt_2pl import (approximate_cut_location,conditional_log_likelihood,eap,expected_score,interpretable_difficulty,
                        item_information,marginal_log_likelihood,normal_quadrature,orient_parameters,response_probability,test_information as summed_test_information)


def test_logistic_probability_known_values():
    result=response_probability([-1,0,1],[0],[2])[:,0]
    np.testing.assert_allclose(result,[1/(1+np.exp(2)),.5,1/(1+np.exp(-2))])


def test_signed_slope_decreases_probability_without_recoding():
    p=response_probability([-2,0,2],[0],[-1])[:,0]
    assert p[0]>.5>p[2] and p[1]==.5


def test_extreme_theta_and_log_probabilities_are_finite():
    probabilities=response_probability([-1e6,0,1e6],[1000,-1000],[6,-6])
    assert np.isfinite(probabilities).all() and np.all((probabilities>=0)&(probabilities<=1))
    values=conditional_log_likelihood([[1,0],[0,1]],[-1000,1000],[10,-10],[6,-6])
    assert np.isfinite(values).all()


def test_conditional_likelihood_matches_manual_products():
    theta=np.array([-1.,.5]);d=np.array([.3,-.2]);a=np.array([1.,-.7])
    x=np.array([[1,0],[0,1]],float)
    expected=[]
    for row in x:
        expected.append([sum(np.log(expit(d[j]+a[j]*t)) if row[j]==1 else np.log(1-expit(d[j]+a[j]*t)) for j in range(2)) for t in theta])
    np.testing.assert_allclose(conditional_log_likelihood(x,theta,d,a),expected,atol=1e-12)


def test_logsumexp_marginal_integration_against_explicit_probability_sum():
    x=np.array([[1,0,np.nan],[0,1,1]],float);d=np.array([.2,-.4,.8]);a=np.array([1.2,-.8,.5])
    raw_nodes,raw_weights=np.polynomial.hermite.hermgauss(31)
    theta=np.sqrt(2)*raw_nodes;weights=raw_weights/np.sqrt(np.pi)
    reference=0.
    for row in x:
        terms=[]
        for t,w in zip(theta,weights):
            probability=1.
            for j,value in enumerate(row):
                if np.isnan(value):continue
                p=expit(d[j]+a[j]*t);probability*=p if value else 1-p
            terms.append(w*probability)
        reference+=np.log(sum(terms))
    assert marginal_log_likelihood(x,d,a,31)==pytest.approx(reference,abs=1e-12)


@pytest.mark.parametrize('nodes',[31,41,81,101])
def test_normal_quadrature_normalization_and_moments(nodes):
    theta,weights=normal_quadrature(nodes)
    assert weights.sum()==pytest.approx(1,abs=1e-14)
    assert weights@theta==pytest.approx(0,abs=1e-13)
    assert weights@(theta**2)==pytest.approx(1,abs=1e-12)


def test_missing_item_likelihood_is_excluded_not_incorrect():
    x=np.array([[1,np.nan,0],[0,np.nan,1]],float);d=np.array([.2,2.,-.4]);a=np.array([1.,2.,-.8])
    actual=marginal_log_likelihood(x,d,a)
    expected=marginal_log_likelihood(x[:,[0,2]],d[[0,2]],a[[0,2]])
    assert actual==pytest.approx(expected,abs=1e-12)
    assert actual!=pytest.approx(marginal_log_likelihood(np.nan_to_num(x,nan=0),d,a))


def test_global_orientation_reflection_keeps_likelihood():
    x=np.array([[1,1,0],[0,1,1]],float);d=np.array([.2,-.4,.1]);a=np.array([-1.,-1.2,.3])
    oriented_d,oriented_a,flipped=orient_parameters(d,a)
    assert flipped and oriented_a.sum()>0
    assert marginal_log_likelihood(x,d,a)==pytest.approx(marginal_log_likelihood(x,oriented_d,oriented_a),abs=1e-12)
    np.testing.assert_equal(x,[[1,1,0],[0,1,1]])


def test_weak_or_uncertain_slope_has_undefined_difficulty():
    b,status=interpretable_difficulty([3.,1.,.2,-.4],[.01,.15,0.,-.8],.10,[.02,.10,.05,.1])
    assert np.isnan(b[:3]).all()
    assert b[3]==pytest.approx(-.5)
    assert 'threshold' in status[0] and 'includes zero' in status[1]


def test_single_item_eap_and_posterior_sd_against_independent_integration():
    density=lambda t:np.exp(-t*t/2)/np.sqrt(2*np.pi)
    denominator=quad(lambda t:expit(t)*density(t),-np.inf,np.inf)[0]
    mean=quad(lambda t:t*expit(t)*density(t),-np.inf,np.inf)[0]/denominator
    second=quad(lambda t:t*t*expit(t)*density(t),-np.inf,np.inf)[0]/denominator
    actual,sd,_=eap([[1]],[0],[1],81)
    assert actual[0]==pytest.approx(mean,abs=1e-10)
    assert sd[0]==pytest.approx(np.sqrt(second-mean*mean),abs=1e-10)


def test_all_missing_eap_is_prior_not_invented_evidence():
    mean,sd,_=eap([[np.nan,np.nan,np.nan]],[1,-1,.2],[1,-1,.4],81)
    assert mean[0]==pytest.approx(0,abs=1e-12) and sd[0]==pytest.approx(1,abs=1e-12)


def test_negative_item_changes_posterior_in_expected_direction():
    mean,_,_=eap([[1,1,1],[1,1,0]],[0,0,0],[1,1,-1],81)
    assert mean[1]>mean[0]


def test_item_information_and_sum_use_signed_slope_squared():
    actual=item_information([0],[0,0],[2,-2])[0]
    np.testing.assert_allclose(actual,[1,1])
    assert summed_test_information(0,[0,0],[2,-2])==pytest.approx(2)
    h=1e-5;p=response_probability([.3],[.2],[-1.3])[0,0]
    derivative=(response_probability([.3+h],[.2],[-1.3])[0,0]-response_probability([.3-h],[.2],[-1.3])[0,0])/(2*h)
    assert item_information(.3,[.2],[-1.3])[0]==pytest.approx(derivative**2/(p*(1-p)),rel=1e-8)


def test_expected_score_is_sum_not_a_passing_probability():
    assert expected_score(0,[0,0,0],[1,2,-1])==pytest.approx(1.5)


def test_expected_score_cut_mapping_known_reference():
    result=approximate_cut_location(np.zeros(60),np.ones(60),42)
    assert result['theta_at_expected_score']==pytest.approx(np.log(.7/.3),abs=1e-10)
    assert result['information_at_location']==pytest.approx(60*.7*.3,abs=1e-10)
    assert result['conditional_sem_at_location']==pytest.approx(1/np.sqrt(12.6),abs=1e-10)
    assert result['n_roots_in_prespecified_range']==1 and result['tcc_monotone_on_grid']


def test_unreachable_cut_is_not_fabricated():
    result=approximate_cut_location([0,0,0],[1,1,1],42)
    assert result['n_roots_in_prespecified_range']==0 and np.isnan(result['theta_at_expected_score'])
