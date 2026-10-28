"""Exact conditional tails and limits from controlled signed models."""
from itertools import product
import numpy as np
import pytest
from scipy.special import expit
from src.classification_probabilities import passing_curves,score_tail,passing_review,limiting_pass_probability,probability_band_intervals


@pytest.mark.parametrize("cut",[0,1,2,3,4])
def test_small_form_tail_matches_direct_pattern_enumeration(cut):
    theta=np.array([-1.3,0.,.7]);d=np.array([.3,-.5,.8]);a=np.array([1.2,-.9,0.])
    tails,pmf,p=passing_curves(theta,d,a,[cut]);expected=[]
    for probabilities in p:
        total=0.
        for response in product([0,1],repeat=3):
            x=np.array(response)
            if x.sum()>=cut:total+=np.prod(np.where(x,probabilities,1-probabilities))
        expected.append(total)
    np.testing.assert_allclose(tails[cut],expected,atol=3e-16)
    np.testing.assert_allclose(tails[cut],pmf[:,cut:].sum(axis=1),atol=0)


@pytest.mark.parametrize("cut",[-1,1.5,5])
def test_invalid_integer_boundary_is_rejected(cut):
    with pytest.raises(ValueError):score_tail(np.array([.2,.5,.3]),cut)


def test_link_direction_preserves_conditional_tails_with_signed_slopes():
    d=np.array([.4,-.3,.8]);a=np.array([1.1,-.7,1e-12]);s=1.3;t=-.4;theta=np.linspace(-3,3,31)
    linked_d=d-a*t/s;linked_a=a/s
    np.testing.assert_allclose(passing_curves(theta,linked_d,linked_a,[2])[0][2],passing_curves((theta-t)/s,d,a,[2])[0][2],atol=3e-16)


def test_signed_limiting_behavior_and_detected_nonmonotonicity():
    d=np.zeros(2);a=np.array([1.,-2.]);theta=np.linspace(-6,6,601)
    p=passing_curves(theta,d,a,[1])[0][1];review=passing_review(theta,p,d,a,1)
    assert not review["monotonic_on_checked_grid"] and review["n_decreasing_intervals_beyond_tolerance"]>0
    assert limiting_pass_probability(d,a,1,True)==1 and limiting_pass_probability(d,a,1,False)==1


def test_probability_regions_use_exact_roots_and_allow_multiple_components():
    grid=np.linspace(-4,4,801);p=expit(grid)
    bands=probability_band_intervals(grid,p,.2,.8,expit)
    assert len(bands)==1
    np.testing.assert_allclose([bands.theta_lower.iloc[0],bands.theta_upper.iloc[0]],[np.log(.2/.8),np.log(.8/.2)],atol=1e-11)
    wave=(np.sin(grid)+1)/2;many=probability_band_intervals(grid,wave,.4,.6,lambda t:(np.sin(t)+1)/2)
    assert len(many)>=3


def test_expected_score_is_not_a_tail_probability():
    p,_,_=passing_curves([0.],np.zeros(6),np.ones(6),[3])
    # TCC=3 but the inclusive P(S>=3) is 42/64, not 0.5.
    assert p[3][0]==pytest.approx(42/64)
