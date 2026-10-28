"""Known forms, signed derivatives, roots and unsupported behavior."""
import numpy as np
import pytest
from scipy.special import expit
from src.irt_true_score_equating import score_roots,tcc,tcc_derivative,tcc_review,true_score_correspondence


def test_known_Rasch_form_difference_and_B_to_A_direction():
    dB=np.full(6,-.4);dA=np.full(6,.2);a=np.full(6,1.2)
    result=true_score_correspondence([3.],dB,a,dA,a,(-2,2))[0]
    assert result["theta"]==pytest.approx(1/3,abs=1e-11)
    assert result["reference_expected_score"]==pytest.approx(6*expit(.6),abs=1e-11)
    reverse=true_score_correspondence([result["reference_expected_score"]],dA,a,dB,a,(-2,2))[0]
    assert reverse["reference_expected_score"]==pytest.approx(3.,abs=1e-11)


def test_identical_form_true_score_identity():
    d=np.array([-.5,.3,.8]);a=np.array([1.2,1.,-.2])
    score=float(tcc(.2,d,a));result=true_score_correspondence([score],d,a,d,a,(-1,1))[0]
    assert result["n_roots_in_scan"]==1 and result["theta"]==pytest.approx(.2,abs=1e-10)
    assert result["reference_expected_score"]==pytest.approx(score,abs=1e-12)


def test_signed_derivative_matches_finite_difference():
    d=np.array([.3,-.5,.8]);a=np.array([1.2,-.8,0.]);theta=np.array([-2.,0.,1.])
    numerical=(tcc(theta+1e-5,d,a)-tcc(theta-1e-5,d,a))/2e-5
    np.testing.assert_allclose(tcc_derivative(theta,d,a),numerical,atol=2e-11)


def test_crossing_TCC_has_multiple_roots_and_is_withheld():
    d=np.zeros(2);a=np.array([1.,-2.]);roots=score_roots(1.1,d,a)
    assert len(roots)==2 and not tcc_review(d,a)["strictly_increasing_on_checked_grid"]
    result=true_score_correspondence([1.1],d,a,[0,0],[1,1],(-2,2))[0]
    assert np.isnan(result["reference_expected_score"]) and "multiple" in result["status"]


def test_no_root_and_tail_support_are_explicit():
    rows=true_score_correspondence([0.,1.,2.],[0,0],[1,1],[0,0],[1,1],(.5,1.5))
    assert rows[0]["n_roots_in_scan"]==0 and rows[2]["n_roots_in_scan"]==0
    assert rows[1]["status"]=="weakly supported tail/model extrapolation"
