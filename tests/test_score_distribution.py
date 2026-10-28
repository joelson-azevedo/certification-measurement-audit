"""Independent enumeration, binomial and numerical-integration fixtures."""
from itertools import product
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import binom,norm
from scipy.special import expit
from src.score_distribution import conditional_pmf,conditional_diagnostics,marginal_distribution,normal_quadrature,reference_quadrature,response_probabilities


def test_two_item_direct_enumeration():
    np.testing.assert_allclose(conditional_pmf([.2,.7]),[.24,.62,.14],atol=1e-16)


@pytest.mark.parametrize("p",[[.1,.3,.5,.8],[0.,1.,.4],[.999,.003,.6,.2,.9]])
def test_recursion_against_all_patterns(p):
    p=np.array(p);expected=np.zeros(len(p)+1)
    for pattern in product([0,1],repeat=len(p)):
        x=np.array(pattern);expected[x.sum()]+=np.prod(np.where(x,p,1-p))
    np.testing.assert_allclose(conditional_pmf(p),expected,atol=2e-16)


@pytest.mark.parametrize("slope",[[1.2,-.9,.0,1e-12],[.8,1.,1.4,.7]])
def test_signed_zero_nearzero_probability_mean_variance(slope):
    theta=np.linspace(-6,6,31);p=response_probabilities(theta,[.4,-.6,.1,.8],slope);pmf=conditional_pmf(p)
    check=conditional_diagnostics(p,pmf)
    assert check["minimum_probability"]>=0 and check["maximum_normalization_error"]<1e-14
    assert check["maximum_mean_TCC_error"]<1e-14 and check["maximum_variance_Bernoulli_error"]<1e-14


def test_sixty_identical_bernoulli_probabilities_is_binomial():
    np.testing.assert_allclose(conditional_pmf(np.full(60,.5)),binom.pmf(np.arange(61),60,.5),atol=1e-15)


@pytest.mark.parametrize("nodes",[21,81,161])
def test_GH_normalization_and_normal_moments(nodes):
    theta,w=normal_quadrature(nodes,.4,1.3)
    assert w.sum()==pytest.approx(1,abs=2e-15)
    assert w@theta==pytest.approx(.4,abs=1e-14)
    assert w@(theta*theta)==pytest.approx(.4**2+1.3**2,abs=1e-13)


def test_marginal_mean_variance_and_integral_by_quad():
    d=np.array([.3,-.7,1.]);a=np.array([1.1,-.9,0.]);theta,w=normal_quadrature(81)
    result=marginal_distribution(d,a,theta,w)
    expected=0.
    for intercept,slope in zip(d,a):expected+=quad(lambda z:expit(intercept+slope*z)*norm.pdf(z),-10,10,epsabs=1e-12)[0]
    assert result["diagnostics"]["marginal_mean"]==pytest.approx(expected,abs=1e-12)
    assert result["diagnostics"]["marginal_variance_total_variance_error"]<1e-13
    assert np.diff(result["cdf"]).min()>=0
    np.testing.assert_array_equal(result["pmf"],marginal_distribution(d,a,theta,w)["pmf"])


def test_mixture_is_equal_weight_distribution_not_pooled_group_estimate():
    theta,w=reference_quadrature("equal_identification_mixture",61,1.2,-.3)
    assert w@theta==pytest.approx(-.15)
    assert w@(theta*theta)-(w@theta)**2==pytest.approx((1+1.2**2)/2+.3**2/4)


@pytest.mark.parametrize("p",[[1.1,.4],[-.01,.4],[np.nan,.4]])
def test_invalid_probabilities_are_rejected_without_clipping(p):
    with pytest.raises(ValueError):conditional_pmf(p)


def test_reference_weights_must_be_common_and_normalized():
    with pytest.raises(ValueError):marginal_distribution([0],[1],[0,1],[.5,.4])
