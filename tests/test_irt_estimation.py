"""MML recovery on fixtures created for testing, never the case's hidden truth."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.special import expit

from src.irt_estimation import EstimationConfig,MarginalLikelihood,detect_boundaries,fit_mml,parameter_uncertainty


@pytest.fixture(scope='module')
def signed_fixture():
    rng=np.random.default_rng(202640)
    n=2500
    d=np.array([-1.1,-.4,0,.5,1.,-.6,.1,.2])
    a=np.array([.8,1.1,.95,1.4,.65,1.2,-1.,.03])
    theta=rng.standard_normal(n)
    x=(rng.random((n,len(a)))<expit(d+theta[:,None]*a)).astype(float)
    x[::29,1]=np.nan
    return x,d,a


@pytest.fixture(scope='module')
def fitted_fixture(signed_fixture):
    return fit_mml(signed_fixture[0])


def test_analytic_gradient_matches_finite_differences(signed_fixture):
    x,d,a=signed_fixture;likelihood=MarginalLikelihood(x[:40],31);parameters=np.r_[d,a]
    _,gradient=likelihood.value_gradient(parameters)
    numerical=[]
    for j in range(len(parameters)):
        step=np.zeros_like(parameters);step[j]=1e-5
        numerical.append((likelihood.value_gradient(parameters+step)[0]-likelihood.value_gradient(parameters-step)[0])/(2e-5))
    np.testing.assert_allclose(gradient,numerical,atol=1e-6,rtol=1e-6)


def test_louis_observed_information_matches_gradient_derivatives(signed_fixture):
    x,d,a=signed_fixture;likelihood=MarginalLikelihood(x[:50],31);parameters=np.r_[d,a]
    actual=likelihood.observed_information(parameters)
    columns=[]
    for j in range(len(parameters)):
        step=np.zeros_like(parameters);step[j]=1e-5
        columns.append((likelihood.value_gradient(parameters+step)[1]-likelihood.value_gradient(parameters-step)[1])/(2e-5))
    np.testing.assert_allclose(actual,np.column_stack(columns),atol=1e-6,rtol=1e-6)


def test_parameter_recovery_on_controlled_signed_fixture(signed_fixture,fitted_fixture):
    _,d,a=signed_fixture;fit=fitted_fixture
    assert fit['diagnostics']['numerically_admissible']
    np.testing.assert_allclose(fit['intercept'],d,atol=.25)
    np.testing.assert_allclose(fit['slope'],a,atol=.35)
    assert fit['slope'][6]<-.5
    assert fit['slope'].sum()>0


def test_deterministic_estimation_and_first_order_check(signed_fixture,fitted_fixture):
    again=fit_mml(signed_fixture[0])
    np.testing.assert_allclose(again['parameters'],fitted_fixture['parameters'],atol=1e-12,rtol=0)
    assert again['diagnostics']['projected_gradient_max_abs_per_candidate']<=1e-6
    assert again['start_diagnostics'].converged.all()


def test_supported_information_and_standard_errors(fitted_fixture):
    fit=fitted_fixture
    assert fit['diagnostics']['uncertainty_available']
    assert fit['diagnostics']['information_minimum_eigenvalue']>0
    assert fit['diagnostics']['maximum_directional_hessian_relative_error']<1e-5
    np.testing.assert_allclose(fit['observed_information']@fit['covariance'],np.eye(len(fit['parameters'])),atol=1e-9)


def test_singular_information_is_not_blindly_inverted(signed_fixture,monkeypatch):
    x,d,a=signed_fixture;likelihood=MarginalLikelihood(x[:50]);parameters=np.r_[d,a]
    monkeypatch.setattr(likelihood,'observed_information',lambda _:np.zeros((len(parameters),len(parameters))))
    covariance,_,diagnostic=parameter_uncertainty(likelihood,parameters,EstimationConfig(),np.zeros(len(parameters),bool))
    assert not diagnostic['uncertainty_available'] and np.isnan(covariance).all()


def test_optimizer_failure_is_reported_not_fabricated(signed_fixture):
    config=replace(EstimationConfig(),max_iterations=0)
    fit=fit_mml(signed_fixture[0][:150],config)
    assert not fit['diagnostics']['converged'] and not fit['diagnostics']['numerically_admissible']
    assert not fit['diagnostics']['uncertainty_available']
    assert 'not established' in fit['diagnostics']['warnings']


def test_boundary_hits_suppress_unconstrained_uncertainty(signed_fixture):
    config=replace(EstimationConfig(),slope_bound=.2)
    fit=fit_mml(signed_fixture[0][:300],config)
    assert fit['diagnostics']['parameter_boundary_hits']>0
    assert not fit['diagnostics']['uncertainty_available'] and not fit['diagnostics']['numerically_admissible']
    mask=detect_boundaries(np.r_[np.zeros(3),[.2,-.2,0]],3,config)
    np.testing.assert_equal(mask,[False,False,False,True,True,False])


def test_no_response_variation_stops_parameter_estimation():
    with pytest.raises(ValueError,match='variation'):
        fit_mml(np.ones((50,4)))


def test_likelihood_does_not_modify_input(signed_fixture):
    x=signed_fixture[0].copy();before=x.copy()
    fit_mml(x)
    np.testing.assert_equal(x,before)
