"""Exact mid-CDF arithmetic, inverse interpolation and support conventions."""
import numpy as np
import pytest
from scipy.stats import binom
from src.irt_observed_score_equating import equate_observed,integer_interpretations,inverse_mid_cdf,mid_cdf,mid_rank_at,support_status
from src.score_distribution import marginal_distribution,reference_quadrature


def test_discrete_midcdf_and_inverse_linear_interpolation():
    pmf=np.array([.2,.5,.3]);np.testing.assert_allclose(mid_cdf(pmf),[.1,.45,.85])
    values,status=inverse_mid_cdf([.1,.275,.45,.65,.85],pmf)
    np.testing.assert_allclose(values,[0,.5,1,1.5,2]);assert (status=="available").all()
    np.testing.assert_allclose(mid_rank_at(values,pmf),[.1,.275,.45,.65,.85])


@pytest.mark.parametrize("p",[.2,.5,.8])
def test_identical_form_observed_identity_including_endpoints(p):
    pmf=binom.pmf(np.arange(7),6,p);pmf/=pmf.sum()
    result=equate_observed(np.arange(7),pmf,pmf)
    np.testing.assert_allclose(result["equated"],np.arange(7),atol=2e-14)


def test_known_mass_shift_form_difference_not_reverse_direction():
    # A = B+1 exactly in this independent PMF fixture; endpoint mass zero explicit.
    source=np.array([.2,.5,.3,0]);reference=np.array([0,.2,.5,.3])
    result=equate_observed([0,1,2],source,reference)
    np.testing.assert_allclose(result["equated"],[1,2,3],atol=1e-14)


def test_tail_mappings_never_silently_clipped():
    result=equate_observed([0,1,2],[.01,.98,.01],[.2,.6,.2])
    assert np.isnan(result["equated"][[0,2]]).all() and result["equated"][1]==pytest.approx(1.)
    assert "outside" in result["status"][0]


def test_flat_percentile_interval_is_nonunique_not_averaged():
    values,status=inverse_mid_cdf([.5],[.5,0.,0.,.5])
    assert np.isnan(values[0]) and "nonunique" in status[0]


def test_common_reference_population_applied_to_both_forms():
    d=np.array([-.5,.3,.7]);a=np.array([1.,.7,1.3])
    for population in ["A_standard_normal","B_identification_on_A","equal_identification_mixture"]:
        theta,w=reference_quadrature(population,81,1.3,.6)
        pmf=marginal_distribution(d,a,theta,w)["pmf"]
        np.testing.assert_allclose(equate_observed(np.arange(4),pmf,pmf)["equated"],np.arange(4))
    t,w=reference_quadrature("A_standard_normal",81);u,v=reference_quadrature("B_identification_on_A",81,1.3,.6)
    wrong=equate_observed([1,2],marginal_distribution(d,a,t,w)["pmf"],marginal_distribution(d,a,u,v)["pmf"])["equated"]
    assert max(abs(wrong-[1,2]))>.05


def test_integer_interpretations_and_tail_support():
    out=integer_interpretations([41.49,41.5,42.,np.nan])
    np.testing.assert_allclose(out["integer_nearest_half_up"],[41,42,42,np.nan],equal_nan=True)
    np.testing.assert_allclose(out["integer_floor"],[41,41,42,np.nan],equal_nan=True)
    assert support_status(.5,0.,0.,(-2,2))=="supported central research region"
    assert "tail" in support_status(.001,0.,0.,(-2,2))
    assert "tail" in support_status(.5,3.,0.,(-2,2))
    assert "unavailable" in support_status(.5,0.,0.,(-2,2),False)


def test_linking_parameter_change_affects_correspondence_deterministically():
    theta,w=reference_quadrature("A_standard_normal",81);d=np.linspace(-1,1,5);a=np.ones(5)
    p=marginal_distribution(d,a,theta,w)["pmf"]
    q=marginal_distribution(d-.2,a/1.1,theta,w)["pmf"]
    mapping=equate_observed(np.arange(6),q,p)["equated"]
    assert mapping[2]>2
    np.testing.assert_array_equal(mapping,equate_observed(np.arange(6),q,p)["equated"])
