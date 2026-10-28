"""Joint probabilities, common-theta integration and independent response simulation."""
import numpy as np
import pytest
from src.decision_consistency import integrated_models,integrate_joint,joint_classification,joint_statistics,monte_carlo_fixture
from src.score_distribution import reference_quadrature


@pytest.mark.parametrize("p",[0.,.1,.5,.9,1.])
def test_same_form_identities_and_extremes(p):
    joint=joint_classification(p)
    assert joint["agreement"]==pytest.approx(p*p+(1-p)**2)
    assert joint["disagreement"]==pytest.approx(2*p*(1-p))
    assert joint["agreement"]+joint["disagreement"]==pytest.approx(1)
    assert joint["both_pass"]==pytest.approx(p*p)
    assert joint["both_fail"]==pytest.approx((1-p)**2)
    if p==.5:assert joint["agreement"]==.5


def test_cross_form_asymmetry_and_directional_sum():
    j=joint_classification(.8,.2)
    assert j["agreement"]==pytest.approx(.32)
    assert j["first_pass_second_fail"]==pytest.approx(.64)
    assert j["first_fail_second_pass"]==pytest.approx(.04)
    assert j["disagreement"]==pytest.approx(.68)
    swapped=joint_classification(.2,.8)
    assert swapped["agreement"]==pytest.approx(j["agreement"])
    assert swapped["first_pass_second_fail"]==pytest.approx(j["first_fail_second_pass"])


def test_common_latent_integral_is_not_product_of_marginal_pass_rates():
    p=np.array([.1,.9]);weights=np.array([.5,.5]);result=integrate_joint(p,p,weights)
    assert result["agreement"]==pytest.approx(.82)
    assert result["both_pass"]==pytest.approx(.41)
    assert result["first_pass"]==pytest.approx(.5)
    assert result["chance_agreement_from_marginals"]==pytest.approx(.5)
    assert result["kappa"]==pytest.approx(.64)
    assert result["both_pass"]!=pytest.approx(result["first_pass"]*result["second_pass"])


def test_positive_negative_agreement_and_kappa_definitions():
    result=joint_statistics(dict(both_pass=.4,first_pass_second_fail=.1,first_fail_second_pass=.2,both_fail=.3))
    assert result["positive_agreement"]==pytest.approx(.8/1.1)
    assert result["negative_agreement"]==pytest.approx(.6/.9)
    assert result["kappa"]==pytest.approx(.4)
    degenerate=joint_statistics(dict(both_pass=1.,first_pass_second_fail=0.,first_fail_second_pass=0.,both_fail=0.))
    assert np.isnan(degenerate["kappa"]) and np.isnan(degenerate["negative_agreement"])
    assert "undefined" in degenerate["kappa_status"]


def test_reference_population_integration_and_identical_form_symmetry():
    d=np.zeros(60);a=np.ones(60);scenarios=[{"scenario":"O","cut_A":42,"cut_B":42}]
    result,details,log=integrated_models(d,a,d,a,1.,0.,"A_standard_normal",scenarios,[81,161,321,641])
    same=result.loc[result.form.eq("A")].iloc[0];cross=result.loc[result.pair_type.eq("cross_form")].iloc[0]
    assert same.agreement==pytest.approx(cross.agreement)
    assert cross.first_pass_second_fail==pytest.approx(cross.first_fail_second_pass)
    assert log.converged.iloc[-1] and details["weights"].sum()==pytest.approx(1)
    t,w=reference_quadrature("B_identification_on_A",81,1.2,.4)
    assert w@t==pytest.approx(.4)


def test_failed_quadrature_does_not_fabricate_a_result():
    d=np.zeros(60);a=np.ones(60)*2
    with pytest.raises(ValueError,match="quadrature"):
        integrated_models(d,a,d-.4,a,1.,0.,"A_standard_normal",[{"scenario":"O","cut_B":42}],[5,7],1e-20,1e-20)


def test_independent_monte_carlo_responses_match_both_theoretical_cases():
    results=monte_carlo_fixture(draws=120000,seed=631)
    assert results.passed.all() and set(results.fixture)=={"same_form","cross_form"}
    assert (results.Monte_Carlo_SE>0).all()
    for row in results.itertuples():assert row.absolute_difference<=row.tolerance
