"""Prespecified rounding/boundary alternatives and transition bookkeeping."""
import numpy as np
import pytest
from src.equating_decision_impact import aggregate_impact,counterfactual_by_score

CONVENTIONS={"continuous":{"reference_score_boundary":2.,"comparison":">="},
             "floor":{"reference_score_boundary":2.,"comparison":">="},
             "nearest":{"reference_score_boundary":1.5,"comparison":">="},
             "ceiling":{"reference_score_boundary":1.,"comparison":">"}}


def test_identity_boundary_and_fractional_rounding_are_explicit():
    p=np.array([.1,.2,.4,.3]);rules=counterfactual_by_score(p,p,CONVENTIONS)
    for _,group in rules.groupby("convention"):
        np.testing.assert_array_equal(group.counterfactual_pass,[False,False,True,True])
    # Mid-percentile thresholds correspond exactly to e>=2, e>=1.5, e>1.
    q=np.array([.2,.4,.3,.1]);rules=counterfactual_by_score(q,p,CONVENTIONS)
    assert rules.loc[rules.convention.eq("continuous"),"counterfactual_pass"].tolist()==[False,False,True,True]
    assert rules.loc[rules.convention.eq("nearest"),"counterfactual_pass"].tolist()==[False,True,True,True]


def test_aggregate_counts_transitions_and_original_decisions_remain_separate():
    raw=np.array([0,1,1,2,3]);recorded=raw>=2
    original=raw.copy();decisions=recorded.copy()
    overall,strata,_=aggregate_impact(raw,recorded,[.2,.4,.3,.1],[.1,.2,.4,.3],CONVENTIONS,cut=2)
    nearest=overall.loc[overall.convention.eq("nearest")].iloc[0]
    assert nearest.original_pass==2 and nearest.counterfactual_pass==4 and nearest.upward_changes==2 and nearest.downward_changes==0
    for row in overall.itertuples():
        assert row.original_pass+row.original_fail==5 and row.counterfactual_pass+row.counterfactual_fail==5
        assert row.counterfactual_pass==row.original_pass+row.upward_changes-row.downward_changes
        group=strata.loc[strata.convention.eq(row.convention)];assert group.upward_changes.sum()==row.upward_changes
    np.testing.assert_array_equal(raw,original);np.testing.assert_array_equal(recorded,decisions)


def test_tail_threshold_ordering_does_not_require_fabricated_tail_scores():
    rules=counterfactual_by_score([.01,.98,.01],[.2,.6,.2],{"primary":{"reference_score_boundary":1.,"comparison":">="}})
    np.testing.assert_array_equal(rules.counterfactual_pass,[False,True,True])


def test_recorded_rule_disagreement_is_rejected():
    with pytest.raises(ValueError):aggregate_impact([0,1],[True,True],[.4,.6],[.4,.6],{"primary":{"reference_score_boundary":1.,"comparison":">="}},cut=1)
