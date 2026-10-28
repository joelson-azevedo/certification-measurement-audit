"""Joint covariance, within-form resampling, stopping and retained failure tests."""
import numpy as np
import pandas as pd
import pytest
from src.irt_estimation import EstimationConfig
from src.linking_stocking_lord import LinkingConfig,normal_grid
from src.linking_uncertainty import covariance_diagnostics,experiment,joint_parameter_draw,percentile_mc_bounds,summarize_experiment


def models():
    d=np.array([-.8,.1,.7]);a=np.array([.8,1.2,1.5]);c=np.array([[.2,0,.04,0,.03,0],[0,.2,0,.01,0,.02],[.04,0,.2,0,.01,0],[0,.01,0,.1,0,.02],[.03,0,.01,0,.1,0],[0,.02,0,.02,0,.1]])*.05
    return {f:{"form":f,"items":["I1","I2","I3"],"d":d,"a":a,"covariance":c,"x":np.zeros((15,3))} for f in ["A","B"]}


def test_full_covariance_seed_exact_and_independent_form_streams():
    m=models()["A"];seed=871;attempt=3
    expected=np.r_[m["d"],m["a"]]+np.linalg.cholesky(m["covariance"])@np.random.default_rng(np.random.SeedSequence([seed,attempt,0])).standard_normal(6)
    np.testing.assert_array_equal(joint_parameter_draw(m,seed,attempt,0),expected)
    assert not np.array_equal(expected,joint_parameter_draw(m,seed,attempt,1))
    assert covariance_diagnostics(m)["covariance_supported"]
    # Cross-item and d/a covariance are both nonzero and affect the exact draw.
    diagonal=np.r_[m["d"],m["a"]]+np.sqrt(np.diag(m["covariance"]))*np.random.default_rng(np.random.SeedSequence([seed,attempt,0])).standard_normal(6)
    assert not np.allclose(expected,diagonal)


def test_covariance_experiment_reproducibility_and_intervals():
    m=models();theta,w=normal_grid(points=81);points={"Stocking-Lord":{"scale_A":1.,"shift_B":0.},"Haebara-directional-A":{"scale_A":1.,"shift_B":0.}}
    args=(m,m["A"]["items"],points,theta,w,LinkingConfig(),EstimationConfig(),"joint_covariance_normal",4,6,12,0.)
    first=experiment(*args);second=experiment(*args)
    np.testing.assert_array_equal(first["arrays"]["parameters_A"],second["arrays"]["parameters_A"])
    np.testing.assert_array_equal(first["draws"].scale_A,second["draws"].scale_A)
    summary=summarize_experiment(first,4)
    assert summary.uncertainty_supported.all() and summary.valid_draws.eq(4).all()
    for row in summary.itertuples():
        values=first["draws"].loc[first["draws"].method.eq(row.method),"scale_A"]
        assert row.scale_A_se==pytest.approx(values.std(ddof=1))
        assert row.scale_A_95_low==pytest.approx(np.quantile(values,.025))


def test_bootstrap_recalibration_failures_and_indices_are_preserved(monkeypatch):
    import src.linking_uncertainty as module
    def fail(*args,**kwargs):raise ValueError("controlled failed calibration")
    monkeypatch.setattr(module,"fit_mml",fail)
    m=models();theta,w=normal_grid(points=31);points={"Stocking-Lord":{"scale_A":1.,"shift_B":0.}}
    result=experiment(m,m["A"]["items"],points,theta,w,LinkingConfig(),EstimationConfig(),"candidate_bootstrap",2,3,19,0.)
    assert len(result["draws"])==3 and not result["draws"].valid_paired_draw.any()
    assert len(result["calibrations"])==6 and result["draws"].failure.str.contains("controlled failed calibration").all()
    expected=np.random.default_rng(np.random.SeedSequence([19,0,0])).integers(0,15,15)
    np.testing.assert_array_equal(result["arrays"]["sample_indices_A"][0],expected)
    assert not np.array_equal(result["arrays"]["sample_indices_A"][0],result["arrays"]["sample_indices_B"][0])
    summary=summarize_experiment(result,2)
    assert not summary.uncertainty_supported.any() and summary.failed_draws.eq(3).all()


def test_link_optimizer_failure_retained_and_not_counted(monkeypatch):
    import src.linking_uncertainty as module
    original=module.fit_link
    def fail(*args,**kwargs):
        fit=original(*args,**kwargs);fit["supported"]=False;return fit
    monkeypatch.setattr(module,"fit_link",fail)
    m=models();theta,w=normal_grid(points=31)
    result=experiment(m,m["A"]["items"],{"Stocking-Lord":{"scale_A":1.,"shift_B":0.}},theta,w,LinkingConfig(),EstimationConfig(),"joint_covariance_normal",2,3,29,0.)
    assert len(result["draws"])==3 and not result["draws"].valid_paired_draw.any()
    assert np.isfinite(result["arrays"]["parameters_A"]).all() and len(result["optimizers"])==6


def test_percentile_mc_bounds_and_no_individual_slope_clipping():
    values=np.arange(1,201)
    low,high=percentile_mc_bounds(values,.025)
    assert low<=np.quantile(values,.025)<=high
    m=models()["A"];m["a"][0]=-1.
    draw=joint_parameter_draw(m,41,0,0)
    assert draw[3]<0


def test_covariance_cache_reuses_empty_calibration_table_on_windows(tmp_path):
    from pathlib import Path
    from src.linking_uncertainty import cached_experiment
    root=Path(__file__).resolve().parents[1]/"analysis";m=models();theta,w=normal_grid(points=31)
    points={"Stocking-Lord":{"scale_A":1.,"shift_B":0.}}
    args=(root,tmp_path,m,m["A"]["items"],points,theta,w,LinkingConfig(),EstimationConfig(),"joint_covariance_normal",2,3,194,0.)
    original=cached_experiment(*args);messages=[];cached=cached_experiment(*args,progress=messages.append)
    assert any("verified exact-input/source cache" in text for text in messages)
    assert cached["calibrations"].empty
    np.testing.assert_array_equal(original["arrays"]["parameters_A"],cached["arrays"]["parameters_A"])
    np.testing.assert_array_equal(original["draws"].scale_A,cached["draws"].scale_A)
