"""ID-based correspondence, broad content and held-out curve independence."""
import json
import numpy as np
import pandas as pd
import pytest
from src.linking_anchor_diagnostics import align_anchors,anchor_sets,content_coverage,leave_one_out,supported_range
from src.linking_stocking_lord import LinkingConfig,fit_link,normal_grid


def fixture():
    ids=[f"ANC_{i:02}" for i in range(1,21)]
    meta=pd.DataFrame({"item_id":ids,"anchor_item":True,"domain":["D"+str(i%5) for i in range(20)]})
    d=np.linspace(-1.3,1.7,20);a=np.linspace(.7,1.5,20);a[13]=-.8
    base={"items":ids,"d":d,"a":a,"eap":np.linspace(-2,2,150),"covariance":np.eye(40)*.01}
    b={**base,"items":ids[::-1],"d":(d+a*.2)[::-1],"a":(a*1.2)[::-1],"eap":(base["eap"]-.2)/1.2}
    return meta,base,b


def test_id_matching_and_prespecified_exclusions():
    meta,a,b=fixture();sets=anchor_sets(meta,a["items"],b["items"])
    assert [len(sets[s]) for s in ["S20","S19","S18"]]==[20,19,18]
    assert "ANC_14" in sets["S20"] and "ANC_14" not in sets["S19"] and "ANC_18" not in sets["S18"]
    da,aa,db,ab=align_anchors(a,b,["ANC_03","ANC_01"])
    np.testing.assert_allclose(db,da+aa*.2);np.testing.assert_allclose(ab,aa*1.2)
    with pytest.raises(ValueError):align_anchors(a,b,["ANC_01","ANC_01"])
    with pytest.raises(ValueError):align_anchors(a,b,["absent"])
    bad=meta.copy();bad.loc[0,"anchor_item"]=False
    with pytest.raises(ValueError):anchor_sets(bad,a["items"],b["items"])


def test_domain_coverage_and_support_intersection():
    meta,a,b=fixture();sets=anchor_sets(meta,a["items"],b["items"])
    coverage=content_coverage(meta,{"A":a,"B":b},sets)
    assert coverage.loc[coverage.anchor_specification.eq("S18")&coverage.form.eq("A"),"anchor_count"].sum()==18
    assert coverage.anchor_count.min()>=2
    np.testing.assert_allclose(supported_range(a,b,1.2,.2),[-1.96,1.96])


def test_leave_one_out_excludes_focal_and_evaluates_heldout_curve():
    meta,a,b=fixture();sets=anchor_sets(meta,a["items"],b["items"]);theta,w=normal_grid(points=101)
    # Perturb one B item: its own held-out link must recover the known transform.
    b["d"]=b["d"].copy();b["d"][b["items"].index("ANC_05")]+=.4
    point=fit_link(*align_anchors(a,b,sets["S18"]),theta,w)
    loo,logs=leave_one_out(a,b,sets["S18"],meta,{"Stocking-Lord":point},theta,w,(-1.5,1.5),LinkingConfig(),{"minimum_anchors":15,"minimum_anchors_per_domain":2})
    assert len(loo)==18 and len(logs)==90
    for row in loo.itertuples():assert row.held_out_item not in json.loads(row.retained_anchor_ids) and row.n_retained_anchors==17
    held=loo.loc[loo.held_out_item.eq("ANC_05")].iloc[0]
    np.testing.assert_allclose([held.scale_A,held.shift_B],[1.2,.2],atol=1e-6)
    from scipy.special import expit
    mask=abs(theta)<=1.5;weight=w[mask]/w[mask].sum();j=a["items"].index("ANC_05")
    delta=expit(a["d"][j]+.4+a["a"][j]*theta[mask])-expit(a["d"][j]+a["a"][j]*theta[mask])
    assert held.held_out_integrated_absolute_probability_difference==pytest.approx(weight@abs(delta),abs=1e-8)
