"""Directional Haebara known-reference tests independent of audit observations."""
import numpy as np
import pytest
from src.linking_haebara import haebara,haebara_objective
from src.linking_stocking_lord import normal_grid,objective_gradient


@pytest.mark.parametrize("scale,shift",[(1.,0.),(1.35,-.45),(.72,.35)])
def test_haebara_affine_recovery_and_labels(scale,shift):
    d=np.linspace(-1.7,1.3,8);a=np.array([.8,-.6,.0,1e-10,1.2,.9,1.6,1.1]);theta,w=normal_grid(points=151)
    fit=haebara(d,a,d+a*shift,a*scale,theta,w)
    assert fit["method"]=="Haebara-directional-A" and fit["supported"]
    np.testing.assert_allclose([fit["scale_A"],fit["shift_B"]],[scale,shift],atol=1e-6)
    assert fit["raw_objective"]<1e-14


def test_haebara_is_item_discrepancy_not_cancelling_tcc_difference():
    theta,w=normal_grid(points=101);d=np.array([-.4,.4]);a=np.array([1.,1.]);db=d[::-1]
    sl=objective_gradient([0,0],d,a,db,a,theta,w,"Stocking-Lord")[0]
    h=haebara_objective([0,0],d,a,db,a,theta,w)[0]
    assert sl<1e-30 and h>.005
