"""Known-affine controlled fixtures, objective scale and analytic gradients."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.special import expit
from src.linking_stocking_lord import LinkingConfig,fit_link,normal_grid,objective_gradient


@pytest.mark.parametrize("scale,shift",[(1.,0.),(1.3,-.4),(.8,.3)])
def test_stocking_lord_recovers_known_affine(scale,shift):
    d=np.linspace(-1.6,1.5,9);a=np.array([.7,1.,1.4,-.8,.9,1.2,1.7,.6,1.1])
    theta,w=normal_grid(points=201)
    fit=fit_link(d,a,d+a*shift,a*scale,theta,w)
    assert fit["supported"] and fit["method"]=="Stocking-Lord"
    np.testing.assert_allclose([fit["scale_A"],fit["shift_B"]],[scale,shift],atol=1e-6)
    assert fit["raw_objective"]<1e-14


@pytest.mark.parametrize("method",["Stocking-Lord","Haebara-directional-A"])
def test_objective_matches_direct_sum_and_gradient(method):
    d=np.array([-.4,.8,1.]);a=np.array([1.,1.5,-.7]);db=d+[.2,-.1,.3];ab=a*[.8,1.1,1.4]
    theta,w=normal_grid(points=151);z=np.array([.12,-.3])
    delta=expit(d+theta[:,None]*a)-expit(db+((theta-z[1])/np.exp(z[0]))[:,None]*ab)
    manual=w@(delta.sum(axis=1)**2)/9 if method=="Stocking-Lord" else np.sum(w[:,None]*delta**2)/3
    value,gradient=objective_gradient(z,d,a,db,ab,theta,w,method)
    assert value==pytest.approx(manual,abs=1e-16)
    finite=[]
    for j in range(2):
        step=np.eye(2)[j]*1e-6
        finite.append((objective_gradient(z+step,d,a,db,ab,theta,w,method)[0]-objective_gradient(z-step,d,a,db,ab,theta,w,method)[0])/2e-6)
    np.testing.assert_allclose(gradient,finite,atol=1e-11,rtol=1e-7)


def test_optimizer_determinism_and_failure_detection():
    d=np.linspace(-1,1,7);a=np.linspace(.7,1.4,7);theta,w=normal_grid(points=101)
    first=fit_link(d,a,d+.4*a,1.3*a,theta,w);second=fit_link(d,a,d+.4*a,1.3*a,theta,w)
    assert first["scale_A"]==second["scale_A"] and first["shift_B"]==second["shift_B"]
    failed=fit_link(d,a,d+.4*a,1.3*a,theta,w,replace(LinkingConfig(),max_iterations=0),starts=[(1.,0.)])
    assert not failed["supported"] and not failed["attempts"].supported.any()
    boundary=fit_link(d,a,d+4*a,a,theta,w)
    assert not boundary["supported"] and boundary["boundary_hit"]


def test_normal_grid_and_malformed_correspondence():
    theta,w=normal_grid(points=101)
    assert w.sum()==pytest.approx(1) and np.all(w>0)
    np.testing.assert_allclose(w,w[::-1]);assert theta[0]==-3 and theta[-1]==3
    with pytest.raises(ValueError):normal_grid(points=2)
    with pytest.raises(ValueError):objective_gradient([0,0],[1],[1,2],[1],[1],theta,w)
    with pytest.raises(ValueError):objective_gradient([0,0],[1],[1],[1],[1],theta,w*2)
