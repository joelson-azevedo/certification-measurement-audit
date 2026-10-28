"""Independent likelihood, DIF-recovery and effect-size reference fixtures."""
import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2

from src.dif import (LogisticConfig,bh_adjust,common_support,designs,effect_sizes,fit_logistic,
                      fit_nested,logistic_loglikelihood,probability_contrasts,separation_diagnostic)


def test_logistic_likelihood_direct_and_extreme_reference():
    x=np.array([[1,-1],[1,0],[1,1]],float);y=np.array([0,1,1]);b=np.array([.3,.8])
    p=expit(x@b)
    assert logistic_loglikelihood(b,x,y)==pytest.approx(np.log((1-p[0])*p[1]*p[2]))
    assert np.isfinite(logistic_loglikelihood(np.array([0,1000.]),x,y))


@pytest.mark.parametrize('group_effect,interaction',[(.8,0),(.1,.9)])
def test_uniform_and_nonuniform_dif_recovery(group_effect,interaction):
    rng=np.random.default_rng(202610260)
    t=rng.normal(size=6000);g=rng.integers(0,2,len(t));t-=t.mean()
    beta=np.array([-.2,.9,group_effect,interaction]);y=rng.binomial(1,expit(designs(t,g)[2]@beta))
    fits,tests=fit_nested(t,g,y)
    assert tests['all_nested_models_supported']
    np.testing.assert_allclose(fits[2]['beta'],beta,atol=.13)
    assert tests['p_omnibus']<1e-8
    key='p_uniform' if interaction==0 else 'p_nonuniform'
    assert tests[key]<1e-8
    assert fits[0]['log_likelihood']<=fits[1]['log_likelihood']+1e-8<=fits[2]['log_likelihood']+1e-8
    assert tests['p_omnibus']==pytest.approx(chi2.sf(2*(fits[2]['log_likelihood']-fits[0]['log_likelihood']),2))


def test_no_dif_identical_conditional_frequencies():
    rows=[]
    for g in [0,1]:
        for t,correct in [(-2,12),(-1,27),(0,50),(1,73),(2,88)]:
            rows.extend([(t,g,1)]*correct+[(t,g,0)]*(100-correct))
    table=np.array(rows);fits,tests=fit_nested(table[:,0],table[:,1],table[:,2])
    assert tests['p_omnibus']==pytest.approx(1,abs=1e-9)
    assert tests['p_uniform']==pytest.approx(1,abs=1e-6)
    np.testing.assert_allclose(fits[2]['beta'][2:],0,atol=1e-8)


def test_mle_against_independent_scipy_objective_optimization():
    rng=np.random.default_rng(710);x=np.column_stack([np.ones(800),rng.normal(size=800),rng.integers(0,2,800)])
    y=rng.binomial(1,expit(x@np.array([.2,.7,-.5])))
    fit=fit_logistic(x,y)
    def objective(b):
        eta=x@b
        return -np.sum(y*(-np.logaddexp(0,-eta))+(1-y)*(-np.logaddexp(0,eta)))
    ref=minimize(objective,np.zeros(3),method='BFGS',options={'gtol':1e-6})
    np.testing.assert_allclose(fit['beta'],ref.x,atol=3e-5)
    assert -fit['log_likelihood']==pytest.approx(ref.fun,abs=1e-8)
    p=expit(x@fit['beta']);h=sum(p[i]*(1-p[i])*np.outer(x[i],x[i]) for i in range(len(x)))
    np.testing.assert_allclose(h@fit['covariance'],np.eye(3),atol=1e-12)


@pytest.mark.parametrize('t,y,status',[
    ([-2,-1,1,2],[0,0,1,1],'complete separation'),
    ([-1,0,0,1],[0,0,1,1],'quasi-complete separation')])
def test_separation_is_detected_and_ml_withheld(t,y,status):
    x=np.column_stack([np.ones(len(t)),t]);y=np.array(y)
    assert separation_diagnostic(x,y)['separation_status']==status
    fit=fit_logistic(x,y)
    assert not fit['numerically_supported'] and np.isnan(fit['beta']).all() and np.isnan(fit['log_likelihood'])


def test_failure_and_rank_deficiency_are_unavailable():
    x=np.column_stack([np.ones(8),[-2,-2,-1,-1,1,1,2,2]])
    y=np.array([0,1,0,1,0,1,0,1])
    fit=fit_logistic(x,y,LogisticConfig(max_iterations=0))
    assert not fit['converged'] and np.isnan(fit['log_likelihood'])
    singular=fit_logistic(np.ones((8,2)),y)
    assert 'rank-deficient' in singular['status']


def test_bh_known_values_and_unavailable_family_slots():
    np.testing.assert_allclose(bh_adjust([.01,.04,.03,.002]),[.02,.04,.04,.008])
    np.testing.assert_allclose(bh_adjust([.01,.04,np.nan]),[.03,.06,np.nan],equal_nan=True)
    np.testing.assert_allclose(bh_adjust([.001,.01],family_size=60),[.06,.3])
    with pytest.raises(ValueError):bh_adjust([-.1,.5])


def test_common_support_local_counts_and_reference_weights():
    rng=np.random.default_rng(442)
    t=np.r_[rng.uniform(-2,2,1000),rng.uniform(-1,3,1000)];g=np.repeat([0,1],1000)
    support=common_support(t,g)
    assert support['available'] and support['lower']>=np.quantile(t[g==1],.05)-1e-12
    assert support['upper']<=np.quantile(t[g==0],.95)+1e-12
    assert support['reference_weights'].sum()==pytest.approx(1)
    assert support['reference_n_south']>0 and support['reference_n_north']>0
    assert support['bins'].loc[support['bins'].supported,['n_south','n_north']].ge(20).all().all()
    absent=common_support(np.r_[np.arange(100),np.arange(100)+200],np.repeat([0,1],100))
    assert not absent['available']


def test_probability_orientation_and_delta_covariance():
    b=np.array([0,1,.7,0]);cov=np.eye(4)*.02;t=np.array([-1.,0,1])
    result=probability_contrasts(b,cov,t)
    assert result.difference_north_minus_south.gt(0).all()
    np.testing.assert_allclose(result.p_south,expit(t))
    np.testing.assert_allclose(result.p_north,expit(t+.7))
    step=1e-5;grad=[]
    for j in range(4):
        plus=b.copy();minus=b.copy();plus[j]+=step;minus[j]-=step
        dp=expit(plus[0]+plus[1]*t+plus[2]+plus[3]*t)-expit(plus[0]+plus[1]*t)
        dm=expit(minus[0]+minus[1]*t+minus[2]+minus[3]*t)-expit(minus[0]+minus[1]*t)
        grad.append((dp-dm)/(2*step))
    grad=np.array(grad).T
    np.testing.assert_allclose(result.difference_se_delta,np.sqrt(np.einsum('ij,jk,ik->i',grad,cov,grad)),atol=1e-10)


def test_unsigned_effect_preserves_crossing_and_nagelkerke_definition():
    beta=np.array([0,1,.5,-1.]);cov=np.eye(4)*.01
    fits=[{'beta':beta[:k],'covariance':cov[:k,:k],'numerically_supported':True,'log_likelihood':ll}
          for k,ll in [(2,-100),(3,-98),(4,-96)]]
    t=np.array([-1.,0,1]);w=np.array([.25,.5,.25])
    support={'available':True,'reference_values':t,'reference_weights':w,
             'grid':pd.DataFrame({'matching_centered':t,'supported':[True]*3}), 'quantile_values':t}
    y=np.r_[np.zeros(100),np.ones(100)]
    effect,grid=effect_sizes(fits,y,support)
    delta=expit(.5+0*t)-expit(t)
    assert effect['signed_reference_mean_difference']==pytest.approx(w@delta)
    assert effect['absolute_reference_mean_difference']==pytest.approx(w@abs(delta))
    assert effect['absolute_reference_mean_difference']>abs(effect['signed_reference_mean_difference'])
    assert effect['crossing_curves_on_supported_grid']
    ll_null=-200*np.log(2)
    r0=(1-np.exp(2*(ll_null+100)/200))/(1-np.exp(2*ll_null/200))
    r2=(1-np.exp(2*(ll_null+96)/200))/(1-np.exp(2*ll_null/200))
    assert effect['delta_nagelkerke_omnibus']==pytest.approx(r2-r0)


def test_crossing_is_not_claimed_across_an_unsupported_gap():
    beta=np.array([0,1,0,1.]);cov=np.eye(4)*.01
    fits=[{'beta':beta[:k],'covariance':cov[:k,:k],'numerically_supported':True,'log_likelihood':ll}
          for k,ll in [(2,-100),(3,-98),(4,-96)]]
    support={'available':True,'reference_values':np.array([-1.,1]),'reference_weights':np.array([.5,.5]),
             'grid':pd.DataFrame({'matching_centered':[-1.,0,1], 'supported':[True,False,True]}),
             'quantile_values':np.array([-1.,1,1])}
    effect,grid=effect_sizes(fits,np.r_[np.zeros(100),np.ones(100)],support)
    assert not effect['crossing_curves_on_supported_grid']
    assert np.isnan(effect['matching_at_supported_curve_crossing'])
    assert grid.loc[1,'difference_north_minus_south'] is np.nan or np.isnan(grid.loc[1,'difference_north_minus_south'])
