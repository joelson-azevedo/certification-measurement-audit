"""Binary reference marginals, reproducibility and controlled factor structures."""

import numpy as np
import pytest

from src.parallel_analysis import parallel_analysis,permute_columns
from src.tetrachoric import tetrachoric_matrix


def test_reference_marginals_and_missing_positions_preserved_exactly():
    rng=np.random.default_rng(10)
    x=rng.binomial(1,np.array([.2,.4,.6,.8]),size=(100,4)).astype(float)
    x[::9,1]=np.nan
    null=permute_columns(x,np.random.default_rng(72))
    np.testing.assert_equal(np.isnan(null),np.isnan(x))
    np.testing.assert_equal(np.nansum(null,axis=0),np.nansum(x,axis=0))
    assert not np.array_equal(x,np.nan_to_num(null))


def test_parallel_analysis_reproducibility():
    rng=np.random.default_rng(11)
    x=rng.binomial(1,.5,size=(300,6))
    observed,_,_=tetrachoric_matrix(x)
    a,ar,logs=parallel_analysis(x,observed,list('abcdef'),replicates=8,seed=41)
    b,br,_=parallel_analysis(x,observed,list('abcdef'),replicates=8,seed=41)
    np.testing.assert_equal(ar,br)
    assert a.equals(b) and logs.valid.all() and logs.margins_exact.all()
    assert a.valid_replications.eq(8).all()
    expected=np.percentile(ar,95,axis=0)
    np.testing.assert_allclose(a.reference_percentile,expected)


@pytest.mark.parametrize("factors",[1,2])
def test_binary_structures_reasonably_recovered_without_case_truth(factors):
    rng=np.random.default_rng(900+factors)
    n,p=2000,12
    abilities=rng.standard_normal((n,factors))
    factor=np.arange(p)//6 if factors==2 else np.zeros(p,dtype=int)
    x=np.column_stack([(.8*abilities[:,factor[j]]+.6*rng.standard_normal(n)>0).astype(int) for j in range(p)])
    r,_,_=tetrachoric_matrix(x)
    table,_,_=parallel_analysis(x,r,[str(j) for j in range(p)],replicates=12,seed=456)
    assert table.suggested_dimensions.iloc[0]==factors


def test_invalid_reference_estimation_stops_instead_of_using_phi():
    x=np.zeros((50,4))
    with pytest.raises(RuntimeError,match="valid tetrachoric references"):
        parallel_analysis(x,np.eye(4),list('abcd'),replicates=2,max_attempts=3)


def test_invalid_replication_request_rejected():
    with pytest.raises(ValueError):
        parallel_analysis(np.zeros((50,4)),np.eye(4),list('abcd'),replicates=0)
