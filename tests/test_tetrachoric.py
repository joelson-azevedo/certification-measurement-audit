"""Numerical reference cases independent of the delivered synthetic examination."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr

from src.tetrachoric import TetrachoricConfig, bivariate_normal_cdf, estimate_tables, tetrachoric_matrix


@pytest.mark.parametrize("table,expected", [
    ([100,100,100,100],0), ([200,100,100,200],.5),
    ([375,125,125,375],np.sqrt(.5)), ([125,375,375,125],-np.sqrt(.5)),
    ([60,140,90,210],0),
])
def test_tetrachoric_known_numerical_references(table,expected):
    # Balanced margins have rho=sin(2*pi*(p00-1/4)); asymmetric example is independence.
    result=estimate_tables([table]).iloc[0]
    assert result.converged
    assert result.tetrachoric_r==pytest.approx(expected,abs=2e-8)
    assert result.absolute_probability_error<1e-9


@pytest.mark.parametrize("a,b,rho",[(0,0,.5),(-1.5,.4,.7),(.8,-.2,-.4),(-2,-1.2,.8)])
def test_plackett_cdf_against_independent_conditional_normal_integration(a,b,rho):
    # Different integral representation, no production CDF implementation reused.
    expected=quad(lambda z: np.exp(-z*z/2)/np.sqrt(2*np.pi)*ndtr((b-rho*z)/np.sqrt(1-rho*rho)),
                  -np.inf,a,epsabs=1e-11)[0]
    assert float(bivariate_normal_cdf(a,b,rho))==pytest.approx(expected,abs=1e-9)


def test_sparse_zero_cell_correction_is_explicit_and_preserves_raw_table():
    table=np.array([[40,0,10,50]])
    before=table.copy()
    row=estimate_tables(table).iloc[0]
    assert row.converged and row.continuity_corrected and row.sparse_table
    assert row.correction_added_to_each_cell==.5
    assert [row.n00,row.n01,row.n10,row.n11]==[40,0,10,50]
    np.testing.assert_equal(table,before)


def test_constant_margin_is_undefined_not_artificially_rescued():
    row=estimate_tables([[20,30,0,0]]).iloc[0]
    assert not row.converged and np.isnan(row.tetrachoric_r)
    assert not row.continuity_corrected and row.status=="constant binary margin"


def test_failure_reporting_and_insufficient_pair_n():
    row=estimate_tables([[40,10,20,30]],TetrachoricConfig(max_iterations=0)).iloc[0]
    assert not row.converged and np.isnan(row.tetrachoric_r) and "nonconvergence" in row.status
    row=estimate_tables([[1,1,1,1]]).iloc[0]
    assert not row.converged and row.status=="insufficient pair N"


def test_matrix_symmetry_diagonal_order_and_pairwise_denominators():
    x=np.array([[0,0,1],[0,1,np.nan],[1,0,0],[1,1,1]],float)
    names=["z_item","a_item","shared"]
    matrix,pairs,n=tetrachoric_matrix(x,names,TetrachoricConfig(minimum_pair_n=2))
    assert pairs.converged.all()
    np.testing.assert_allclose(matrix,matrix.T)
    np.testing.assert_equal(np.diag(matrix),np.ones(3))
    assert pairs.item_i.tolist()==["z_item","z_item","a_item"]
    assert pairs.item_j.tolist()==["a_item","shared","shared"]
    assert n[0,1]==4 and n[0,2]==3 and n[2,2]==3
    assert pairs.n_pair.tolist()==[4,3,3]


@pytest.mark.parametrize("x", [[[0,2],[1,0]],[0,1],[[0,np.inf],[1,0]]])
def test_invalid_binary_data_rejected(x):
    with pytest.raises(ValueError):
        tetrachoric_matrix(x)


def test_duplicate_item_order_rejected():
    with pytest.raises(ValueError):
        tetrachoric_matrix([[0,0],[1,1]],["same","same"])
