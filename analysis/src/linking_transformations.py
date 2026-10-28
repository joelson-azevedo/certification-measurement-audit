"""Affine B-to-A transformations using the signed logistic intercept-slope form."""
from __future__ import annotations

import numpy as np
from scipy.special import expit

LABEL="PROVISIONAL LINKED SCALE — NOT OPERATIONAL CERTIFICATION SCORES"


def validate_constants(scale,shift):
    if not np.isfinite(scale) or not np.isfinite(shift) or scale<=0:
        raise ValueError("A finite positive scale and finite shift are required")


def transform_parameters(intercept,slope,scale,shift):
    validate_constants(scale,shift)
    d,a=np.asarray(intercept,float),np.asarray(slope,float)
    if d.shape!=a.shape or not np.isfinite(d).all() or not np.isfinite(a).all():
        raise ValueError("Finite matching intercept/slope arrays required")
    transformed_slope=a/scale
    return d-transformed_slope*shift,transformed_slope


def inverse_constants(scale,shift):
    validate_constants(scale,shift)
    return 1/scale,-shift/scale


def transform_posterior(mean,sd,scale,shift):
    validate_constants(scale,shift)
    mean,sd=np.asarray(mean,float),np.asarray(sd,float)
    if mean.shape!=sd.shape or not np.isfinite(mean).all() or not np.isfinite(sd).all() or (sd<0).any():
        raise ValueError("Finite posterior mean and nonnegative SD arrays required")
    return scale*mean+shift,scale*sd


def stable_difficulty(d,a,slope_se=None,minimum_absolute_slope=.10):
    d,a=np.asarray(d,float),np.asarray(a,float)
    valid=np.abs(a)>=minimum_absolute_slope
    if slope_se is not None:
        se=np.asarray(slope_se,float)
        valid&=np.isfinite(se)&(np.abs(a)>1.959963984540054*se)
    return np.divide(-d,a,out=np.full_like(d,np.nan),where=valid)


def fixed_link_covariance(covariance,scale,shift):
    """Transform a joint [all d, all a] covariance CONDITIONAL on fixed constants.

    This does not propagate estimated linking-constant uncertainty/dependence.
    """
    validate_constants(scale,shift)
    covariance=np.asarray(covariance,float);n=covariance.shape[0]
    if covariance.ndim!=2 or covariance.shape!=(n,n) or n%2:
        raise ValueError("Even square [all intercepts, all slopes] covariance required")
    k=n//2
    jacobian=np.block([[np.eye(k),-shift/scale*np.eye(k)],[np.zeros((k,k)),np.eye(k)/scale]])
    return jacobian@covariance@jacobian.T


def probabilities(theta,d,a):return expit(np.asarray(theta,float)[...,None]*np.asarray(a,float)+np.asarray(d,float))


def information(theta,d,a):
    p=probabilities(theta,d,a)
    return np.sum(np.asarray(a,float)**2*p*(1-p),axis=-1)
