"""Signed-slope MML with analytic gradients and Louis observed information."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import time

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.special import expit, logsumexp, logit

from .irt_2pl import normal_quadrature, orient_parameters, validate_responses


@dataclass(frozen=True)
class EstimationConfig:
    primary_nodes: int = 31
    sensitivity_nodes: int = 41
    slope_bound: float = 6.0
    intercept_bound: float = 10.0
    max_iterations: int = 1500
    optimizer_gtol: float = 1e-8
    optimizer_ftol: float = 1e-14
    first_order_tolerance_per_candidate: float = 1e-6
    boundary_tolerance: float = 1e-4
    difficulty_minimum_abs_slope: float = .10
    information_minimum_eigenvalue: float = 1e-8
    information_maximum_condition_number: float = 1e10
    directional_hessian_tolerance: float = 1e-5

    def as_dict(self):
        return asdict(self)


class MarginalLikelihood:
    def __init__(self, x, nodes=31):
        self.x = validate_responses(x)
        self.n, self.k = self.x.shape
        if self.n < 2 or self.k < 3:
            raise ValueError("At least two candidates and three items needed for calibration")
        self.observed = np.isfinite(self.x).astype(float)
        self.correct = np.nan_to_num(self.x, nan=0.0)
        self.theta, self.weights = normal_quadrature(nodes)
        self.nodes = nodes
        self.log_weights = np.log(self.weights)

    def components(self, parameters):
        d, a = parameters[:self.k], parameters[self.k:]
        eta = self.theta[:,None]*a+d
        log_p, log_q = -np.logaddexp(0,-eta), -np.logaddexp(0,eta)
        conditional = self.correct@log_p.T+(self.observed-self.correct)@log_q.T
        joint = conditional+self.log_weights[None,:]
        log_marginal = logsumexp(joint,axis=1)
        posterior = np.exp(joint-log_marginal[:,None])
        return eta, posterior, log_marginal, conditional

    def value_gradient(self, parameters):
        eta, posterior, log_marginal, _ = self.components(parameters)
        probabilities = expit(eta)
        expected_correct = posterior.T@self.correct
        expected_observed = posterior.T@self.observed
        residual = expected_correct-expected_observed*probabilities
        gradient = -np.concatenate([residual.sum(axis=0), (residual*self.theta[:,None]).sum(axis=0)])
        return -float(log_marginal.sum()), gradient

    def average_value_gradient(self, parameters):
        value, gradient = self.value_gradient(parameters)
        return value/self.n, gradient/self.n

    def observed_information(self, parameters):
        """Louis identity: E[-H_complete|X] - Cov(score_complete|X), summed over people."""
        eta, posterior, _, _ = self.components(parameters)
        probabilities = expit(eta)
        residual = self.correct[:,None,:]-self.observed[:,None,:]*probabilities[None,:,:]
        scores = np.concatenate([residual,residual*self.theta[None,:,None]],axis=2)
        expected_score = np.einsum("nq,nqp->np",posterior,scores,optimize=True)
        weighted = (scores*np.sqrt(posterior)[:,:,None]).reshape(-1,2*self.k)
        covariance_sum = weighted.T@weighted-expected_score.T@expected_score
        expected_observed = posterior.T@self.observed
        weight = expected_observed*probabilities*(1-probabilities)
        d2 = weight.sum(axis=0)
        da = (weight*self.theta[:,None]).sum(axis=0)
        a2 = (weight*self.theta[:,None]**2).sum(axis=0)
        complete_information = np.block([[np.diag(d2),np.diag(da)],[np.diag(da),np.diag(a2)]])
        information = complete_information-covariance_sum
        return (information+information.T)/2


def initialization(x):
    x = validate_responses(x)
    observed = np.isfinite(x)
    correct = np.nan_to_num(x,nan=0.0)
    totals, counts = correct.sum(axis=1), observed.sum(axis=1)
    p = np.divide(correct.sum(axis=0), observed.sum(axis=0), out=np.full(x.shape[1],.5),where=observed.sum(axis=0)>0)
    slopes=[]
    for j in range(x.shape[1]):
        mask=observed[:,j]
        rest_count=counts[mask]-1
        rest=np.divide(totals[mask]-correct[mask,j],rest_count,out=np.zeros(mask.sum()),where=rest_count>0)
        xx=correct[mask,j]
        r=np.corrcoef(xx,rest)[0,1] if len(xx)>1 and np.var(xx)>0 and np.var(rest)>0 else 0.0
        slopes.append(np.clip(2*r/np.sqrt(max(1-r*r,1e-6)),-2.0,2.0))
    a=np.asarray(slopes)
    d=logit(np.clip(p,.001,.999))*np.sqrt(1+np.pi*a*a/8)
    return np.concatenate([d,a])


def detect_boundaries(parameters,k,config):
    limits=np.concatenate([np.repeat(config.intercept_bound,k),np.repeat(config.slope_bound,k)])
    return np.abs(parameters)>=limits-config.boundary_tolerance


def projected_gradient(parameters,gradient,k,config):
    limits=np.concatenate([np.repeat(config.intercept_bound,k),np.repeat(config.slope_bound,k)])
    projected=gradient.copy()
    projected[(parameters<=-limits+config.boundary_tolerance)&(gradient>0)]=0
    projected[(parameters>=limits-config.boundary_tolerance)&(gradient<0)]=0
    return projected


def parameter_uncertainty(likelihood,parameters,config,boundaries):
    information=likelihood.observed_information(parameters)
    eigen=np.linalg.eigvalsh(information)
    condition=float(np.linalg.cond(information))
    errors=[]
    for phase in [.0,.7,1.3]:
        direction=np.sin(np.arange(len(parameters))+phase+.2)
        direction/=np.linalg.norm(direction)
        step=1e-4
        numerical=(likelihood.value_gradient(parameters+step*direction)[1]-likelihood.value_gradient(parameters-step*direction)[1])/(2*step)
        analytic=information@direction
        errors.append(float(np.linalg.norm(numerical-analytic)/max(1,np.linalg.norm(numerical))))
    supported=bool(np.isfinite(information).all() and eigen[0]>config.information_minimum_eigenvalue
                   and condition<config.information_maximum_condition_number
                   and max(errors)<=config.directional_hessian_tolerance and not boundaries.any())
    covariance=np.full_like(information,np.nan)
    status="available: positive-definite observed information; model conditional"
    if supported:
        try:
            covariance=cho_solve(cho_factor(information,lower=True),np.eye(len(parameters)))
        except np.linalg.LinAlgError:
            supported=False
    if not supported:
        status="unavailable: boundary, singular/indefinite/ill-conditioned information or derivative check failure"
    return covariance,information,{"uncertainty_available":supported,"uncertainty_status":status,
                                  "information_minimum_eigenvalue":float(eigen[0]),"information_maximum_eigenvalue":float(eigen[-1]),
                                  "information_condition_number":condition,"maximum_directional_hessian_relative_error":max(errors)}


def fit_mml(x,config=None,nodes=None,warm_start=None,calculate_uncertainty=True):
    config=config or EstimationConfig()
    nodes=nodes or config.primary_nodes
    likelihood=MarginalLikelihood(x,nodes)
    k=likelihood.k
    if np.any(likelihood.observed.sum(axis=0)<2) or np.any(likelihood.correct.sum(axis=0)==0) or np.any((likelihood.observed-likelihood.correct).sum(axis=0)==0):
        raise ValueError("An item lacks observed response variation; calibration unsupported")
    first=initialization(x)
    second=first.copy();second[k:]*=.65
    starts=[first,second] if warm_start is None else [np.asarray(warm_start,float),first]
    bounds=[(-config.intercept_bound,config.intercept_bound)]*k+[(-config.slope_bound,config.slope_bound)]*k
    attempts=[];candidates=[];started=time.perf_counter()
    for number,start in enumerate(starts,1):
        start=np.asarray(start,float).copy()
        start[:k]=np.clip(start[:k],-config.intercept_bound*.95,config.intercept_bound*.95)
        start[k:]=np.clip(start[k:],-config.slope_bound*.95,config.slope_bound*.95)
        result=minimize(likelihood.average_value_gradient,start,method="L-BFGS-B",jac=True,bounds=bounds,
                        options={"maxiter":config.max_iterations,"gtol":config.optimizer_gtol,"ftol":config.optimizer_ftol,"maxls":50,"maxcor":20})
        d,a,flipped=orient_parameters(result.x[:k],result.x[k:])
        parameters=np.concatenate([d,a])
        value,gradient=likelihood.value_gradient(parameters)
        projected=projected_gradient(parameters,gradient,k,config)
        gradient_size=float(np.max(np.abs(projected))/likelihood.n)
        valid=bool(result.success and np.isfinite(value) and gradient_size<=config.first_order_tolerance_per_candidate)
        attempts.append({"start":number,"optimizer_success":bool(result.success),"first_order_verified":gradient_size<=config.first_order_tolerance_per_candidate,
                         "converged":valid,"message":str(result.message),"iterations":int(result.nit),"function_evaluations":int(result.nfev),
                         "marginal_log_likelihood":-value,"projected_gradient_max_abs_per_candidate":gradient_size,"global_orientation_flipped":flipped})
        candidates.append((valid,value,parameters,result,gradient))
    eligible=[c for c in candidates if c[0]] or candidates
    chosen=min(eligible,key=lambda c:c[1])
    converged,value,parameters,result,gradient=chosen
    boundaries=detect_boundaries(parameters,k,config)
    covariance=np.full((2*k,2*k),np.nan);information=np.full_like(covariance,np.nan)
    uncertainty={"uncertainty_available":False,"uncertainty_status":"not estimated for nonconverged model or disabled",
                 "information_minimum_eigenvalue":np.nan,"information_maximum_eigenvalue":np.nan,
                 "information_condition_number":np.nan,"maximum_directional_hessian_relative_error":np.nan}
    if converged and calculate_uncertainty:
        covariance,information,uncertainty=parameter_uncertainty(likelihood,parameters,config,boundaries)
    trustworthy=bool(converged and not boundaries.any() and (uncertainty["uncertainty_available"] if calculate_uncertainty else True))
    warnings=[]
    if not converged:warnings.append("optimizer/first-order convergence not established")
    if boundaries.any():warnings.append("parameter safeguards reached; point/uncertainty unstable")
    if calculate_uncertainty and not uncertainty["uncertainty_available"]:warnings.append(uncertainty["uncertainty_status"])
    valid_ll=[r["marginal_log_likelihood"] for r in attempts if r["converged"]]
    diagnostic={"n_candidates":likelihood.n,"n_items":k,"n_observed_responses":int(likelihood.observed.sum()),
                "n_missing_responses":int((1-likelihood.observed).sum()),"n_all_missing_candidates":int((likelihood.observed.sum(axis=1)==0).sum()),
                "quadrature_nodes":nodes,"optimizer":"L-BFGS-B; analytic gradient of average negative marginal LL",
                "converged":converged,"numerically_admissible":trustworthy,"optimizer_success":bool(result.success),
                "iterations":int(result.nit),"marginal_log_likelihood":-value,
                "maximum_absolute_gradient":float(np.max(np.abs(gradient))),
                "projected_gradient_max_abs_per_candidate":float(np.max(np.abs(projected_gradient(parameters,gradient,k,config)))/likelihood.n),
                "parameter_boundary_hits":int(boundaries.sum()),"negative_slopes":int((parameters[k:]<0).sum()),
                "near_zero_slopes":int((np.abs(parameters[k:])<config.difficulty_minimum_abs_slope).sum()),
                "orientation":"sum of item slopes positive; symmetric global orientation only",
                "converged_start_log_likelihood_range":float(max(valid_ll)-min(valid_ll)) if valid_ll else np.nan,
                "elapsed_seconds":time.perf_counter()-started,"warnings":"; ".join(warnings),**uncertainty}
    return {"intercept":parameters[:k],"slope":parameters[k:],"parameters":parameters,"covariance":covariance,
            "observed_information":information,"boundary_mask":boundaries,"diagnostics":diagnostic,"start_diagnostics":pd.DataFrame(attempts)}
