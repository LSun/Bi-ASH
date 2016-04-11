# This file contains methods for optimizing the mixture proportions given data
# from known mixture components.
# The penalized likelihood being maximized
# is \sum_j \log \sum_k\sum_l \pi_k\pi_l f_{jkl} + \sum_j (prior_j-1) \log \pi_k + 
# Adapted from Matthew Stephens's ashr package.



#' @title Estimate mixture proportions of a mixture model by EM algorithm
#'
#' @description Given the individual component likelihoods for a 3-way mixture model, estimates the mixture proportions by an EM algorithm.
#'
#' @details Fits a 3-way (k,l) component mixture model \deqn{f(x|\pi)= \sum_k\sum_l \pi_{k}\pi_l f_{kl}(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates mixture proportions \eqn{\pi_k} and \eqn{\pi_l} by maximum likelihood, or by maximum a posteriori (MAP) estimation for a Dirichlet prior on \eqn{\pi} 
#' (if a prior is specified).  Uses the SQUAREM package to accelerate convergence of EM. Used by the bimixEM main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' 
#' @param array_lik, a 3-way n by k by L matrix with (j,k,l)th element equal to \eqn{f_{kl}(x_j)}.
#' @param prior, a k+l vector of the parameters of the Dirichlet prior on \eqn{\pi_k} (first k) and \eqn{\pi_l} (last l). Recommended to be rep(1,k+l)
#' @param pi_init, the initial value of \eqn{\pi_k} and \eqn{\pi_l} to use. If not specified defaults to (1/k,...,1/k,1/l,...,1/l).
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 
#' 
#' @return A list, including the estimates (pihat)=c(\eqn{\pi_k},\eqn{\pi_l}), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
#' 

bimixEM = function (array_lik, prior, pi_init = NULL, control = list()) 
{
  control.default = list(K = 1, method = 3, square = TRUE, 
                         step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1, 
                         tol = 1e-07, maxiter = 5000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  K = dim(array_lik)[1]
  L = dim(array_lik)[2]
  if (is.null(pi_init)) {
    pi_init = c(rep(1/K, K),rep(1/L,L))
  }
  res = squarem(par = pi_init, fixptfn = bifixpoint, objfn = binegpenloglik, 
                array_lik = array_lik, prior = prior, control = controlinput)
  return(list(pihat = 2* normalize(pmax(0, res$par)), B = res$value.objfn, 
              niter = res$iter, converged = res$convergence))
}

# helper functions used by mixEM

bifixpoint = function(pi, array_lik, prior){  
  K = dim(array_lik)[1]
  L = dim(array_lik)[2]
  pi_k = pi[1:K]
  pi_l = pi[-(1:K)]
  prior_k = prior[1:K]
  prior_l = prior[-(1:K)]
  matrix_lik = apply(aperm(array_lik,c(2,1,3))*pi_l,2,colSums)
  pi_knew = fixpoint(pi_k,matrix_lik,prior_k)
  matrix_lik = apply(array_lik*pi_knew,2,colSums)
  pi_lnew = fixpoint(pi_l,matrix_lik,prior_l)
  return(c(pi_knew,pi_lnew))
}

fixpoint = function (pi, matrix_lik, prior) 
{
  pi = normalize(pmax(0, pi))
  m = t(pi * t(matrix_lik))
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}

normalize = function (x) 
{
  return(x/sum(x))
}

binegpenloglik = function (pi, array_lik, prior) 
{
  return(-bipenloglik(pi, array_lik, prior))
}

bipenloglik = function (pi, array_lik, prior) 
{
  K = dim(array_lik)[1]
  pi_k = pi[1:K]
  pi_l = pi[-(1:K)]
  loglik = sum(log(colSums(t(apply(pi_k*array_lik,2,colSums))*pi_l)))
  subset = (prior != 1)
  priordens = sum((prior - 1)[subset] * log(pi[subset]))
  return(loglik + priordens)
}
