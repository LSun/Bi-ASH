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

bilik = function (betahat0, sebetahat0, betahat1, sebetahat1, sigmaalpha, sigmagamma) {
  N1 = length(betahat0)
  N2 = length(betahat1)
  K = length(sigmaalpha)
  L = length(sigmagamma)
  J = N1 + N2
  array_lik = array(0, dim = c(K,L,J))
  for (j in 1:N1){
    for (k in 1:K){
      for (l in 1:L){
        array_lik[k,l,j] = dnorm(betahat0[j],0,sqrt(sebetahat0[j]^2+sigmaalpha[k]^2))
      }
    }
  }
  for (j in (N1+1):(N1+N2)){
    for (k in 1:K){
      for (l in 1:L){
        array_lik[k,l,j] = dnorm(betahat1[j-N1],0,sqrt(sebetahat1[j-N1]^2+sigmaalpha[k]^2+sigmagamma[l]^2))
      }
    }
  }
  return(array_lik)
}
