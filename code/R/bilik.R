#' @title bilik
#' @description compute the density of the components of the normal mixture m when, in group0, convoluted with a normal with standard deviation s, whereas in group1, convoluted with another normal mixture and a normal with standard deviation s, the density is evaluated at x
#' @param betahat0 the observed \hat\beta in group0
#' @param sebetahat0 the observed standard deviation \hat\s in group0
#' @param betahat1 the observed \hat\beta in group1
#' @param sebetahat1 the observed standard deviation \hat\s in group1
#' @param sigmaalpha a grid of \sigma for the mixture model of the effect in both groups
#' @param sigmagamma a grid of \sigma for the mixture model of theeffect only in group1
#' @return a n by k by l array of densities

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
