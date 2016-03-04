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
