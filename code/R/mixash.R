library(SQUAREM)
library(gaussquad)
library(ashr)
#source('mix.R')
#source('added.R')

#' @title Function to compute the local false sign rate
#'
#' @param NegativeProb A vector of posterior probability that beta is negative.
#' @param ZeroProb A vector of posterior probability that beta is zero.
#' @return The local false sign rate.
compute_lfsr = function(NegativeProb,ZeroProb){
  ifelse(NegativeProb> 0.5*(1-ZeroProb),1-NegativeProb,NegativeProb+ZeroProb)
}  


# If x is a n-column vector, turn it into n by 1 matrix
# If x is a matrix, keep it
tomatrix = function(x){
  if(is.vector(x)){
    x = as.matrix(x)
  }
  return(x)
}

#estimate mixture proportions of beta's prior by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)
EMest_mean = function(betahat,sebetahat,pilik,g,prior,null.comp=1,nullcheck=FALSE, df, control=list()){   
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  pi.init = g$pi
  k = ncomp(g)
  n = length(betahat)
  l = dim(tomatrix(sebetahat))[2]
  group=rep(1:k,l)
  
  controlinput$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples
  
  if(controlinput$trace==TRUE){tic()}
  
  matrix_lik_raw = t(compdens_conv_mixlik(g,betahat,sebetahat,df,pilik))
  matrix_lik = t(rowsum(t(matrix_lik_raw),group))
  
  EMfit = mixEM(matrix_lik,prior,pi.init, control=controlinput)
  
  pi = EMfit$pihat     
  loglik = EMfit$B # actually return log lower bound not log-likelihood! 
  converged = EMfit$converged
  niter = EMfit$niter
  loglik.final = EMfit$B[length(EMfit$B)]
  
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  
  if(nullcheck==TRUE){ 
    if(null.loglik > loglik.final){ #check whether exceeded "null" likelihood where everything is null
      pi=rep(0,k)
      pi[null.comp]=1
      m  = t(pi * t(matrix_lik)) 
      m.rowsum = rowSums(m)
      loglik = sum(log(m.rowsum))
    }
  }
  
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}
  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g))
}

# Approximate non-standard mixture t-likelihood by normal-mixture
# component i: ~scale[i]*T(df[i]), w.p. pi[i]
approxlik_gq=function(params,q,appsigma,appweight){  
  pi=params[1:q]
  scale=params[(q+1):(2*q)]
  df=params[(2*q+1):(3*q)]
  
  fi=numeric(0)
  sigma=numeric(0)
  for (i in 1:q){
    fi = c(fi,pi[i]*appweight[i,])
    sigma = c(sigma, scale[i]*appsigma[i,])
  }
  fi=fi[order(sigma)]
  sigma=sort(sigma)
  return(c(fi,sigma))
}

# Approximate t-distribution (of df) by r-components normal mixture
approxt = function(df, r){
  alpha=df/2-1
  rules=glaguerre.quadrature.rules(r,alpha,normalized=TRUE)
  sigma=sqrt(df/(2*rules[[r]]$x))
  weight=rules[[r]]$w/sum(rules[[r]]$w)
  return(list(sigma=sigma,weight=weight))
}

# Approximate mixture t likelihood (with l components) 
# by mixture normal (with q components)
# pi, alpha, beta are n by l matrices
# component i: ~(scale[n,i])*T(df[n,i]), w.p. pi[n,i]
mixlik_sd=function(pi,scale,df){
  q=dim(tomatrix(pi))[2]
  ll=max(5,floor(20/q))
  params=cbind(pi,scale,df)
  
  appweight=matrix(rep(0,ll*q),nrow=q)
  appsigma=matrix(rep(0,ll*q),nrow=q)
  for (i in 1:q){
    app = approxt(df=df[i],r=ll) # l components for approximating each t-distribution
    appweight[i,]=app$weight
    appsigma[i,]=app$sigma
  }
  
  results = t(apply(params,1,approxlik_gq,q=q,appsigma=appsigma,appweight=appweight))
  pilik=results[,1:(dim(results)[2]/2)]
  selik=results[,(dim(results)[2]/2+1):dim(results)[2]]
  return(list(pilik=pilik,selik=selik))
}


#compute posterior shape (alpha1) and rate (beta1)
post.igmix = function(m,betahat,sebetahat,v){
  n = length(sebetahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2)
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing]=m$beta
  return(list(alpha=alpha1,beta=t(beta1)))
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd = function(betahat,sebetahat,mult){
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  } else {
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}



#' @title Main Adaptive SHrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat, a p vector of estimates 
#' @param sebetahat, a p vector of corresponding standard errors
#' @param method: specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist: distribution of components in mixture ("normal", "uniform" or "halfuniform")
#'
#' @param lambda1: multiplicative "inflation factor" for standard errors (like Genomic Control)
#' @param lambda2: additive "inflation factor" for standard errors (like Genomic Control)
#' @param nullcheck: whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component
#' @param df: appropriate degrees of freedom for (t) distribution of betahat/sebetahat
#' @param randomstart: bool, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param pointmass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR: bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, lfdr...
#' @param singlecomp: bool, indicating whether to use a single inverse-gamma distribution as the prior distribution of the variances
#' @param SGD: bool, indicating whether to use the stochastic gradient descent method to fit the prior distribution of the variances
#' @param unimodal: unimodal constraint for the prior distribution of the variances ("variance") or the precisions ("precision")
#' @param prior: string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param mixsd: vector of sds for underlying mixture components 
#' @param gridmult: the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
#' @param minimal_output: if TRUE, just outputs the fitted g and the lfsr (useful for very big data sets where memory is an issue) 
#' @param g: the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). User may supply changes to this list of parameter, say, control=list(maxiter=10000,trace=TRUE)
#' 
#'
#' @return a list with elements fitted.g is fitted mixture
#' logLR : logP(D|mle(pi)) - logP(D|null)
#' 
#' @export
#' 
#' @examples 
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' summary(beta.ash)
#' plot(betahat,beta.ash$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
mixash = function(betahat,sebetahat,df,pilik,
                  method = c("shrink","fdr"), 
                  mixcompdist = c("normal","uniform","halfuniform"),
                  lambda1=1,lambda2=0,nullcheck=FALSE,randomstart=FALSE, 
                  pointmass = TRUE, 
                  onlylogLR = FALSE, 
                  singlecomp = FALSE,
                  SGD = TRUE,
                  prior=c("uniform","nullbiased"), 
                  mixsd=NULL,gridmult=sqrt(2),
                  minimaloutput=FALSE,
                  g=NULL,
                  control=list()){
  
  #method provides a convenient interface to set a particular combinations of parameters for prior an
  #If method is supplied, use it to set up specific values for these parameters; provide warning if values
  #are also specified by user
  #If method is not supplied use the user-supplied values (or defaults if user does not specify them)
  if(length(pilik)==1){
    pilik=rep(1,length(betahat))
  }else if(dim(tomatrix(pilik))[1]!=length(betahat)){
    stop("Error: pilik must be 1, or in same shape as sebetahat.")
    
  }
  
  #   if(length(df)==1){
  #     df=sebetahat/sebetahat*df
  #   }else if(dim(tomatrix(sebetahat))[2]>1 & length(df)==dim(tomatrix(sebetahat))[2]){
  #     df=matrix(rep(df,each=dim(sebetahat)[1]),ncol=length(df))
  #   }else if(dim(tomatrix(sebetahat))[1]>1 & length(df)==dim(tomatrix(sebetahat))[1]){
  #     df=matrix(rep(df,dim(sebetahat)[2]),nrow=length(df))
  #   }else{
  #     stop("Error: df must have length 1, or same length as betahat, or same as dim(sebetahat)[2].")
  #   }
  
  
  if(!missing(method)){
    method = match.arg(method) 
    if(method=="shrink"){
      if(missing(prior)){
        prior = "uniform"
      } else {
        warning("Specification of prior overrides default for method shrink")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        warning("Specification of pointmass overrides default for method shrink")
      }
    }
    
    if(method=="fdr"){
      if(missing(prior)){
        prior = "nullbiased"
      } else {
        warning("Specification of prior overrides default for method fdr")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        warning("Specification of pointmass overrides default for method fdr")
      }
    }  
  }
  
  if(onlylogLR){
    pointmass = TRUE  
  }
  
  mixcompdist = match.arg(mixcompdist)
  
  if(!is.numeric(prior)){
    prior = match.arg(prior)
  }  
  
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(length(pilik)==dim(tomatrix(sebetahat))[2]){
    pilik = t(rep(pilik,length(betahat)),ncol=length(betahat))
  }
  if(dim(tomatrix(sebetahat))[1] != length(betahat)){
    stop("Error: sebetahat must have length 1, or same length as betahat")
  }
  
  completeobs = (!is.na(betahat) & !is.na(apply(tomatrix(sebetahat),1,sum)) & 
                   !is.na(apply(tomatrix(pilik),1,sum)))
  n=sum(completeobs)
  
  if(n==0){
    if(onlylogLR){
      return(list(pi=NULL, logLR = 0))
    }
    else{
      stop("Error: all input values are missing")
    }
  }  
  
  pilik = tomatrix(pilik)
  sebetahat = tomatrix(sebetahat)
  if(mixcompdist=="normal"){
    appnorm = mixlik_sd(pilik[completeobs,],sebetahat[completeobs,],df)
    pilik = matrix(rep(NA,length(betahat)*dim(appnorm$pilik)[2]),ncol=dim(appnorm$pilik)[2])
    pilik[completeobs,] = appnorm$pilik
    selik = matrix(rep(NA,length(betahat)*dim(pilik)[2]),ncol=dim(pilik)[2])
    selik[completeobs,] = appnorm$selik
    moddf = NULL
  }else if(mixcompdist=="uniform" | mixcompdist=="halfuniform"){
    pilik = pilik
    selik = sebetahat
    moddf = df
  }
  selik = tomatrix(selik)
  
  l = dim(pilik)[2]
  
  #Handling control variables
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if(n>50000){control.default$trace=TRUE}
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  if(!is.null(g)){
    controlinput$maxiter = 1 # if g is specified, don't iterate the EM
    prior = rep(1,ncomp(g)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    null.comp=1 #null.comp also not used, but required 
  } else {
    if(is.null(mixsd)){
      mixsd = autoselect.mixsd(betahat[completeobs],apply(tomatrix(pilik[completeobs,]*selik[completeobs,]),1,sum),gridmult)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }
    
    null.comp = which.min(mixsd) #which component is the "null"
    
    k = length(mixsd)
    if(!is.numeric(prior)){
      if(prior=="nullbiased"){ # set up prior to favour "null"
        prior = rep(1,k)
        prior[null.comp] = 10 #prior 10-1 in favour of null
      }else if(prior=="uniform"){
        prior = rep(1,k)
      }
    }
    
    if(length(prior)!=k | !is.numeric(prior)){
      stop("invalid prior specification")
    }
    
    if(randomstart){
      pi = rgamma(k,1,1)
    } else {
      if(k<n){
        pi=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
        pi[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
      } else {
        pi=rep(1,k)/k
      }
    }
    
    pi=pi/sum(pi)
    if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist")
    if(mixcompdist=="normal") g=normalmix(pi,rep(0,k),mixsd)
    if(mixcompdist=="uniform") g=unimix(pi,-mixsd,mixsd)
    if(mixcompdist=="halfuniform"){
      g = unimix(c(pi,pi)/2,c(-mixsd,rep(0,k)),c(rep(0,k),mixsd))
      prior = rep(prior, 2)
      pi = rep(pi, 2)
    }
  }
  
  
  pi.fit=EMest_mean(betahat[completeobs],lambda1*selik[completeobs,]+lambda2,pilik[completeobs,],g,prior,null.comp=null.comp,
                    nullcheck=nullcheck,df=moddf[completeobs],control=controlinput)  
  
  if(onlylogLR){
    logLR = tail(pi.fit$loglik,1) - pi.fit$null.loglik
    return(list(fitted.g=pi.fit$g, logLR = logLR))
  } else if(minimaloutput){
    n=length(betahat)
    ZeroProb = rep(0,length=n)
    NegativeProb = rep(0,length=n)
    
    #print("normal likelihood")
    ZeroProb[completeobs] = colSums(comppostprob_mixlik(pi.fit$g,betahat[completeobs],sebetahat[completeobs,],moddf,pilik[completeobs,])[comp_sd(pi.fit$g)==0,,drop=FALSE])     
    NegativeProb[completeobs] = cdf_post_mixlik(pi.fit$g, 0, betahat[completeobs],sebetahat[completeobs,],moddf,pilik[completeobs,]) - ZeroProb[completeobs]
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0) 
    
    lfsr = compute_lfsr(NegativeProb,ZeroProb)
    result = list(fitted.g=pi.fit$g,lfsr=lfsr,fit=pi.fit)
    return(result) 
  } else{
    
    
    #     post = posterior_dist(pi.fit$g,betahat,sebetahat)
    n=length(betahat)
    ZeroProb = rep(0,length=n)
    NegativeProb = rep(0,length=n)
    PosteriorMean = rep(0,length=n)
    PosteriorSD=rep(0,length=n)
    
    pilikco = tomatrix(pilik[completeobs,])
    selikco = tomatrix(selik[completeobs,])
    
    ZeroProb[completeobs] = colSums(comppostprob_mixlik(pi.fit$g, betahat[completeobs], selikco, moddf,  pilikco)[comp_sd(pi.fit$g)==0,,drop=FALSE])    
    NegativeProb[completeobs] = cdf_post_mixlik(pi.fit$g, 0, betahat[completeobs],selikco,moddf, pilikco) - ZeroProb[completeobs]
    PosteriorMean[completeobs] = postmean_mixlik(pi.fit$g,betahat[completeobs],selikco,moddf, pilikco)
    PosteriorSD[completeobs] =postsd_mixlik(pi.fit$g,betahat[completeobs],selikco,moddf, pilikco) 
    
    #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0) 
    PosteriorMean[!completeobs] = calc_mixmean(pi.fit$g)
    PosteriorSD[!completeobs] =calc_mixsd(pi.fit$g)  
    PositiveProb =  1- NegativeProb-ZeroProb    
    
    lfsr = compute_lfsr(NegativeProb,ZeroProb)
    #lfsra =  compute_lfsra(PositiveProb,NegativeProb,ZeroProb) 
    
    lfdr = ZeroProb
    qvalue = qval.from.lfdr(lfdr)
    
    result = list(fitted.g=pi.fit$g,logLR =tail(pi.fit$loglik,1) - pi.fit$null.loglik,
                  PosteriorMean = PosteriorMean,PosteriorSD=PosteriorSD,
                  PositiveProb =PositiveProb,NegativeProb=NegativeProb, ZeroProb=ZeroProb,lfsr = lfsr,
                  #lfsra=lfsra, 
                  lfdr=lfdr,qvalue=qvalue,
                  fit=pi.fit,lambda1=lambda1,lambda2=lambda2,call=match.call(),data=list(betahat = betahat, sebetahat=sebetahat))
    class(result)= "mixash"
    return(result)
    
  }
  
}


