library(SQUAREM)
library(methods)
source("R/bi-ash.R")
source("R/bilik.R")

N0 = 100
N1 = 100

NITER=5
converged = rep(0,NITER)
LLRTres = rep(0,NITER)
iter_time = rep(0,NITER)


for (w in 1:NITER){
t<-proc.time()

I = runif(N0)
sebetahat0 = rep(1,N0)
betahat0 = (I <= 0.5) * 0 + (I > 0.5) * rnorm(N0,0,2) + rnorm(N0,0,1)
I1 = runif(N1)
I2 = runif(N1)
sebetahat1 = rep(1,N1)
betahat1 = (I1 <= 0.5) * 0 + (I1 > 0.5) * rnorm(N1,0,2) + (I2 <= 0.9) * 0 + (I2 > 0.9) * rnorm(N1,0,5) + rnorm(N1,0,1)

Array_lik = bilik(betahat0, sebetahat0, betahat1, sebetahat1, c(0,2), c(0,5))
Pi=c(0.5,0.5,0.9,0.1)
Prior=rep(1,4)

sigmaalpha = ashr:::autoselect.mixsd(betahat0,sebetahat0,sqrt(2))
sigmagamma = ashr:::autoselect.mixsd(betahat1,sebetahat1,sqrt(2))
array_lik = bilik(betahat0, sebetahat0, betahat1, sebetahat1, sigmaalpha, sigmagamma)
prior = rep(1,length(sigmaalpha)+length(sigmagamma))

bioutput=bimixEM(array_lik,prior)
LLRT=2*((-bioutput$B)-bipenloglik(Pi,Array_lik,Prior))

t<-proc.time()-t

converged[w]=bioutput$converged
LLRTres[w] = LLRT
iter_time[w] = as.vector(t)[3]
}

system("mkdir ../output/biASHSim")
save(converged,LLRTres,iter_time,file="../output/biASHSim/res.RData")