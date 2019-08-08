#' @title Fit Empirical Bayes Normal Means with Correlated Noise
#'
#' @description This is the main interface for fitting correlated EBNM models
#'   based on algorithms proposed by Sun and Stephens.  The default
#'   behaviour is simply to run the biconvex optimization and return
#'   the result and posterior calculations.
#'
#' @param x A p vector of observations
#' @param s A scalar or a p vector of standard deviations.
#' @param deltaAt0 Logical, indicating whether to use a point mass at zero as one of components for a mixture distribution of the prior.
#' @param gd.order
#' @param omega.lambda
#' @param omega.rho
#' @param omega.pen
#' @param mixsd.mult
#' @param gd.priority Logical, indicating whether to optimizer over prior
#' @param control
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom utils capture.output modifyList
#'
biashr = function (x1, s1 = 1,
                 x2, s2 = 1,
                 pi.at.0 = TRUE,
                 pi.mixsd.mult = sqrt(2),
                 pi.null.weight = 10,
                 omega.at.0 = TRUE,
                 omega.mixsd.mult = sqrt(2),
                 omega.null.weight = 10,
                 pi.first = FALSE,
                 control = list(maxiter = 50)) {

  if (all(s1 > 0) & length(s1) == 1L) {
    s1 = rep(s1, length(x1))
  } else if (length(x1) != length(s1) | !all(s1 > 0)) {
    stop("s1 should either be a positive number or a vector of positive numbers with the same length of x")
  }
  
  if (all(s2 > 0) & length(s2) == 1L) {
    s2 = rep(s2, length(x2))
  } else if (length(x2) != length(s2) | !all(s2 > 0)) {
    stop("s2 should either be a positive number or a vector of positive numbers with the same length of x")
  }

  ## setting a dense grid of sd for pi and omega
  if (pi.at.0) {
    sd1 = c(0, autoselect.mixsd(x1, s1, mult = pi.mixsd.mult))
  } else {
    sd1 = autoselect.mixsd(x1, s1, mult = pi.mixsd.mult)
  }
  
  if (omega.at.0) {
    sd2 = c(0, autoselect.mixsd(x2, s2, mult = omega.mixsd.mult))
  } else {
    sd2 = autoselect.mixsd(x2, s2, mult = omega.mixsd.mult)
  }
  
  sd.array <- sqrt(outer(outer(s1^2, sd1^2, FUN = "+"), sd2^2, FUN = "+"))
  mlik.array <- dnorm(x1, 0, sd.array)
  mlik.mat <- dnorm(x2, 0, sqrt(outer(s2^2, sd2^2, FUN = "+")))
  K <- dim(mlik.array)[2]
  L <- dim(mlik.array)[3]

  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                         tol = 1e-07, maxiter = 1000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)

  res = SQUAREM::squarem(par = c(1, rep(0, K - 1), 1, rep(0, L - 1)),
                         fixptfn = bifixpoint, objfn = negpenloglik,
                         mlik.array = mlik.array, mlik.mat = mlik.mat,
                         pi.at.0 = pi.at.0, pi.null.weight = pi.null.weight,
                         omega.at.0, omega.null.weight = omega.null.weight,
                         pi.first = pi.first,
                         control = controlinput)

  pi.hat = normalize(res$par[1 : K])
  omega.hat = normalize(res$par[-(1 : K)])
  g.fitted = normalmix(pi = pi.hat, mean = 0, sd = sd1)
  f.fitted = normalmix(pi = omega.hat, mean = 0, sd = sd2)
  penloglik = -res$value.objfn
  niter = res$fpevals
  converged = res$convergence
  
  theta.post.mean.array <- aperm(aperm(x1 / sd.array^2, c(2, 1, 3)) * sd1^2, c(2, 1, 3))
  theta.post.weight.array <- mlik.array / colSums(apply(aperm(mlik.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  theta.post.mean <- colSums(apply(aperm(theta.post.mean.array * theta.post.weight.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  
  if (pi.at.0) {
    theta.lfdr <- pi.hat[1] * colSums(t(dnorm(x1, 0, sqrt(outer(s1^2, sd2^2, FUN = "+")))) * omega.hat)/ colSums(apply(aperm(mlik.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  } else {
    theta.lfdr <- rep(0, length(x1))
  }

  theta.qvalue = qval.from.lfdr(theta.lfdr)
  
  output <- list(g.fitted = g.fitted,
                 f.fitted = f.fitted,
                 penloglik = penloglik,
                 niter = res$niter,
                 converged = res$converged,
                 theta.postmean = theta.post.mean,
                 theta.lfdr = theta.lfdr,
                 theta.qvalue = theta.qvalue
  )
  class(output) <- "biashr"
  
  return(output)
}

autoselect.mixsd = function (x, s, mult) {
  s = s[s != 0]
  sigmaamin = min(s)/10
  if (all(x^2 <= s^2)) {
    sigmaamax = 8 * sigmaamin
  } else {
    sigmaamax = 2 * sqrt(max(x^2 - s^2))
  }
  if (mult == 0) {
    return(c(0, sigmaamax/2))
  } else {
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint) : 0) * sigmaamax)
  }
}

bifixpoint <- function (pi.omega.current,
                        mlik.array, mlik.mat,
                        pi.at.0, pi.null.weight,
                        omega.at.0, omega.null.weight,
                        pi.first, control) {
  K <- dim(mlik.array)[2]
  L <- dim(mlik.array)[3]
  m <- dim(mlik.array)[1]
  n <- dim(mlik.mat)[1]
  pi.current <- pi.omega.current[1 : K]
  omega.current <- pi.omega.current[-(1 : K)]
  if (pi.first) {
    mlik.mat.pi <- apply(aperm(mlik.array, c(3, 2, 1)) * omega.current, 2, colSums)
    if (pi.at.0 & pi.null.weight != 0) {
      mlik.mat.pi.ext <- rbind(mlik.mat.pi, c(1, rep(0, K - 1)))
      optim.weight <- c(rep(1, m), pi.null.weight)
    } else {
      mlik.mat.pi.ext <- mlik.mat.pi
      optim.weight <- rep(1, m)
    }
    optim.pi <- mixsqp::mixsqp(L = mlik.mat.pi.ext,
                               w = optim.weight,
                               x0 = pi.current,
                               control = list(verbose = FALSE))
    pi.new <- optim.pi$x
    mlik.mat.omega <- rbind(apply(aperm(mlik.array, c(2, 3, 1)) * pi.new, 2, colSums), mlik.mat)
    if (omega.at.0 & omega.null.weight != 0) {
      mlik.mat.omega.ext <- rbind(mlik.mat.omega, c(1, rep(0, L - 1)))
      optim.weight <- c(rep(1, m + n), omega.null.weight)
    } else {
      mlik.mat.omega.ext <- mlik.mat.omega
      optim.weight <- rep(1, m + n)
    }
    optim.omega <- mixsqp::mixsqp(L = mlik.mat.omega.ext,
                                  w = optim.weight,
                                  x0 = omega.current,
                                  control = list(verbose = FALSE))
    omega.new <- optim.omega$x
  } else {
    mlik.mat.omega <- rbind(apply(aperm(mlik.array, c(2, 3, 1)) * pi.current, 2, colSums), mlik.mat)
    if (omega.at.0 & omega.null.weight != 0) {
      mlik.mat.omega.ext <- rbind(mlik.mat.omega, c(1, rep(0, L - 1)))
      optim.weight <- c(rep(1, m + n), omega.null.weight)
    } else {
      mlik.mat.omega.ext <- mlik.mat.omega
      optim.weight <- rep(1, m + n)
    }
    optim.omega <- mixsqp::mixsqp(L = mlik.mat.omega.ext,
                                  w = optim.weight,
                                  x0 = omega.current,
                                  control = list(verbose = FALSE))
    omega.new <- optim.omega$x
    mlik.mat.pi <- apply(aperm(mlik.array, c(3, 2, 1)) * omega.new, 2, colSums)
    if (pi.at.0 & pi.null.weight != 0) {
      mlik.mat.pi.ext <- rbind(mlik.mat.pi, c(1, rep(0, K - 1)))
      optim.weight <- c(rep(1, m), pi.null.weight)
    } else {
      mlik.mat.pi.ext <- mlik.mat.pi
      optim.weight <- rep(1, m)
    }
    optim.pi <- mixsqp::mixsqp(L = mlik.mat.pi.ext,
                               w = optim.weight,
                               x0 = pi.current,
                               control = list(verbose = FALSE))
    pi.new <- optim.pi$x
  }
  return(c(pi.new, omega.new))
}

negpenloglik = function (pi.omega.current,
                         mlik.array, mlik.mat,
                         pi.at.0, pi.null.weight,
                         omega.at.0, omega.null.weight,
                         pi.first, control) {
  K = dim(mlik.array)[2]
  pi.current = pi.omega.current[1 : K]
  omega.current = pi.omega.current[-(1 : K)]
  loglik = sum(log(pmax(0, colSums(t(apply(aperm(mlik.array, c(2, 3, 1)) * pi.current, 2, colSums)) * omega.current)))) +
    sum(log(pmax(0, colSums(t(mlik.mat) * omega.current))))
  if (pi.at.0 & pi.null.weight != 0) {
    penloglik = loglik + pi.null.weight * log(pi.current[1])
  } else {
    penloglik = loglik
  }
  if (omega.at.0 & omega.null.weight != 0) {
    penloglik = penloglik + omega.null.weight * log(omega.current[1])
  } else {
    penloglik = penloglik
  }
  return(-penloglik)
}

normalize = function (x) {
  return(x/sum(x))
}

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

qval.from.lfdr <- function (lfdr) {
  if (sum(!is.na(lfdr)) == 0) {
    return(rep(NA, length(lfdr)))
  }
  o = order(lfdr)
  qvalue = rep(NA, length(lfdr))
  qvalue[o] = (cumsum(sort(lfdr))/(1 : sum(!is.na(lfdr))))
  return(qvalue)
}





















set_control_squarem=function(control,nobs){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if (nobs > 50000) control.default$trace = TRUE
  control.default$tol = min(0.1/nobs,1.e-7) # set default convergence criteria to be more stringent for larger samples
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

set_control_mixIP=function(control){
  control.default=list(rtol=1e-6)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}




#' Title
#'
#' @param output
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.cash <- function (output, ...) {
  output[1 : 6]
}

#' Title
#'
#' @param output
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.cash <- function (output, ...) {
  print(summary.cash(output, ...))
}


#' Title
#'
#' @param output
#'
#' @return
#' @export
#'
#' @examples
get_svalue <- function (output) {
  array_PP <- array_PosProb(output$x, output$s, output$deltaAt0, output$fitted_g$sd, gd.ord = output$gd.order, gd.normalized = TRUE)
  array_PP = aperm(array_PP, c(2, 3, 1))
  theta_PosProb <- colSums(t(apply(output$fitted_g$pi * array_PP, 2, colSums)) * output$omega) / colSums(t(apply(output$fitted_g$pi * output$array_F, 2, colSums)) * output$omega)
  lfsr <- compute_lfsr(1 - output$lfdr - theta_PosProb, output$lfdr)
  svalue <- qval.from.lfdr(lfsr)
  return(list(lfsr = lfsr,
              svalue = svalue))
}

compute_lfsr <- function (NegativeProb, ZeroProb) {
  ifelse(NegativeProb > 0.5 * (1 - ZeroProb), 1 - NegativeProb, NegativeProb + ZeroProb)
}


# w.cvxr.uncns = function (matrix_lik_w, w.init = NULL) {
#   FF <- matrix_lik_w[, -1]
#   f <- matrix_lik_w[, 1]
#   p <- ncol(FF)
#   w <- CVXR::Variable(p)
#   objective <- CVXR::Maximize(CVXR::SumEntries(CVXR::Log(FF %*% w + f)))
#   prob <- CVXR::Problem(objective)
#   if (is.null(w.init)) {
#     capture.output(result <- solve(prob), file = "/dev/null")
#   } else {
#     capture.output(result <- solve(prob, warm_start = w.init[-1]), file = "/dev/null")
#   }
#   return(result)
# }

w.mosek = function (matrix_lik_w, matrix_lik_z, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[, 1]
  B = matrix_lik_z[, -1]
  b = matrix_lik_z[, 1]
  m = ncol(A)
  nA = nrow(A)
  nB = nrow(B)
  AB <- rbind(A, B)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    g.init <- as.vector(matrix_lik_w %*% w.init)
    v.init <- c(1 / g.init, rep(0, nB))
    v.init.list <- list(xx = v.init)
    P$sol <- list(itr = v.init.list, bas = v.init.list)
  }
  P$c <- c(a, b)
  P$A <- Matrix::Matrix(t(AB), sparse = TRUE)
  if (is.null(w_prior) | all(w_prior == 0) | missing(w_prior)) {
    P$bc <- rbind(rep(0, m), rep(0, m))
  } else {
    P$bc <- rbind(-w_prior, w_prior)
  }
  P$bx <- rbind(rep(0, nA + nB), rep(Inf, nA + nB))
  opro <- matrix(list(), nrow = 5, ncol = nA)
  rownames(opro) <- c("type", "j", "f", "g", "h")
  opro[1, ] <- as.list(rep("log", nA))
  opro[2, ] <- as.list(1 : nA)
  opro[3, ] <- as.list(rep(-1, nA))
  opro[4, ] <- as.list(rep(1, nA))
  opro[5, ] <- as.list(rep(0, nA))
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$suc - z$sol$itr$slc
  list(w = w, status = status)
}

w.mosek.primal = function (matrix_lik_w, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[,1]
  m = ncol(A)
  n = nrow(A)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    w.init.list <- list(xx = w.init)
    P$sol <- list(bas = w.init.list)
  }
  P$c <- rep(0, n + m)
  P$A <- Matrix::Matrix(cbind(diag(n), -A), sparse = TRUE)
  P$bc <- rbind(a, a)
  P$bx <- rbind(c(rep(0, n), rep(-Inf, m)),
                c(rep(Inf, n), rep(Inf, m)))
  if (missing(w_prior) | is.null(w_prior) | all(w_prior == 0)) {
    opro <- matrix(list(), nrow = 5, ncol = n)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(rep("log", n))
    opro[2, ] <- as.list(1 : n)
    opro[3, ] <- as.list(rep(-1, n))
    opro[4, ] <- as.list(rep(1, n))
    opro[5, ] <- as.list(rep(0, n))
  } else {
    opro <- matrix(list(), nrow = 5, ncol = n + m)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(c(rep("log", n), rep("pow", m)))
    opro[2, ] <- as.list(1 : (n + m))
    opro[3, ] <- as.list(c(rep(-1, n), w_prior))
    opro[4, ] <- as.list(c(rep(1, n), rep(2, m)))
    opro[5, ] <- as.list(rep(0, n + m))
  }
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$xx[-(1 : n)]
  list(w = w, status = status)
}


mixIP = function (matrix_lik, prior, pi_init = NULL, control = list()) {
  if(!requireNamespace("REBayes", quietly = TRUE)) {
    stop("mixIP requires installation of package REBayes")}
  control = set_control_mixIP(control)
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,rep(1,n))
  A = A[w!=0,] #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  #w = rep(1,n+k)
  res = REBayes::KWDual(A, rep(1,k), normalize(w), control=control)
  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL"), control=control))
}




## this function is vectorized for x
## more efficient if let it run for x at once
gauss.deriv = function(x, ord) {
  return(dnorm(x) * EQL::hermite(x, ord))
}

Hermite = function (gd.ord) {
  x <- PolynomF::polynom()
  H <- PolynomF::polylist(x, - 1 + x^2)
  if (gd.ord >= 3) {
    for(n in 2 : (gd.ord - 1))
      H[[n+1]] <- x * H[[n]] - n * H[[n-1]]
  }
  return(H)
}

array_pm = function (betahat, sebetahat, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  temp2_0 = dnorm(beta.std.mat)
  hermite = Hermite(gd.ord + 1)
  if (gd.normalized) {
    for (i in 0 : gd.ord) {
      temp2[, , i + 1] = temp2_0 * hermite[[i + 1]](beta.std.mat) * (-1)^(i + 1) / sqrt(factorial(i))
    }
  } else {
    for (i in 1 : (gd.ord + 1)) {
      temp2[, , i + 1] = temp2_0 * hermite[[i + 1]](beta.std.mat) * (-1)^(i + 1)
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  se.std.mat = sebetahat / sd.mat
  sd.std.mat2 = t(sd / t(sd.mat))^2
  temp1 = exp(outer(log(se.std.mat), 0 : gd.ord, FUN = "*"))
  for (ii in 0 : gd.ord) {
    temp1[, , ii + 1] = temp1[, , ii + 1] * sd.std.mat2
  }
  array_pm = (-1) * temp1 * temp2
  rm(temp1)
  rm(temp2)
  return(array_pm)
}

lfdr_top = function (pi0, w, betahat, sebetahat, gd.ord) {
  hermite = Hermite(gd.ord)
  gd.mat = matrix(0, ncol = gd.ord + 1, nrow = length(betahat))
  gd.mat[, 1] = dnorm(betahat / sebetahat)
  for (l in 1 : gd.ord) {
    gd.mat[, l + 1] = (-1)^l * hermite[[l]](betahat / sebetahat) * gd.mat[, 1] / sqrt(factorial(l))
  }
  gd.std.mat = gd.mat / sebetahat
  return(gd.std.mat %*% w * pi0)
}

array_PosProb = function (betahat, sebetahat, deltaAt0, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  se.std.mat <- sebetahat / sd.mat
  sd.std.mat <- t(sd / t(sd.mat))
  sd.se.mat <- 1 / outer(sebetahat, sd, FUN = "/")
  beta.std.sd.se.mat <- beta.std.mat * sd.se.mat
  pdf.beta.std.mat <- dnorm(beta.std.mat)
  pdf.beta.std.sd.se.mat <- dnorm(beta.std.sd.se.mat)
  cdf.beta.std.sd.se.mat <- pnorm(beta.std.sd.se.mat)
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  hermite = Hermite(gd.ord)
  temp2[, , 1] <- cdf.beta.std.sd.se.mat * pdf.beta.std.mat / sd.mat
  temp2[, , 2] <- cdf.beta.std.sd.se.mat * hermite[[1]](beta.std.mat) * (-1) * pdf.beta.std.mat * se.std.mat^2 / sebetahat +
    sd.se.mat * pdf.beta.std.sd.se.mat * pdf.beta.std.mat * se.std.mat^2 / sebetahat
  if (gd.normalized) {
    temp2[, , 3] <- (cdf.beta.std.sd.se.mat * hermite[[2]](beta.std.mat) * (-1)^2 * pdf.beta.std.mat * se.std.mat^3 / sebetahat +
                       2 * sd.se.mat * se.std.mat^3 / sebetahat * pdf.beta.std.sd.se.mat * hermite[[1]](beta.std.mat) * (-1) * pdf.beta.std.mat +
                       sd.se.mat^2 * se.std.mat^3 / sebetahat * hermite[[1]](beta.std.sd.se.mat) * (-1) * pdf.beta.std.sd.se.mat * pdf.beta.std.mat) / sqrt(factorial(2))
    for (i in 3 : gd.ord) {
      temp2[, , (i + 1)] <- cdf.beta.std.sd.se.mat * hermite[[i]](beta.std.mat) / sqrt(factorial(i)) * (-1)^i * pdf.beta.std.mat * se.std.mat^(i + 1) / sebetahat +
        i * sd.se.mat * se.std.mat^(i + 1) / sebetahat * pdf.beta.std.sd.se.mat * hermite[[i - 1]](beta.std.mat) / sqrt(factorial(i)) * (-1)^(i - 1) * pdf.beta.std.mat +
        sd.se.mat^i * se.std.mat^(i + 1) / sebetahat * hermite[[i - 1]](beta.std.sd.se.mat) / sqrt(factorial(i)) * (-1)^(i - 1) * pdf.beta.std.sd.se.mat * pdf.beta.std.mat +
        Reduce('+', lapply(2 : (i - 1), function (m) {
          sqrt(choose(i, m)) * hermite[[m - 1]](beta.std.sd.se.mat) / sqrt(factorial(m)) * pdf.beta.std.sd.se.mat * sd.std.mat^m *
            hermite[[i - m]](beta.std.mat) / sqrt(factorial(i - m)) * pdf.beta.std.mat * se.std.mat^(i - m) / sd.mat * (-1)^(i - 1)
        }))
    }
  } else {
    temp2[, , 3] <- temp3_1 * hermite[[2]](beta.std.mat) * (-1)^2 * temp2_0 +
      2 * temp3 * temp3_0 * hermite[[1]](beta.std.mat) * (-1) * temp2_0 +
      temp3^2 * hermite[[1]](beta.std.mat * temp3) * (-1) * temp3_0 * temp2_0
    for (i in 3 : gd.ord) {
      temp2[, , i + 1] <- ((temp3_1 * hermite[[i]](beta.std.mat) * (-1)^i +
                              i * temp3 * temp3_0 * hermite[[i - 1]](beta.std.mat) * (-1)^(i - 1) +
                              temp3^i * hermite[[i - 1]](beta.std.mat * temp3) * temp3_0) * temp2_0 +
                             sum(sapply(2 : (i - 1), function (m) {
                               exp(lchoose(i, m) + m * temp3) * hermite[[m - 1]](beta.std.mat * temp3) *
                                 hermite[[i - m]](beta.std.mat)
                             })) * temp3_0 * temp2_0) / sqrt(factorial(i))
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  if (deltaAt0) {
    array_PosProb = temp2
    array_PosProb[, 1, ] <- 0
  } else {
    array_PosProb = temp2
  }
  rm(temp2)
  return(array_PosProb)
}


plot.ghat = function (fitted.g, mixcompdist = "normal", xlim = c(-10, 10)) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  x = seq(xlim[1], xlim[2], 0.01)
  if (mixcompdist == "normal") {
    y = sum(pi * pnorm(x, mean, sd))
  }
}

ghat.cdf = function (x, fitted.g) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  return(sum(pi * pnorm(x, mean, sd)))
}

