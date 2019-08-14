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

  control.default = list(tol = 1e-07, maxiter = 1000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)

  res = SQUAREM::fpiter(par = c(1, rep(0, K - 1), 1, rep(0, L - 1)),
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
  
  theta.postweight.array <- mlik.array / colSums(apply(aperm(mlik.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  theta.postmean.array <- aperm(aperm(x1 / sd.array^2, c(2, 1, 3)) * sd1^2, c(2, 1, 3))
  theta.postsd.array <- sqrt(aperm(aperm(s1^2 / sd.array^2, c(2, 1, 3)) * sd1^2, c(2, 1, 3)) + aperm(aperm(aperm(1 / sd.array^2, c(2, 1, 3)) * sd1^2, c(3, 2, 1)) * sd2^2, c(2, 3, 1)))

  theta.postmean <- colSums(apply(aperm(theta.postmean.array * theta.postweight.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  
  if (pi.at.0) {
    theta.lfdr <- pi.hat[1] * colSums(t(dnorm(x1, 0, sqrt(outer(s1^2, sd2^2, FUN = "+")))) * omega.hat)/ colSums(apply(aperm(mlik.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  } else {
    theta.lfdr <- rep(0, length(x1))
  }
  theta.qvalue = qval.from.lfdr(theta.lfdr)
  
  theta.postposprob.array <- pnorm(0, theta.postmean.array, theta.postsd.array, lower.tail = FALSE)
  theta.postposprob <- colSums(apply(aperm(theta.postposprob.array * theta.postweight.array, c(2, 1, 3)) * pi.hat, 2, colSums) * omega.hat)
  theta.lfsr <- compute_lfsr(1 - theta.lfdr - theta.postposprob, theta.lfdr)
  theta.svalue <- qval.from.lfdr(theta.lfsr)
  
  output <- list(g.fitted = g.fitted,
                 f.fitted = f.fitted,
                 penloglik = penloglik,
                 converged = converged,
                 niter = niter,
                 theta.postmean = theta.postmean,
                 theta.lfdr = theta.lfdr,
                 theta.qvalue = theta.qvalue,
                 theta.lfsr = theta.lfsr,
                 theta.svalue = theta.svalue
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

compute_lfsr <- function (NegativeProb, ZeroProb) {
  ifelse(NegativeProb > 0.5 * (1 - ZeroProb), 1 - NegativeProb, NegativeProb + ZeroProb)
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
summary.biashr <- function (output, ...) {
  output[1 : 5]
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
print.biashr <- function (output, ...) {
  print(summary.biashr(output, ...))
}

compute_cdf = function (x, normalmix) {
  pi = normalmix$pi
  mean = normalmix$mean
  sd = normalmix$sd
  return(sum(pi * pnorm(x, mean, sd)))
}
