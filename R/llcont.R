#' Individual Log-Likelihoods
#'
#' Obtain log-likelihood values associated with individual observations, evaluated at the ML estimates.
#'
#' This is a S3 generic function.
#' Currently, the method is defined for \code{lm}, \code{glm}, \code{glm.nb},
#' \code{clm}, \code{hurdle}, \code{zeroinfl}, \code{mlogit}, \code{nls},
#' \code{polr}, \code{rlm}, \code{lavaan}, \code{vglm}, and \code{mirt} objects.
#'
#' @param x a model object
#' @param \dots arguments passed to specific methods
#'
#' @author Ed Merkle, Dongjun You, and Lennart Schneider
#' 
#' @return An object of class \code{numeric} containing individuals' contributions to the log-likelihood.  The sum of these contributions equals the model log-likelihood.
#'
#' @examples
#' ## Fit gamma glm, check that sum of llcont() equals
#' ## the model loglikelihood:
#' clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
#'                        lot1 = c(118,58,42,35,27,25,21,19,18),
#'                        lot2 = c(69,35,26,21,18,16,13,12,12))
#' gam1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
#' sum(llcont(gam1))
#' logLik(gam1)
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dbinom dgamma dnbinom dnorm dpois
#' @importFrom stats model.frame model.matrix model.response model.weights
#' @importFrom stats weights deviance logLik
#' @importFrom lavaan lavInspect
#'
#' @export
llcont <- function(x, ...) UseMethod("llcont")

#' @export
################################################################
## Getting log-likelihood of glm objects for individual cases
################################################################
llcont.glm <- function(x, ...){
    fam <- x$family$family
    y <- x$y
    wt <- weights(x)
    dev <- deviance(x)
    disp <- dev/sum(wt)
    ## Model-predicted values of y
    #mpreds <- predict(x, type="response")
    mpreds <- fitted(x)

    ## Calculate Log Likelihood associated with this family
    switch(fam,
           binomial = {
             if(is.matrix(y)) {
               n <- apply(y, 1, sum)
               y <- ifelse(n == 0, 0, y[, 1]/n)
             } else {
               n <- rep.int(1, length(y))
             }
             m <- if (any(n > 1)) n else wt
             wt <- ifelse(m > 0, (wt/m), 0)
             dbinom(round(m * y), round(m), mpreds, log = TRUE) * wt
           },
           quasibinomial = {
             NA
           },
           poisson = {
             dpois(y, mpreds, log=TRUE) * wt
           },
           quasipoisson = {
             NA
           },
           gaussian = {
             nobs <- length(y)
             -((log(dev/nobs * 2 * pi) + 1) - log(wt)) / 2
           },
           inverse.gaussian = {
             -((log(disp * 2 * pi) + 1) + 3 * log(y)) / 2
           },
           Gamma = {
             dgamma(y, shape=1/disp, scale=mpreds*disp, log=TRUE) * wt
           })
}

################################################################
## Getting log-likelihood of glm.nb objects for individual cases
################################################################
#' @export
llcont.negbin <- function(x, ...){
  ## Casewise log-likelihood function, adapted from glm.nb()
  loglik <- function(n, th, mu, y, w) w * (lgamma(th +
        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
        log(mu + (y == 0)) - (th + y) * log(th + mu))

  loglik(length(x$y), x$theta, x$fitted.values,
         x$y, x$prior.weights)
}

################################################################
## Getting log-likelihood of clm objects for individual cases
################################################################
#' @export
llcont.clm <- function(x, ...) {
  Class <- class(x)
  Call <- x$call
  mpreds <- fitted(x)
  mf <- model.frame(x)

  if ("(weights)" %in% names(mf)) {
    wts <- mf[["(weights)"]]
  } else {
    wts <- rep(1, nrow(mf))
  }

  if (is.finite(logLik(x))) wts * log(mpreds)
  else Inf
}

################################################################
## Getting log-likelihood of hurdle objects for individual cases
################################################################
#' @export
llcont.hurdle <- function(x, ...) {
  ## If 'model' argument is missing, model.frame needs to be created.
  ## requesting a re-run with 'model=TRUE' would be much easier for us
  if (is.null(x$model)) {
    stop("Please run the model again with the argument 'model=TRUE'")
  } else {
    X <- model.matrix(x, model="count")
    Z <- model.matrix(x, model="zero")
  }

  Y <- x$y
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- Y <= 0
  Y1 <- Y > 0
  weights <- if (is.null(w <- weights(x))) rep.int(1L, n) else w
  offsetx <- x$offset$count
  offsetx <- if (is.null(offsetx)) rep.int(0, n) else offsetx
  offsetz <- x$offset$zero
  offsetz <- if (is.null(offsetz)) rep.int(0, n) else offsetz

  dist <- x$dist$count
  zero.dist <- x$dist$zero
  linkinv <- x$linkinv

  ## individual LL functions
  zeroPoisson <- function(parms) {
    mu <- as.vector(exp(Z %*% parms + offsetz))
    loglik0 <- -mu
    Y0 * weights * loglik0 + ifelse(Y1, weights * log(1 - exp(loglik0)), 0)
  }

  countPoisson <- function(parms) {
    mu <- Y1 * as.vector(exp(X %*% parms + offsetx))
    loglik0 <- -mu
    loglik1 <- Y1 * dpois(Y, lambda = mu, log = TRUE)
    Y1 * weights * loglik1 - ifelse(Y1, weights * log(1 - exp(loglik0)), 0)
  }

  zeroNegBin <- function(parms) {
    mu <- as.vector(exp(Z %*% parms[1:kz] + offsetz))
    theta <- exp(parms[kz + 1])
    loglik0 <- suppressWarnings(dnbinom(0, size = theta,
                                        mu = mu, log = TRUE))
    Y0 * weights * loglik0 +
        ifelse(Y1, weights * log(1 - exp(loglik0)), 0)
  }

  countNegBin <- function(parms) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    theta <- exp(parms[kx + 1])
    loglik0 <- suppressWarnings(dnbinom(0, size = theta,
                                        mu = mu, log = TRUE))
    loglik1 <- suppressWarnings(dnbinom(Y, size = theta,
                                        mu = mu, log = TRUE))
    ifelse(Y1, weights * loglik1 - weights * log(1 - exp(loglik0)), 0)
    }

  zeroGeom <- function(parms) zeroNegBin(c(parms, 0))

  countGeom <- function(parms) countNegBin(c(parms, 0))

  zeroBinom <- function(parms) {
    mu <- as.vector(linkinv(Z %*% parms + offsetz))
    Y0 * weights * log(1 - mu) + Y1 * weights * log(mu)
  }

  ## calculating individual LL
  countDist <- switch(dist,
                      poisson = countPoisson,
                      geometric = countGeom,
                      negbin = countNegBin)
  zeroDist <- switch(zero.dist,
                     poisson = zeroPoisson,
                     geometric = zeroGeom,
                     negbin = zeroNegBin,
                     binomial = zeroBinom)

  loglikfun <- function(parms) {
    countDist(parms[1:(kx + (dist == "negbin"))]) +
        zeroDist(parms[(kx + (dist == "negbin") + 1):
                       (kx + kz + (dist == "negbin") +
                        (zero.dist == "negbin"))])
  }

  if (x$separate) {
    countDist(c(x$coefficients$count,
                if (dist == "negbin") log(x$theta["count"]) else NULL)) +
    zeroDist(c(x$coefficients$zero,
               if (zero.dist == "negbin") log(x$theta["zero"]) else NULL))
  } else {
    loglikfun(c(x$coefficients$count,
                if (dist == "negbin") log(x$theta["count"]) else NULL,
                x$coefficients$zero,
                if (zero.dist == "negbin") log(x$theta["zero"]) else NULL))
  }
}

################################################################
## Getting log-likelihood of zeroinfl objects for individual cases
################################################################
#' @export
llcont.zeroinfl <- function(x, ...) {
  if (is.null(x$model)) {
    stop("Please run the model again with the argument 'model=TRUE'")
  } else {
    X <- model.matrix(x, model="count")
    Z <- model.matrix(x, model="zero")
  }

  Y <- x$y
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- (Y <= 0) * 1
  Y1 <- (Y > 0) * 1
  weights <- if (is.null(w <- weights(x))) rep.int(1L, n) else w
  offsetx <- x$offset$count
  offsetx <- if (is.null(offsetx)) rep.int(0, n) else offsetx
  offsetz <- x$offset$zero
  offsetz <- if (is.null(offsetz)) rep.int(0, n) else offsetz

  dist <- x$dist
  linkinv <- x$linkinv

  ## individual LL functions
  ziPoisson <- function(parms) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    loglik0 <- log(phi + exp(log(1 - phi) - mu))
    loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
    Y0 * weights * loglik0 + Y1 * weights * loglik1
  }

  ziNegBin <- function(parms) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) +
                   suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) +
               suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    Y0 * weights * loglik0 + Y1 * weights * loglik1
  }

  ziGeom <- function(parms) ziNegBin(c(parms, 0))

  ## calculating individual LL
  loglikfun <- switch(dist,
                      poisson = ziPoisson,
                      geometric = ziGeom,
                      negbin = ziNegBin)

  loglikfun(c(x$coefficients$count,
              x$coefficients$zero,
              if (dist == "negbin") log(x$theta) else NULL))
} ## end of llcont.zeroinfl

################################################################
## Getting log-likelihood of lm objects for individual cases
################################################################
#' @export
llcont.lm <- function (x, ...) {
  if (inherits(x, "mlm"))
    stop("'logLik' does not support multiple responses")
  res <- x$residuals
  p <- x$rank
  N <- length(res)
  s <- summary(x)$sigma
  ## calculate s2 for ML method
  sml2 <- (s * sqrt((N - p) / N))^2
  if (is.null(w <- x$weights)) {
    w <- rep.int(1, N)
  } else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }

  0.5 * (log(w) - (log(2 * pi) + log(sml2) + (w * res^2)/sml2))
}

################################################################
## Getting log-likelihood of mlogit objects for individual cases
################################################################
#' @export
llcont.mlogit <- function(x, ...) log(fitted(x))

################################################################
## Getting log-likelihood of nls objects for individual cases
################################################################
#' @export
llcont.nls <- function (x, ...) {
  res <- x$m$resid()
  N <- length(res)
  s <- summary(x)$sigma
  N_p <- summary(x)$df[2]
  sml2 <- (s * sqrt((N_p) / N))^2
  if (is.null(w <- x$weights))
    w <- rep_len(1, N)
  zw <- w == 0

  -0.5 * (log(2 * pi) + log(sml2) - log(w + zw) + (w * res^2)/sml2)
}

################################################################
## Getting log-likelihood of polr  objects for individual cases
################################################################
#' @export
llcont.polr <- function(x, ...) {
  m <- x$model
  y <- unclass(model.response(m))
  wherey <- matrix(c(as.numeric(names(y)), y), ncol=2)
  idx <- matrix(0, nrow=length(y), ncol=length(x$lev))
  idx[wherey] <- 1

  model.weights(m) * log(apply(idx * x$fitted.values, 1, sum))
}

################################################################
## Getting log-likelihood of rlm objects for individual cases
################################################################
#' @export
llcont.rlm <- function (x, ...) {
  res <- x$residuals
  p <- x$rank
  N <- length(res)
  if (is.null(w <- x$weights)) {
    w <- rep.int(1, N)
  } else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }

  s2 <- sum(w*res^2)/N
  0.5 * (log(w) - (log(2 * pi) + log(s2) + (w * res^2)/s2))
}

################################################################
## Getting log-likelihood of lavaan objects for individual cases
################################################################
#' @export
llcont.lavaan <- function(x, ...){
  ## make sure this is multivariate normal likelihood
  if (lavInspect(x, "options")$estimator != "ML" |
      lavInspect(x, "options")$se != "standard"){
      stop("nonnest2 only works for lavaan models fit via ML\n  (assuming multivariate normality, with no robust SEs).")
  }
  mispatts <- lavInspect(x, "patterns")
  if(any(class(mispatts) == "list")){
    npatts <- max(sapply(mispatts, nrow))
  } else {
    npatts <- nrow(mispatts)
  }
  if (lavInspect(x, "options")$missing != "ml" & npatts > 1){
      stop("nonnest2 does not work with pairwise/listwise deletion.\n  refit the model with missing='ml'.")
  }
  samplestats <- x@SampleStats
  ntab <- lavInspect(x, "norig")
  llvec <- rep(NA, sum(lavInspect(x, "norig")))
  ngroups <- lavInspect(x, "ngroups")

  for(g in 1:ngroups) {
    if (ngroups > 1){
      moments <- fitted(x)[[g]]
      grpind <- lavInspect(x, "case.idx")[[g]]
    } else {
      moments <- fitted(x)
      grpind <- lavInspect(x, "case.idx")
    }
    Sigma.hat <- unclass(moments$cov)
    ## To ensure it is symmetric; needed?
    Sigma.hat <- (Sigma.hat + t(Sigma.hat))/2

    if(!samplestats@missing.flag) { # complete data
      if(x@Model@meanstructure) { # mean structure
        Mu.hat <- unclass(moments$mean)
      } else {
        ## set mean structure to sample estimates
        Mu.hat <- apply(x@Data@X[[g]], 2, mean)
      }
      llvec[grpind] <- dmvnorm(x@Data@X[[g]], Mu.hat, Sigma.hat, log=TRUE)

      ## subtract logl.X
      x.idx <- samplestats@x.idx[[g]]
      if(!is.null(x.idx) && length(x.idx) > 0L){
        Mu.X <- samplestats@mean.x[[g]]
        Sigma.X <- samplestats@cov.x[[g]]
        if(is.null(Mu.X)){
          Mu.X <- Mu.hat[x.idx]
        }
        if(is.null(Sigma.X)){
          Sigma.X <- Sigma.hat[x.idx, x.idx, drop=FALSE]
        }

        if(length(x.idx) == 1){
          tmpll.x <- dnorm(x@Data@X[[g]][,x.idx], Mu.X, sqrt(Sigma.X), log=TRUE)
        } else {
          tmpll.x <- dmvnorm(x@Data@X[[g]][,x.idx], Mu.X, Sigma.X, log=TRUE)
        }
        if(inherits(tmpll.x, "try-error")) tmpll.x <- NA
        llvec[grpind] <- llvec[grpind] - tmpll.x
      }
      
    } else { # incomplete data
      nsub <- ntab[g]
      M <- samplestats@missing[[g]]
      Mp <- x@Data@Mp[[g]]
      tmpll <- rep(NA, nsub)

      Mu.hat <- unclass(moments$mean)
      nvar <- ncol(moments$cov)

      for(p in 1:length(M)) {
        ## Data
        case.idx <- Mp$case.idx[[p]]
        var.idx <- M[[p]][["var.idx"]]
        X <- x@Data@X[[g]][case.idx, var.idx, drop = FALSE]

        ## avoid fail for one observed variable
        if(sum(var.idx) == 1){
          tmpll[case.idx] <- dnorm(X, Mu.hat[var.idx], sqrt(Sigma.hat[var.idx,var.idx]), log=TRUE)
        } else {
          tmpll[case.idx] <- dmvnorm(X, Mu.hat[var.idx], Sigma.hat[var.idx, var.idx], log=TRUE)
        }

        ## subtract logl.X
        x.idx <- samplestats@x.idx[[g]]
        obsx <- which(x.idx %in% var.idx)
        x.idx <- x.idx[obsx]
        if(!is.null(x.idx) && length(x.idx) > 0L){
          Mu.X <- Mu.hat[x.idx]
          Sigma.X <- Sigma.hat[x.idx, x.idx, drop=FALSE]

          if(length(x.idx) == 1){
            tmpll.x <- dnorm(x@Data@X[[g]][,x.idx], Mu.X, sqrt(Sigma.X), log=TRUE)
          } else {
            tmpll.x <- dmvnorm(X[,x.idx], Mu.X, Sigma.X, log=TRUE)
          }
          if(inherits(tmpll.x, "try-error")) tmpll.x <- NA
          tmpll[case.idx] <- tmpll[case.idx] - tmpll.x
        }
      }

      llvec[grpind] <- tmpll
    } # incomplete
  } # group
  ## remove obs with no variables
  empidx <- unlist(lavInspect(x, 'empty.idx'))
  if(length(empidx) > 0){
    llvec <- llvec[-empidx]
  }
  
  llvec
}

################################################################
## Getting log-likelihood of vglm objects for individual cases
################################################################
#' @export
llcont.vglm <- function(x, ...){
  logLik(x, summation = FALSE)
}

########################################################################
## individual log-likelihood of MultipleGroupClass objects (mirt function)
########################################################################
#' @export
llcont.MultipleGroupClass <- function(x, ...) {
  dat <- apply(x@Data$data, 1, paste, collapse = "")
  tab <- apply(x@Data$tabdata, 1, paste, collapse = "")
  ind <- match(dat, tab)
  g  <- x@Data$group
  gN <- x@Data$groupNames
  llcont <- numeric(length(g))
  for(a in seq_along(gN)) {
    llcont[g == gN[a]] <- log(x@Internals$Pl[[a]][ind[g == gN[[a]]]])
  }
  return(llcont)
}

#' @export
llcont.SingleGroupClass <- function(x, ...) {
  dat <- apply(x@Data$data, 1, paste, collapse = "")
  tab <- apply(x@Data$tabdata, 1, paste, collapse = "")
  ind <- match(dat, tab)
  llcont <- log(x@Internals$Pl[ind])
  return(llcont)
}
