#' Vuong Tests for Model Comparison
#'
#' \code{vuongtest} formally tests pairs of models using theory provided
#' by Vuong (1989).  Functionality is available for models of classes
#' lm, glm, glm.nb, clm, hurdle, zeroinfl, mlogit, nls, polr, rlm,
#' and lavaan.
#'
#' Users should take care to ensure that the two models have
#' the same dependent variable (or, for lavaan objects, identical
#' modeled variables), with observations ordered identically within
#' each model object.
#'
#' @param object1 a model object
#' @param object2 a model object
#'
#' @author Ed Merkle and Dongjun You
#'
#' @return an object of class \code{vuongtest} containing test results.
#'
#' @references
#'
#' Vuong, Q. H. (1989).  Likelihood ratio tests for model selection and non-nested hypotheses.  \emph{Econometrica, 57}, 307-333.
#'
#' Merkle, E. C., You, D., & Preacher, K. (2014). Testing non-nested structural equation models.  \emph{Manuscript under review}.
#'
#' @examples
#' \dontrun{
#' require(MASS)
#' house1 <- glm(Freq ~ Infl + Type + Cont, family=poisson, data=housing)
#' house2 <- glm(Freq ~ Infl + Sat, family=poisson, data=housing)
#' house3 <- glm(Freq ~ Infl, family=poisson, data=housing)
#' ## house3 is nested within house1 and house2
#' anova(house3, house1, test="Chisq")
#' anova(house2, house1, test="Chisq")
#'
#' ## house 2 is not nested in house1, so this test is invalid
#' anova(house2, house1, test="Chisq")
#'
#' ## Use vuongtest() instead
#' vuongtest(house2, house1)
#'
#' ## Application to models with different distributional assumptions
#' require(pscl)
#' bio1 <- glm(art ~ fem + mar + phd + ment, family=poisson, data=bioChemists)
#' bio2 <- hurdle(art ~ fem + mar + phd + ment, data=bioChemists)
#' bio3 <- zeroinfl(art ~ fem + mar + phd + ment, data=bioChemists)
#' vuongtest(bio2, bio1)
#' vuongtest(bio3, bio1)
#' vuongtest(bio1, bio2)
#' vuongtest(bio1, bio3)
#' vuongtest(bio3, bio2)

#' ## Application to latent variable models
#' require(lavaan)
#' HS.model <- 'visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit1 <- cfa(HS.model, data=HolzingerSwineford1939)
#' fit2 <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
#' vuongtest(fit1, fit2)
#' }
#'
#' @importFrom sandwich estfun
#' @importFrom CompQuadForm imhof
#' @export
vuongtest <- function(object1, object2) {
  classA <- class(object1)[1L]
  classB <- class(object2)[1L]
  callA <- if (isS4(object1)) object1@call else object1$call
  callB <- if (isS4(object2)) object2@call else object2$call

  llA <- llcont(object1)
  llB <- llcont(object2)

  if (!isTRUE(all.equal(sum(llA), as.numeric(logLik(object1)))))
    stop("The individual log-likelihoods do not sum up to the log-likelihood. Please report your model and object to the maintainer.")

  if (!isTRUE(all.equal(sum(llB), as.numeric(logLik(object2)))))
    stop("The individual log-likelihoods do not sum up to the log-likelihood. Please report your model and object to the maintainer.")

  ## Eq (4.2)
  n <- NROW(llA)
  omega.hat.2 <- (n-1)/n * var(llA - llB)

  ## Get p-value of weighted chi-square dist
  ## Need to install the dr package
  lamstar2 <- calcLambda(object1, object2, n)
  ## tmp <- dr.pvalue(lamstar2, n * omega.hat.2)
  ## pOmega <- tmp[[4]]
  ## Alternative:
  ## library(CompQuadForm)
  pOmega <- imhof(n * omega.hat.2, lamstar2)[[1]]

  ## Calculate and test LRT; Eq (6.4)
  lr <- sum(llA - llB)
  teststat <- (1/sqrt(n)) * lr/sqrt(omega.hat.2)

  ## Two 1-tailed p-values:
  pLRTA <- pnorm(teststat, lower.tail=FALSE)
  pLRTB <- pnorm(teststat)

  ## Interval for BIC differences
  #bicA <- BIC(object1)
  #bicB <- BIC(object2)
  #bicdiff <- bicA - bicB

  ## BIC CI
  #ci <- bicdiff + qnorm(c(alpha/2,(1-alpha/2)))*sqrt(n * 4 * omega.hat.2)

  rval <- list(omega = omega.hat.2, p_omega = pOmega,
               z_LRT = teststat,
               p_LRT = list(A=pLRTA, B=pLRTB),
               class = list(class1=classA, class2=classB),
               call = list(call1=callA, call2=callB))#, BICCI=ci)
  class(rval) <- "vuongtest"
  return(rval)
}

################################################################
## A, B as defined in Vuong Eq (2.1) and (2.2)
################################################################
calcAB <- function(object, n){
  ## Eq (2.1)
  A <- chol2inv(chol(n * vcov(object)))

  ## Eq (2.2)
  sc <- estfun(object)
  sc.cp <- crossprod(sc)/n
  B <- matrix(sc.cp, nrow(A), nrow(A))

  list(A=A, B=B, sc=sc)
}

## a function to get the cross-product from Eq (2.7)
calcBcross <- function(sc1, sc2, n){
  ## Get Eq (2.7)
  crossprod(sc1, sc2)/n
}


################################################################
## Calculating W, Vuong Eq (3.6)
################################################################
calcLambda <- function(object1, object2, n) {
  AB1 <- calcAB(object1, n)
  AB2 <- calcAB(object2, n)
  Bc <- calcBcross(AB1$sc, AB2$sc, n)

  W <- cbind(rbind(-AB1$B %*% chol2inv(chol(AB1$A)),
                   t(Bc) %*% chol2inv(chol(AB1$A))),
             rbind(-Bc %*% chol2inv(chol(AB2$A)),
                   AB2$B %*% chol2inv(chol(AB2$A))))

  lamstar <- eigen(W, only.values=TRUE)$values
  ## Discard imaginary part, as it only occurs for tiny eigenvalues?
  ## using Mod
  Mod(lamstar)^2
}

################################################################
## print method for vuongtest (under construction)
################################################################
#' @method print vuongtest
#' @export
print.vuongtest <- function(x, ...) {
  cat("\nModel 1 \n")
  cat(" Class:", x$class$class1, "\n")
  cat(" Call:", deparse(x$call$call1), "\n", fill=TRUE)
  cat("Model 2 \n")
  cat(" Class:", x$class$class2, "\n\n")
  cat(" Call:", deparse(x$call$call2), "\n\n", fill=TRUE)

  cat("Variance test \n")
  cat("  H0: Model 1 and Model 2 are equivalent for the focal population", "\n")
  cat("  w2 = ", formatC(x$omega, digits=3L, format="f"), ",   ",
      "p = ", format.pval(x$p_omega, digits=3L), "\n\n", sep="")

  cat("Non-nested likelihood ratio test \n")
  cat("  H0: Model fits are equal for the focal population \n")
  cat("  H1A: Model 1 fits better than Model 2 \n")
  cat("  z = ", formatC(x$z_LRT, digits=3L, format="f"), ",   ",
      "p = ", format.pval(x$p_LRT[[1]], digits=3L), "\n", sep="")
  cat("  H1B: Model 2 fits better than Model 1 \n")
  cat("  z = ", formatC(x$z_LRT, digits=3L, format="f"), ",   ",
      "p = ", format.pval(x$p_LRT[[2]], digits=4L), "\n", sep="")
}


.onAttach <- function(...) {
  packageStartupMessage("  This is nonnest2 0.1
  nonnest2 is BETA software! Please report any bugs.")
}
