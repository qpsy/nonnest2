#' Information Criteria Confidence Intervals
#'
#' \code{icci} calculates confidence intervals of BIC.
#' Functionality is available for models of classes
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
#' @param conf.level confidence level of the interval
#'
#' @author Ed Merkle and Dongjun You
#'
#' @return an object of class \code{icci} containing test results.
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
#'
#' ## CI for BIC
#' icci(house2, house1)
#'
#' require(pscl)
#' bio1 <- glm(art ~ fem + mar + phd + ment, family=poisson, data=bioChemists)
#' bio2 <- hurdle(art ~ fem + mar + phd + ment, data=bioChemists)
#' bio3 <- zeroinfl(art ~ fem + mar + phd + ment, data=bioChemists)
#' icci(bio2, bio1)
#' icci(bio3, bio1)
#' icci(bio3, bio2)

#' require(lavaan)
#' HS.model <- 'visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit1 <- cfa(HS.model, data=HolzingerSwineford1939)
#' fit2 <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
#' icci(fit1, fit2, 0.05)
#' }
#'
#' @export
icci <- function(object1, object2, conf.level=.95) {
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

  bicA <- BIC(object1)
  bicB <- BIC(object2)

  aicA <- AIC(object1)
  aicB <- AIC(object2)

  ## Interval for BIC differences
  ## if (classA %in% c("hurdle", "zeroinfl", "mlogit")) {
  ##   lls <- logLik(object1)
  ##   nos <- length(fitted(object1))
  ##   bicA <- -2 * as.numeric(lls) + log(nos) * attr(lls, "df")
  ## } else {
  ##   bicA <- BIC(object1)
  ## }

  ## if (classB %in% c("hurdle", "zeroinfl", "mlogit")) {
  ##   lls <- logLik(object2)
  ##   nos <- length(fitted(object2))
  ##   bicB <- -2 * as.numeric(lls) + log(nos) * attr(lls, "df")
  ## } else {
  ##   bicB <- BIC(object2)
  ## }

  bicdiff <- bicA - bicB
  aicdiff <- aicA - aicB
  alpha <- 1 - conf.level

  ## BIC CI
  BICci <- bicdiff + qnorm(c(alpha/2,(1-alpha/2)))*sqrt(n * 4 * omega.hat.2)
  ## AIC CI
  AICci <- aicdiff + qnorm(c(alpha/2,(1-alpha/2)))*sqrt(n * 4 * omega.hat.2)

  rval <- list(class = list(class1=classA, class2=classB),
               call = list(call1=callA, call2=callB),
               BIC = list(BIC1=bicA, BIC2=bicB),
               BICci = BICci,
               AIC = list(AIC1=aicA, AIC2=aicB),
               AICci = AICci,
               confLevel = conf.level)
  class(rval) <- "icci"
  return(rval)
}

################################################################
## print method for icci
################################################################
#' @method print icci
#' @export
print.icci <- function(x, ...) {
  cat("\nModel 1 \n")
  cat(" Class:", x$class$class1, "\n")
  cat(" Call:", deparse(x$call$call1), fill=TRUE)
  cat(" AIC:", formatC(x$AIC$AIC1, digits=3L, format="f"), "\n")
  cat(" BIC:", formatC(x$BIC$BIC1, digits=3L, format="f"), "\n")
  cat("Model 2 \n")
  cat(" Class:", x$class$class2, "\n")
  cat(" Call:", deparse(x$call$call2), "\n", fill=TRUE)
  cat(" AIC:", formatC(x$AIC$AIC2, digits=3L, format="f"), "\n")
  cat(" BIC:", formatC(x$BIC$BIC2, digits=3L, format="f"), "\n\n")

  if (any(c(x$class$class1, x$class$class2) %in% c("hurdle", "zeroinfl", "mlogit"))) {
    warning("\n Currently, BIC cannot be calculated for the objects of hurdle, zeroinfl and mlogit.", call.=FALSE)

    cat(x$confLevel * 100,
        "% Confidence Interval of AIC difference (AICdiff = AIC1 - AIC2) \n", sep="")
    cat("  ", formatC(x$AICci[1], digits=3L, format="f"), " < ", "AICdiff",
        " < ", formatC(x$AICci[2], digits=3L, format="f"), "\n", sep="")
  } else {
    cat(x$confLevel * 100,
        "% Confidence Interval of AIC difference (AICdiff = AIC1 - AIC2) \n", sep="")
    cat("  ", formatC(x$AICci[1], digits=3L, format="f"), " < ", "AICdiff",
        " < ", formatC(x$AICci[2], digits=3L, format="f"), "\n\n", sep="")

    cat(x$confLevel * 100,
        "% Confidence Interval of BIC difference (BICdiff = BIC1 - BIC2) \n", sep="")
    cat("  ", formatC(x$BICci[1], digits=3L, format="f"), " < ", "BICdiff",
        " < ", formatC(x$BICci[2], digits=3L, format="f"), "\n", sep="")
  }
}
