#' Vuong Tests for Model Comparison
#'
#' Test pairs of models using Vuong's (1989) <DOI:10.2307/1912557> theory.  This includes
#' a test of model distinguishability and a test of model fit.
#'
#' For non-nested models, the test of distinguishability indicates whether or
#' not the models can possibly be distinguished on the basis of the observed
#' data.  The LRT then indicates whether or not one model fits better
#' than another.
#'
#' For nested models (\code{nested=TRUE}), both tests serve as robust
#' alternatives to the classical likelihood ratio tests.  In this case,
#' the \code{adj} argument is ignored.
#'
#' Users should take care to ensure that the two models have
#' the same dependent variable (or, for lavaan objects, identical
#' modeled variables), with observations ordered identically within
#' each model object.  Assuming the same data matrix is used to fit each model,
#' observation ordering should generally be identical.  There are currently
#' no checks for this, however.
#'
#'
#' @param object1 a model object
#' @param object2 a model object
#' @param nested if \code{TRUE}, models are assumed to be nested
#' @param adj Should an adjusted test statistic be calculated?  Defaults to \dQuote{none}, with possible adjustments being \dQuote{aic} and \dQuote{bic}
#'
#' @author Ed Merkle and Dongjun You
#'
#' @return an object of class \code{vuongtest} containing test results.
#'
#' @references
#'
#' Vuong, Q. H. (1989).  Likelihood ratio tests for model selection and non-nested hypotheses.  \emph{Econometrica, 57}, 307-333. <DOI:10.2307/1912557>
#'
#' Merkle, E. C., You, D., & Preacher, K. (2016). Testing non-nested structural equation models.  \emph{Psychological Methods, 21}, 151-163. <DOI:10.1037/met0000038>
#'
#' @examples
#' \dontrun{
#' ## Count regression comparisons
#' require(MASS)
#' house1 <- glm(Freq ~ Infl + Type + Cont, family=poisson, data=housing)
#' house2 <- glm(Freq ~ Infl + Sat, family=poisson, data=housing)
#' house3 <- glm(Freq ~ Infl, family=poisson, data=housing)
#' ## house3 is nested within house1 and house2
#' anova(house3, house1, test="Chisq")
#' anova(house3, house2, test="Chisq")
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
#'
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
#' @importFrom stats coef pnorm var vcov
#' @importMethodsFrom lavaan coef fitted logLik vcov
#' @export
vuongtest <- function(object1, object2, nested=FALSE, adj="none") {

  ## check objects, issue warnings/errors, get classes/calls
  obinfo <- check.obj(object1, object2)
  callA <- obinfo$callA; classA <- obinfo$classA
  callB <- obinfo$callB; classB <- obinfo$classB
    
  llA <- llcont(object1)
  llB <- llcont(object2)

  ## If nested==TRUE, decide which is the full model by seeing
  ## which has the larger log-likelihood.  From here on,
  ## object1 is the full model
  if(nested){
      if(sum(llB, na.rm = TRUE) > sum(llA, na.rm = TRUE)){
          tmp <- object1
          object1 <- object2
          object2 <- tmp
          tmp <- llA
          llA <- llB
          llB <- tmp
          tmp <- callA
          callA <- callB
          callB <- tmp
      }
  }

  ## Eq (4.2)
  nmis <- sum(is.na(llA)) # (missing all data)
  n <- NROW(llA) - nmis
  omega.hat.2 <- (n-1)/n * var(llA - llB, na.rm = TRUE)

  ## Get p-value of weighted chi-square dist
  lamstar <- calcLambda(object1, object2, n)

  ## Note: dr package requires non-negative weights, which
  ##       does not help when nested==TRUE
  ## tmp <- dr.pvalue(lamstar2, n * omega.hat.2)
  ## pOmega <- tmp[[4]]
  pOmega <- imhof(n * omega.hat.2, lamstar^2)$Qq

  ## Calculate and test LRT; Eq (6.4)
  lr <- sum(llA - llB, na.rm = TRUE)
  teststat <- (1/sqrt(n)) * lr/sqrt(omega.hat.2)

  ## Adjustments to test statistics
  ## FIXME lavaan equality constraints; use df instead?
  if(classA %in% c("SingleGroupClass", "MultipleGroupClass")){
    nparA <- mirt::extract.mirt(object1, "nest")
  } else {
    nparA <- length(coef(object1))
  }
  if(classB %in% c("SingleGroupClass", "MultipleGroupClass")){
    nparB <- mirt::extract.mirt(object2, "nest")
  } else {
    nparB <- length(coef(object2))
  }
  
  if(adj=="aic"){
    teststat <- teststat - (nparA - nparB)
  }
  if(adj=="bic"){
    teststat <- teststat - (nparA - nparB) * log(n)/2
  }

  ## Null distribution and test stat depend on nested
  if(nested){
      teststat <- 2 * lr
      ## lamstar is negative due to the ordering of objects in calcLambda();
      ## we could get the same thing via calcLambda(object2, object1, n)
      pLRTA <- imhof(teststat, -lamstar)[[1]]

      ## Compare different approximations
      ## Empirical
      ## tmp <- matrix(rchisq(1000 * length(lamstar), 1), 1000, length(lamstar))
      ## edist <- apply(tmp, 1, function(x) sum(lamstar*x))
      ## print(summary(edist))
      ## imhof without negative weights
      ## print(imhof(teststat, lamstar)[[1]])
      ## davies without negative weights
      ## print(davies(teststat, lamstar)$Qq)

      pLRTB <- NA
  } else {
      ## Two 1-tailed p-values from a normal:
      pLRTA <- pnorm(teststat, lower.tail=FALSE)
      pLRTB <- pnorm(teststat)
  }

  rval <- list(omega = omega.hat.2, p_omega = pOmega,
               LRTstat = teststat,
               p_LRT = list(A=pLRTA, B=pLRTB),
               class = list(class1=classA, class2=classB),
               call = list(call1=callA, call2=callB),
               nested = nested)
  class(rval) <- "vuongtest"
  return(rval)
}

################################################################
## A, B as defined in Vuong Eq (2.1) and (2.2)
################################################################
calcAB <- function(object, n){
  ## Eq (2.1)
  if(class(object)[1] == "lavaan"){
    tmpvc <- vcov(object)
    dups <- duplicated(colnames(tmpvc))
    tmpvc <- n * tmpvc[!dups,!dups]
    ## to throw error if complex constraints
    ## (NB we should eventually just use this instead of dups)
    if(nrow(object@Model@ceq.JAC) > 0){
      vcerr <- vcov(object, remove.duplicated=TRUE)
    }
    #if(nrow(object@Model@ceq.JAC) > 0){
    #  A <- vcov(object, remove.duplicated=TRUE)
    #} else {
    #  A <- vcov(object)
    #}
  } else if(class(object) %in% c("lm", "glm", "nls")){
    scaling <- summary(object)$sigma
    if(is.null(scaling)){
      scaling <- 1
    } else {
      scaling <- scaling^2
    }
    tmpvc <- n * vcov(object)
  } else {
    tmpvc <- n * vcov(object)
    ## in case mirt vcov was not estimated
    if(nrow(tmpvc) == 1 & is.na(tmpvc[1,1])) stop("Please re-estimate the mirt model with SE=TRUE")
  }
  A <- chol2inv(chol(tmpvc))

  ## Eq (2.2)
  if(class(object)[1] == "lavaan"){
    sc <- estfun(object, remove.duplicated=TRUE)
  } else if (class(object)[1] %in% c("SingleGroupClass", "MultipleGroupClass")){
    sc <- mirt::estfun.AllModelClass(object)
  } else if(class(object) %in% c("lm", "glm", "nls")){
    sc <- (1/scaling) * estfun(object)
  } else {
    sc <- estfun(object)
  }
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
  ## using Re
  Re(lamstar)
}

################################################################
## print method for vuongtest
################################################################
#' @method print vuongtest
#' @export
print.vuongtest <- function(x, ...) {
  cat("\nModel 1 \n")
  cat(" Class:", x$class$class1, "\n")
  model1call <- deparse(x$call$call1)
  cat(" Call: ", model1call[1], if (length(model1call) > 1) "...\n" else "\n", sep="")
  cat("\nModel 2 \n")
  cat(" Class:", x$class$class2, "\n")
  model2call <- deparse(x$call$call2)
  cat(" Call: ", model2call[1], if (length(model2call) > 1) "...\n" else "\n", sep="")

  cat("\nVariance test \n")
  cat("  H0: Model 1 and Model 2 are indistinguishable", "\n")
  cat("  H1: Model 1 and Model 2 are distinguishable", "\n")
  cat("    w2 = ", formatC(x$omega, digits=3L, format="f"), ",   ",
      "p = ", format.pval(x$p_omega, digits=3L), "\n\n", sep="")

  if(x$nested){
      cat("Robust likelihood ratio test of distinguishable models \n")
      cat("  H0: Model 2 fits as well as Model 1 \n")
      cat("  H1: Model 1 fits better than Model 2 \n")
      cat("    LR = ", formatC(x$LRTstat, digits=3L, format="f"), ",   ",
          "p = ", format.pval(x$p_LRT[[1]], digits=3L), "\n", sep="")
  } else {
      cat("Non-nested likelihood ratio test \n")
      cat("  H0: Model fits are equal for the focal population \n")
      cat("  H1A: Model 1 fits better than Model 2 \n")
      cat("    z = ", formatC(x$LRTstat, digits=3L, format="f"), ",   ",
          "p = ", format.pval(x$p_LRT[[1]], digits=3L), "\n", sep="")
      cat("  H1B: Model 2 fits better than Model 1 \n")
      cat("    z = ", formatC(x$LRTstat, digits=3L, format="f"), ",   ",
          "p = ", format.pval(x$p_LRT[[2]], digits=4L), "\n", sep="")
  }
}

###################################################################
## check objects prior to doing computations and possibly warn/stop
## also return calls/classes
###################################################################
check.obj <- function(object1, object2) {
  classA <- class(object1)[1L]
  classB <- class(object2)[1L]
  
  if(isS4(object1)){
    if(classA %in% c("SingleGroupClass", "MultipleGroupClass")){
      callA <- object1@Call
      ## recommended vcov type for mirt models:
      if(object1@Options$SE.type != "Oakes") warning("SE.type='Oakes' is recommended for mirt models")
    } else {
      callA <- object1@call
    }
  } else {
    callA <- object1$call
  }
  if(isS4(object2)){
    if(classB %in% c("SingleGroupClass", "MultipleGroupClass")){
      callB <- object2@Call
      if(object2@Options$SE.type != "Oakes") warning("SE.type='Oakes' is recommended for mirt models")
    } else {
      callB <- object2@call
    }
  } else {
    callB <- object2$call
  }

  list(classA = classA, classB = classB, callA = callA, callB = callB)
}  
  

.onAttach <- function(...) {
  version <- read.dcf(file=system.file("DESCRIPTION", package="nonnest2"), fields="Version")
  packageStartupMessage(paste0("This is nonnest2 ", version, ".\nnonnest2 has not been tested with all combinations of model classes."))
}
