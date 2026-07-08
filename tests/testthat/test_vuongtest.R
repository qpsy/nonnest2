context("Vuong (1989) statistic integrity")

## These tests pin the numeric core of vuongtest() and icci() to the equations
## in Vuong (1989, Econometrica 57, 307-333) so that future refactors cannot
## silently alter a formula.  Each package value is checked two ways:
##   (1) against an INDEPENDENT recomputation from the published equations
##       (guards the vuongtest()/icci() arithmetic), and
##   (2) against hard-coded reference anchors computed once from this fixture
##       (guards llcont() and the eigenvalue/imhof path).
## The fixture uses only base-R glm() so no heavy model-fitting packages are
## required.  See docs/papers/SOURCES.md for the equation-by-equation mapping.

make_fixture <- function() {
  set.seed(42)
  n <- 200
  x <- rnorm(n)
  z <- rnorm(n)
  y <- rpois(n, exp(0.3 + 0.5 * x))
  d <- data.frame(y = y, x = x, z = z)
  list(
    d  = d,
    m1 = glm(y ~ x,     family = poisson, data = d),  # 2 parameters
    m2 = glm(y ~ z,     family = poisson, data = d),  # 2 parameters
    m3 = glm(y ~ x + z, family = poisson, data = d)   # 3 parameters
  )
}

test_that("non-nested Vuong statistic matches Vuong (1989) Eq (4.2) and the z-form", {
  fx <- make_fixture()
  res <- vuongtest(fx$m1, fx$m2)

  ## Independent recomputation straight from the paper -----------------------
  m  <- llcont(fx$m1) - llcont(fx$m2)          # individual LR contributions
  n  <- length(m)
  omega2 <- (n - 1) / n * var(m)               # Eq (4.2): ML variance of m_i
  lr     <- sum(m)                             # Eq (6.4): sum of LR contributions
  z      <- (1 / sqrt(n)) * lr / sqrt(omega2)  # non-nested statistic ~ N(0,1)

  expect_equal(res$omega,      omega2)
  expect_equal(res$LRTstat,    z)
  expect_equal(res$p_LRT$A,    pnorm(z, lower.tail = FALSE))
  expect_equal(res$p_LRT$B,    pnorm(z))

  ## Hard-coded anchors (catch damage upstream of the arithmetic, e.g. llcont)
  expect_equal(res$omega,   0.4077911447, tolerance = 1e-6)
  expect_equal(res$LRTstat, 2.4396085692, tolerance = 1e-6)
  expect_equal(res$p_LRT$A, 0.0073515919, tolerance = 1e-5)
  ## variance / distinguishability test (weighted chi-square, Vuong Thm 4.1)
  expect_true(res$p_omega < 1e-5)
})

test_that("AIC/BIC adjustments follow Vuong (1989) Eq (5.7) corrections", {
  fx <- make_fixture()
  ## m3 (3 params) vs m2 (2 params): parameter-count difference is non-zero,
  ## so the correction term actually bites.
  ## suppressWarnings() only silences benign CompQuadForm::imhof numerical
  ## notes ("Qq + abserr is positive"); it does not affect the statistics.
  base <- suppressWarnings(vuongtest(fx$m3, fx$m2))
  aic  <- suppressWarnings(vuongtest(fx$m3, fx$m2, adj = "aic"))
  bic  <- suppressWarnings(vuongtest(fx$m3, fx$m2, adj = "bic"))

  m  <- llcont(fx$m3) - llcont(fx$m2)
  n  <- length(m)
  omega2 <- (n - 1) / n * var(m)
  lr     <- sum(m)
  dpar   <- length(coef(fx$m3)) - length(coef(fx$m2))   # = 1
  denom  <- sqrt(n) * sqrt(omega2)

  ## AIC correction subtracts (p - q); BIC subtracts (p - q) * log(n) / 2
  expect_equal(base$LRTstat, lr / denom)
  expect_equal(aic$LRTstat,  (lr - dpar) / denom)
  expect_equal(bic$LRTstat,  (lr - dpar * log(n) / 2) / denom)

  ## adjustments must strictly shrink the statistic here (dpar > 0, lr > 0)
  expect_true(bic$LRTstat < aic$LRTstat)
  expect_true(aic$LRTstat < base$LRTstat)
})

test_that("nested Vuong statistic equals 2 * LR (Vuong 1989 Eq 3.6)", {
  fx <- make_fixture()
  res <- vuongtest(fx$m3, fx$m1, nested = TRUE)   # m1 (y~x) nested in m3 (y~x+z)

  lr <- sum(llcont(fx$m3) - llcont(fx$m1))
  expect_equal(res$LRTstat, 2 * lr)               # robust LR statistic
  expect_equal(res$LRTstat, 0.5772371622, tolerance = 1e-6)
  expect_true(res$p_LRT$A >= 0 && res$p_LRT$A <= 1)
})

test_that("icci AIC/BIC interval half-width equals 2*sqrt(n)*omega (Merkle, You & Preacher 2016)", {
  fx <- make_fixture()
  ic <- icci(fx$m1, fx$m2)

  m  <- llcont(fx$m1) - llcont(fx$m2)
  n  <- length(m)
  omega2  <- (n - 1) / n * var(m)
  se      <- sqrt(n * 4 * omega2)                 # SD of the IC difference
  aicdiff <- AIC(fx$m1) - AIC(fx$m2)
  bicdiff <- AIC(fx$m1, k = log(n)) - AIC(fx$m2, k = log(n))

  expect_equal(as.numeric(ic$AICci), aicdiff + qnorm(c(.025, .975)) * se)
  expect_equal(as.numeric(ic$BICci), bicdiff + qnorm(c(.025, .975)) * se)
  ## AIC and BIC CIs share the same half-width (only the centre differs)
  expect_equal(diff(ic$AICci), diff(ic$BICci))
  expect_equal(as.numeric(ic$AICci), c(-79.46472283, -8.663301608), tolerance = 1e-5)
})
