# Canonical source equations for nonnest2

This file documents the published equations that the numeric core of
`nonnest2` implements, so that the code can be audited against its sources and
so that the regression tests in `tests/testthat/test_vuongtest.R` have a
traceable provenance. Equation numbers refer to Vuong (1989).

## Primary sources

- **Vuong, Q. H. (1989).** Likelihood ratio tests for model selection and
  non-nested hypotheses. *Econometrica, 57*(2), 307–333.
  DOI: <https://doi.org/10.2307/1912557>.
  *(Econometric Society / JSTOR — copyrighted, not open access; cited by DOI in
  place of a bundled PDF.)*
- **Merkle, E. C., You, D., & Preacher, K. J. (2016).** Testing non-nested
  structural equation models. *Psychological Methods, 21*(2), 151–163.
  DOI: <https://doi.org/10.1037/met0000038>.
  Open-access author copy (arXiv): <https://arxiv.org/abs/1402.6720>.

## Equation-by-equation mapping to the code

| Quantity | Source equation | Code location | Formula as implemented |
| --- | --- | --- | --- |
| Individual LR contribution `m_i = log(f1_i / f2_i)` | Vuong §5–6 | `llcont()`, `llA - llB` in `vuongtest()` | `m <- llA - llB` |
| Variance estimator `omega_hat^2` | Vuong Eq (4.2) | `vuongtest()`, `icci()` | `(n-1)/n * var(llA - llB)` (ML variance, denominator `n`) |
| Likelihood-ratio statistic `LR_n` | Vuong Eq (6.4) | `vuongtest()` | `lr <- sum(llA - llB)` |
| AIC correction of `LR_n` | Vuong Eq (5.7), `K = p - q` | `vuongtest(adj="aic")` | `lr - (nparA - nparB)` |
| BIC (Schwarz) correction of `LR_n` | Vuong Eq (5.7), `K = (p - q)·log(n)/2` | `vuongtest(adj="bic")` | `lr - (nparA - nparB) * log(n)/2` |
| Non-nested test statistic (asymptotically `N(0,1)`) | Vuong Thm 5.1 | `vuongtest()` | `(1/sqrt(n)) * lr / sqrt(omega_hat^2)` |
| Nested/overlapping LR statistic (weighted `chi^2`) | Vuong Eq (3.6), Thm 3.3 | `vuongtest(nested=TRUE)` | `2 * lr`, tail via `imhof(2*lr, -lambda)` |
| Distinguishability ("variance") test (weighted `chi^2`, weights `lambda^2`) | Vuong Thm 4.1 | `vuongtest()` | `imhof(n * omega_hat^2, lambda^2)` |
| Matrices `A`, `B` | Vuong Eq (2.1), (2.2) | `calcAB()` | `A = (n·vcov)^-1`, `B = crossprod(scores)/n` |
| Cross-product `B_{f,g}` | Vuong Eq (2.7) | `calcBcross()` | `crossprod(sc1, sc2)/n` |
| Matrix `W` (eigenvalues = weights `lambda`) | Vuong Eq (3.6) | `calcLambda()` | block matrix `[[-B1 A1^-1, -Bc A2^-1], [Bc' A1^-1, B2 A2^-1]]` |
| AIC/BIC difference confidence interval | Merkle, You & Preacher (2016), Eq (7)–(8) | `icci()` | `diff ± qnorm(alpha/2, 1-alpha/2) · sqrt(n·4·omega_hat^2)` |

## Audit status (2026-07)

All of the above were verified against an independent re-derivation from the
source equations on a base-R `glm` fixture; package output matched to machine
precision. No formula damage was found. Recent commits touching the numeric
files were:

- `c9ef3d0` — adds `"DiscreteClass"` to the mirt class dispatch lists in
  `vuongtest()`/`icci()` (routing only; no formula change). **Verified equivalent.**
- `908cac8` — counts free parameters for lavaan models with equality
  constraints via `attr(logLik(object), "df")` instead of `length(coef(object))`
  in the AIC/BIC correction. This makes `nparA`/`nparB` the true free-parameter
  counts `p`/`q` of Vuong Eq (5.7); it is a correctness fix, not damage.
  **Verified correct.**

`tests/testthat/test_vuongtest.R` pins each equation above to both an
independent recomputation and hard-coded reference anchors, so any future
refactor that alters a margin, drops a term, changes a constant, or reorders a
`log`/`sum` in a non-equivalent way will fail the suite.
