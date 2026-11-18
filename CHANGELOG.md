# Changelog

All notable changes to nwalaam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-11-18

### Changed
- **The default estimator is now stochastic approximation (`method(sa)`)**; the equilibrium-expectation update of ALAAMEE is `method(ee)`. Reason: compared with the exact maximum likelihood estimate (computed by enumerating all outcome vectors), the EE update is biased. On the karate club (34 actors, 2^34 outcome vectors enumerated: Density 0.535, Activity -0.832, Contagion 1.461) stochastic approximation gives 0.52 / -0.83 / 1.46 and EE gives -0.03 / -0.75 / 1.53 (as does ALAAMEE's own Python implementation); on a random 16-actor network the EE estimates lie 0.15 to 0.24 exact standard errors below the exact estimate, stochastic approximation 0.00 to 0.01. The standard errors of both agree with the exact ones (within the Monte-Carlo error).
- Stochastic approximation is Snijders' Robbins-Monro procedure with a derivative matrix from the simulated statistics; it repeats (up to 20 times) until every convergence t-ratio is below 0.1. The covariance matrix is the inverse of the batch-means covariance of the simulated statistics at the estimate. `saiter(#)` sets its length.

### Added
- `nwalaam simulate`: draws outcome vectors from the fitted model or from given coefficients (`theta()`), stored as variables (`generate()`) or in a frame (`frame()`), with `seed()`.
- `estat gof`: simulation-based goodness of fit (t-ratios and the Mahalanobis distance of the observed statistics from the simulated ones) for the statistics of the model and a default set of other statistics (as in ALAAMEE), with `extra()`, `noextra`, `nsim()`, `plot`.
- `estat mcmcdiag`: convergence diagnostics per run (estimate, convergence t-ratios, autocorrelation) and trace plots.
- Effects: `samepartneractivity`, `diffpartneractivity`, `sameindirect`, `diffindirect`, `alterbinarytwostar1`, `alterbinarytwostar2` (undirected); `gwcontagion`, `logcontagion`, `powercontagion`, `same`, `diff`, `samecontagion`, `sendermatch`, `sendermismatch`, `receivermatch`, `receivermismatch`, `reciprocitymatch`, `reciprocitymismatch` (directed). In directed networks a mutual tie counts twice in the contagion statistics, as in `contagion`.
- Examples on the Glasgow and Lazega data in `examples/`.
- The EE traces of the first run and per-run diagnostics are stored (`e(trace1)`, `e(run_diag)`).

## [0.1.0] - 2026-11-18

### Added
- First version. `nwalaam outcome netname, effects(...)` fits an autologistic actor attribute model (ALAAM) for a binary node outcome on a fixed undirected or directed one-mode network, with the equilibrium-expectation (EE) estimator of ALAAMEE (several independent runs, Algorithm S starting values, batch-means covariance, convergence t-ratios) and a maximum pseudo-likelihood estimate (`mple`).
- 21 undirected and 25 directed effects (see `help nwalaam`), all in the native plugin. Not yet: two-mode networks, setting and geographic homophily, distance effects, simulation and goodness of fit.
