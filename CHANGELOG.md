# Changelog

All notable changes to nwalaam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `nwalaam bootstrap, reps(#)`: parametric bootstrap after a fit. Draws outcome vectors from the fitted model (as `nwalaam simulate`), refits each with the same method and effects from the fitted coefficients, and reports model SE, bootstrap SD, bias and percentile interval; failed refits are counted, `e()` and the data are untouched; `r(b_boot)`, `r(sd_boot)`, `r(V_boot)`, `r(ci_boot)`, `r(bias_boot)`, `r(nfail)`; `saving()` writes the refits to a dataset. New test `cscripts/test_nwalaam_bootstrap.do`.

### Help
- States the limits of `method(ee)` and `mple`: the equilibrium-expectation estimates and standard errors on very small networks and its standard errors on a large one, the pseudo-likelihood estimate (standard errors too small, no estimate under separation, r(430)), the advice to report stochastic approximation, and and that the package has no option for the exact likelihood.

## [1.3.1] - 2026-11-20

### Fixed
- A missing value in a covariate of an effect (`ego(x)`, `alter(x)` and the other covariate effects) was silently counted as 0 in the statistics. The fit now stops with r(416), naming the variable and the number of missing values. Fits whose covariates have no missing values are unchanged.
- The stochastic approximation stopped with r(430) ("diverged") for some seeds although other seeds of the same model converge (the book's smoking model with `ego(alc_c) ego(sport_c)`, seed 111). Cause: the start of Algorithm S can contain a huge negative coefficient (-3e25 for a centred covariate whose statistic is nearly 0), because the plugin limits only the upward step; every repetition then began from that start and diverged. The command now restarts from the same start limited to -10 to 10 and, if needed, from zero; `e(b_start)` shows the start used. The estimates of every seed that converged before are unchanged (checked for seeds 1 to 30 of that model, coefficient by coefficient to 15 digits). The error message for a fit that still diverges advises another `seed()`.

### Changed
- Self-loops of a one-mode network were ignored without a word; the command now prints a one-line note with their number (the statistics are unchanged).
- Help: states the missing-covariate rule, the self-loop rule and the starting values; the comparison of `method(ee)` with the exact maximum likelihood estimate is given as measured in checks on small networks (estimates up to 0.5 standard errors from the exact ones, standard errors 0.7 to 1.2 times the exact ones).

The compiled plugin is the one of 1.2.0 (unchanged). The manual keeps its version.

Tested on macOS with Stata 17 BE: the tests of the package including the new `test_nwalaam_patch131` (nwalaam has no dialog test). On Linux (Ubuntu 22.04) and Windows 11 (both StataNow 19.5 MP): all test files pass; on Linux the chapter ego model was also compared for seeds 1 to 30 against 1.3.0: `e(b)` identical to 15 digits, seed 111 now works. Stata 18 not run.

## [1.3.0] - 2026-11-20

### Added
- `estat gof, distplot` draws the histogram of the simulated values of every statistic of the model (`all`: every statistic of the table) with the observed value as a red line. ALAAMEE reports the t-ratios only; this is an extension.
- `saving()`, `export()` and graph options for the plots of `estat gof` and `estat mcmcdiag`; the plots have the default names `nwalaam_gof`, `nwalaam_gofdist` and `nwalaam_trace`. `name()`, `saving()`, `export()` or graph options without a plot option are an error (r(198)).

### Changed
- The plots of `estat gof` and `estat mcmcdiag` use the shared graph helper of the core (`_nw_gfinish`): the plots need nwcommands (core) 3.29.0 or later and stop with r(199) and a message when the core is older; the other commands do not need it. The pictures drawn by the default calls are unchanged (the goodness-of-fit bar plot is identical to 1.2.0 byte for byte; the trace plot differs only by the fix below). `saving()` and `name()` accept a trailing `, replace`; the titles of `distplot` are black.

### Fixed
- `estat mcmcdiag, plot` drew the wrong columns of the stored trace in the panels of the differences dzA (after `method(sa)` the panels were empty, after `method(ee)` they repeated the coefficient trace).

Tested on macOS with Stata 17 BE (plot test, ado test and threshold test, against core 3.29.0); on Linux and on Windows with StataNow 19.5 MP (the plot test and the package tests of this version). Stata 18 on the Mac was not run for this version. The compiled plugin is the one of 1.2.0 (unchanged).

## [1.2.0] - 2026-11-19

### Added
- Three effects for exposure thresholds, for undirected networks and directed networks: `thresholdcontagion(k)`, the number of actors with the outcome that have at least `k` neighbours with the outcome (default `k` = 2); `fractioncontagion(q)`, the number of actors with the outcome of which at least the share `q` of the neighbours have the outcome (0 < `q` <= 1, default 0.5); and `sharecontagion`, the sum over the actors with the outcome of the share of their neighbours that have the outcome. In a directed network the neighbours are the ties in and out (a mutual tie counts twice). With them a specification can represent complex contagion; ALAAMEE has no such effects. They are not defined for two-mode networks. The change statistics are exact (the threshold effects are step functions) and are tested against the definitions of the statistics.

- Three effects for the actors that are exposed but do not have the outcome: `thresholdexposure(k)`, the number of actors without the outcome that have at least `k` neighbours with the outcome (default 2); `fractionexposure(q)`, the number of actors without the outcome of which at least the share `q` of the neighbours have the outcome (default 0.5); and `shareexposure`, the sum over the actors without the outcome of the share of their neighbours that have the outcome. With these and the three effects above the package has 33 effects for undirected and 43 for directed networks. Same neighbour definition as above (a mutual tie counts twice); not defined for two-mode networks. The change statistics are exact and tested against the definitions.

### Changed
- The plugin needs the new version (protocol 7); a plugin of version 1.1.x is refused with a message to update.

Tested on macOS with Stata 17 BE: the change-statistic test (every effect against its definition, including the six new effects), the threshold and exposure test and the main test file, run after the last change; the full suite has not been rerun since the exposure effects were added. The same three tests passed on Linux (Ubuntu 22.04, StataNow 19 MP, plugin built there with `make unix`). Windows: not tested for this version.

## [1.1.1] - 2026-11-18

### Fixed
- Two-mode networks whose nodes are not stored with the first level first (for example the Davis southern women network, whose nodes are in label order): the actors of the outcome level are now found from the modes of the nodes. Before, the first nodes of the network were taken as the outcome level, so the fit used the wrong actors (the header showed too few actors) and could stop with r(430). Networks stored with the first level first are unchanged.

### Changed
- Help file: the equilibrium-expectation standard errors can be too small on larger networks (stochastic approximation's agree with the exact ones); check them against a parametric bootstrap. No change to any command or result.

Tested on macOS with Stata 17 BE (the full test suite). Windows and Linux: not tested for this version; the package has no native or Mata change since 1.1.0.

## [1.1.0] - 2026-11-18

### Added
- `threads(#)`: the independent runs of `method(ee)` are computed by several threads at the same time (macOS, Linux and Windows). The default is one thread per core, at most one per run. The estimates, standard errors and diagnostics are exactly the same for every number of threads and equal those of version 1.0 for the same seed; `test_nwalaam_threads.do` checks this.

### Changed
- The sampler skips the exponential function when a proposal is accepted with certainty (same draws, same results, a little faster).
- Plugin protocol 6.

## [1.0.0] - 2026-11-18

First public release. Autologistic actor attribute models for a binary actor outcome on a fixed network (undirected, directed and two-mode): the effects of ALAAMEE (Stivala), stochastic approximation as the default estimator (checked against the exact maximum likelihood on the Karate Club and on a 16-actor network), the equilibrium-expectation update of ALAAMEE as `method(ee)`, maximum pseudolikelihood (`mple`), `simulate`, `estat gof`, `estat mcmcdiag`, the zero-or-one and snowball samplers, a native plugin for macOS, Linux and Windows. The earlier entries below record the development steps.

## [0.3.1] - 2026-11-18

### Fixed
- StataNow's console binary on a Mac reports `c(os)` = "Unix"; the plugin was then looked up as the Linux file. The platform is now taken as the core does: Unix with a machine type starting "Mac" is macOS. Test `test_nwalaam_platform.do`.

## [0.3.0] - 2026-11-18

### Added
- Two-mode (bipartite) networks: `mode(1|2)` names the level on which the outcome is defined; the effects of ALAAMEE for two-mode networks, including the new `fourcycle1` and `fourcycle2`.
- `settinghomophily(netname)`, `geographichomophily(matname)` and `contagiondist(matname)` (a second network, a distance matrix).
- The zero-or-one sampler (`sampler(zoo)`) and the snowball-conditional sampler (`zone(varname)`: the outermost zone keeps its observed outcome), in estimation, `simulate` and `estat gof`.
- Tests: change statistics against ALAAMEE's for six random networks (undirected, directed, two-mode; missing values in outcome and covariates), `test_nwalaam_alaamee.do`; the exact maximum likelihood (by enumeration) for the zero-or-one sampler with setting and distance effects, a two-mode network and the snowball sampler, `test_nwalaam_exact2.do`.

### Changed
- The directed `gwcontagion`, `logcontagion` and `powercontagion` use ALAAMEE's change statistic (the in-neighbours of an out-neighbour, the out-neighbours of an in-neighbour).
- The covariance matrix of the stochastic-approximation estimate is the inverse of the sample covariance of the simulated statistics (before: of their batch-means covariance); on the 16-actor test network the standard errors are within 2 percent of the exact ones (before: 15 to 20 percent Monte-Carlo error).
- Plugin protocol 5.

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
