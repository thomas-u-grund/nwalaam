# Changelog

All notable changes to nwalaam are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project uses [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-18

### Added
- First version. `nwalaam outcome netname, effects(...)` fits an autologistic actor attribute model (ALAAM) for a binary node outcome on a fixed undirected or directed one-mode network, with the equilibrium-expectation (EE) estimator of ALAAMEE (several independent runs, Algorithm S starting values, batch-means covariance, convergence t-ratios) and a maximum pseudo-likelihood estimate (`mple`).
- 21 undirected and 25 directed effects (see `help nwalaam`), all in the native plugin. Not yet: two-mode networks, setting and geographic homophily, distance effects, simulation and goodness of fit.
