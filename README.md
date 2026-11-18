# nwalaam

<!--version-->**Version 0.2.0** (released 18 November 2026)<!--/version-->

Autologistic actor attribute models (ALAAMs) for Stata: a binary outcome on the nodes of a fixed network, explained by the network
structure around the nodes (contagion, activity, triangles, covariates ...). Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows ALAAMEE (Stivala, Wang and Lomi) and MPNet.

```stata
. nwalaam faction karate, effects(density activity contagion) seed(1)     // stochastic approximation (maximum likelihood)
. nwalaam faction karate, effects(density activity contagion) method(ee)  // the equilibrium-expectation update of ALAAMEE
. estat mcmcdiag
. estat gof
. nwalaam simulate, nsim(100) generate(sim) seed(2)
```

The default estimator is stochastic approximation: it reproduces the exact maximum likelihood estimate where that can be computed
by enumeration (the karate club), while the equilibrium-expectation update of ALAAMEE stops at slightly different estimates; see
`help nwalaam` and `CHANGELOG.md`. Examples on the Glasgow and Lazega data are in `examples/`.

Status: version 0.2.0, not released; a local repository without a remote. Tested on macOS with Stata 17. Linux and Windows plugins are
built by `.github/workflows/build-plugins.yml` when the repository is published; until then they are not tested. See `help nwalaam`,
`CHANGELOG.md`, `NOTICE.md` and `RELEASING.md`.
