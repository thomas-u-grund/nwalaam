# nwalaam

<!--version-->**Version 1.3.1** (released 20 November 2026)<!--/version-->

Autologistic actor attribute models (ALAAMs) for Stata: a binary outcome on the nodes of a fixed network, explained by the network
structure around the nodes (contagion, activity, triangles, covariates ...). Part of the nwcommands family; requires the core package
`nwcommands` (3.29.0 or later). Free software, GPL-3. Follows ALAAMEE (Stivala, Wang and Lomi) and MPNet.

```stata
. nwalaam faction karate, effects(density activity contagion) seed(1)     // stochastic approximation (maximum likelihood)
. nwalaam faction karate, effects(density activity contagion) method(ee)  // the equilibrium-expectation update of ALAAMEE
. estat mcmcdiag
. estat gof
. nwalaam simulate, nsim(100) generate(sim) seed(2)
. nwalaam y affiliations, effects(density activity fourcycle1 fourcycle2) mode(1)       // two-mode network
. nwalaam y g, effects(density contagion settinghomophily(school)) sampler(zoo)        // setting network, zero-or-one sampler
. nwalaam faction karate, effects(density activity thresholdcontagion(2))                  // exposure threshold (complex contagion); thresholdexposure(2) counts the exposed actors without the outcome
. nwalaam y g, effects(density contagion) zone(zone)                                   // snowball sample (conditional sampler)
```

The effects, the three samplers (toggle, zero-or-one, snowball-conditional), undirected, directed and two-mode networks and the setting-network and
distance effects of ALAAMEE are available; the change statistics agree with ALAAMEE's to rounding error on random networks of each kind
(`cscripts/test_nwalaam_alaamee.do`).

The default estimator is stochastic approximation: it reproduces the exact maximum likelihood estimate where that can be computed
by enumeration (the karate club), while the equilibrium-expectation update of ALAAMEE stops at slightly different estimates; see
`help nwalaam` and `CHANGELOG.md`. Examples on the Glasgow and Lazega data are in `examples/`.

Status: version 0.3.0, not released; a local repository without a remote. Tested on macOS with Stata 17. Linux and Windows plugins are
built by `.github/workflows/build-plugins.yml` when the repository is published; until then they are not tested. See `help nwalaam`,
`CHANGELOG.md`, `NOTICE.md` and `RELEASING.md`.

## Manual

A user manual (PDF in `manual/`, with the LaTeX source and the do-files of its examples) is archived at Zenodo: <https://doi.org/10.5281/zenodo.24397901> (all versions). Cite the manual by that address together with the package version printed by `nwalaam version`.
