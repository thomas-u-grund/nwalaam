# nwalaam

<!--version-->**Version 0.1.0** (released 18 November 2026)<!--/version-->

Autologistic actor attribute models (ALAAMs) for Stata: a binary outcome on the nodes of a fixed network, explained by the network
structure around the nodes (contagion, activity, triangles, covariates ...). Part of the nwcommands family; requires the core package
`nwcommands` (3.26.0 or later). Free software, GPL-3. Follows ALAAMEE (Stivala, Wang and Lomi) and MPNet.

```stata
. nwalaam faction karate, effects(density activity contagion) runs(30) seed(1)
```

Status: version 0.1.0, not released; a local repository without a remote. Tested on macOS with Stata 17. Linux and Windows plugins are
built by `.github/workflows/build-plugins.yml` when the repository is published; until then they are not tested. See `help nwalaam`,
`CHANGELOG.md`, `NOTICE.md` and `RELEASING.md`.
