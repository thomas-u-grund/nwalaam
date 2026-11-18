# TODO (nwalaam 0.2.0)

Done: undirected (27) and directed (37) effects in the plugin, checked against their definitions; maximum pseudo-likelihood (= logit);
stochastic approximation (default; agrees with the exact maximum likelihood) and the EE update of ALAAMEE; `nwalaam simulate`, `estat gof`,
`estat mcmcdiag`; tests (see `RELEASING.md`); examples on Glasgow and Lazega.

Open
1. Two-mode (bipartite) networks and their effects (ALAAMEE `changeStatisticsALAAMbipartite.py`); test data: ALAAMEE's `tiny` bipartite example with the MPNet change statistics file.
2. Setting homophily, geographic homophily and contagion by distance (need a second network or a distance matrix).
3. The zoo / conditional sampler of ALAAMEE for rare outcomes (the basic toggle sampler wastes proposals when few nodes have the outcome).
4. Snowball-sampled data (zones) and the structurally fixed outcome of ALAAMEE.
5. Linux and Windows plugins: build with `.github/workflows/build-plugins.yml` once the repository has a remote; test on Linux and Windows (not tested at all so far).
6. Integration into the core: `nwinstall` registration (`modelpkg`), `nwversion`, `nwtopical`, dialog, website pages; a core release.
7. Comparison with the ALAAMEE Python run (it was run once on the karate club by the coordinator: both EE implementations agree); a systematic comparison needs the owner's go.
8. Stochastic approximation: the standard errors have about 15 percent Monte-Carlo error at the default `saiter(500)`; a better estimator of the information matrix (or more phases) would steady them.
