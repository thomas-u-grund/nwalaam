# TODO (nwalaam 0.1.0)

Done: undirected (21) and directed (25) effects in the plugin, checked against their definitions; maximum pseudo-likelihood (= logit);
EE estimator with Algorithm S start, several runs, batch-means information matrix, convergence t-ratios; display rules of the family;
tests (see `RELEASING.md`).

Open
1. Two-mode (bipartite) networks and their effects (ALAAMEE `changeStatisticsALAAMbipartite.py`); test data: ALAAMEE's `tiny` bipartite example with the MPNet change statistics file.
2. More effects: setting homophily, geographic homophily and contagion by distance (need a second network or a distance matrix), same/diff partner-activity two-paths, binary-covariate alter two-stars (`alterbinarytwostar1/2`), directed `gwcontagion`, directed `same`/`diff`/`samecontagion`.
3. `nwalaam simulate` (draw outcome vectors from the fitted model) and `estat gof` (the goodness-of-fit t-ratios of ALAAMEE), `estat` for the convergence diagnostics (plots of theta and dzA).
4. The zoo / conditional sampler of ALAAMEE for rare outcomes (the basic toggle sampler wastes proposals when few nodes have the outcome).
5. Snowball-sampled data (zones) and the structurally fixed outcome of ALAAMEE; stochastic approximation (`estimateALAAMSA`) as a second estimator.
6. Linux and Windows plugins: build with `.github/workflows/build-plugins.yml` once the repository has a remote; test on Linux and Windows with Stata 18/19 (not tested at all so far).
7. Test under Stata 18 and 19 (only 17 so far).
8. Integration into the core: `nwinstall` registration (`modelpkg`), `nwversion`, website pages; a core release.
9. Comparison with the ALAAMEE Python run: not possible in the first session (the sandbox refused to execute the downloaded repository); needs the owner's go.
10. The EE point estimate differs from the exact maximum likelihood estimate by about 0.2 exact standard errors on the 16-node test network, with a learning rate that does not change this (see `dev/lr_check.do`); understand whether this is a property of the algorithm (the same in ALAAMEE) or a bias to remove.
