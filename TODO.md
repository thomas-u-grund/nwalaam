# TODO (nwalaam 0.3.0)

Done: undirected (30) and directed (40) effects in the plugin, checked against their definitions; maximum pseudo-likelihood (= logit);
stochastic approximation (default; agrees with the exact maximum likelihood) and the EE update of ALAAMEE; `nwalaam simulate`, `estat gof`,
`estat mcmcdiag`; tests (see `RELEASING.md`); examples on Glasgow and Lazega.

Done in 0.3.0: two-mode networks and effects, setting / geographic / contagion-distance effects, zero-or-one and snowball-conditional samplers, change statistics tested against ALAAMEE, better SA standard errors.

Open
1. Linux and Windows plugins: build with `.github/workflows/build-plugins.yml` once the repository has a remote; test on Linux and Windows (not tested at all so far).
2. Integration into the core: `nwinstall` registration (`modelpkg`), `nwversion`, `nwtopical`, dialog, website pages; a core release.
3. Check ALAAMEE's zero-or-one acceptance ratio against the exact distribution (by reading, it uses Dy / (Dmax - Dy) with Dy counted after the toggle, which differs from the exact ratio by one actor); nwalaam uses the exact ratio.
4. ALAAMEE's outcome-start and simulation conventions in `estat gof` (it starts from random outcomes; nwalaam from the observed ones); ALAAMEE's input and output files and R scripts are not reproduced (deliberate: Stata datasets and `e()` instead).
5. Deliberate non-matches: continuous covariates are not centred (as ALAAMEE); the MPNet change-statistics file of ALAAMEE's `tiny` bipartite example lists C4X-1 values twice those of ALAAMEE's Python (we follow the Python, which agrees with the MPNet observed statistic).
6. EE standard errors on larger networks: on ALAAMEE's simulated 500-node example (generating values -7.2, 0.55, 1.0, 1.2, 1.15) a parametric bootstrap (40 outcomes simulated from the SA fit, refitted by SA) gives standard deviations 0.70, 0.110, 0.178, 0.314, 0.210; the SA standard errors are 0.82, 0.117, 0.165, 0.326, 0.189 (within 20 percent) but the EE ones are 0.43, 0.038, 0.17, 0.27, 0.21, too small for density and activity (estimates: SA -7.94, 0.67, 0.98, 1.34, 1.11; EE -7.60, 0.61, 1.02, 1.28, 1.08). Look at the EE information part (autocorrelation of the dzA batch means) and at advising SA for standard errors.
