* Default calls (no name, saving, export or graph options) of the plots of estat gof and estat mcmcdiag; each default picture is exported to
* the directory in the global NWP_OUT. Run with the released nwalaam (1.2.0) and with the changed one: the PNG files must be identical,
* except the traces (trace_sa, trace_ee), whose dzA panels were drawn from the wrong columns in 1.2.0 and are fixed.
* Used by test_nwalaam_plots.do (reference directory: environment variable NWALAAM_PLOT_REF). Produce the reference after `cscript`.
set graphics off
capture mkdir "$NWP_OUT"
set seed 123
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
gen double cc = age - 30
nwalaam faction karate, effects(density contagion ego(cc)) saiter(100) seed(9)
estat gof, nsim(100) seed(3) plot
quietly graph export "$NWP_OUT/gof.png", replace width(600)
graph drop _all
estat mcmcdiag, plot
quietly graph export "$NWP_OUT/trace_sa.png", replace width(600)
graph drop _all
nwalaam faction karate, effects(density contagion) method(ee) runs(3) iterations(4000) thin(50) burnin(500) seed(5)
estat mcmcdiag, plot
quietly graph export "$NWP_OUT/trace_ee.png", replace width(600)
graph drop _all
