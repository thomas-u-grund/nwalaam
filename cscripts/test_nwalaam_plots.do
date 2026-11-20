* nwalaam result plots: estat gof (plot, distplot) and estat mcmcdiag (plot) with name(), saving(), export() and graph options; error paths.
* Run from the repository root (cscripts/run_tests.sh). Writes its graph files to the scratch directory in c(tmpdir).
set seed 123
local tmp `"`c(tmpdir)'"'
local out `"`tmp'/nwalaam_plots_test"'
capture mkdir `"`out'"'
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
gen double cc = age - 30
nwalaam faction karate, effects(density contagion ego(cc)) saiter(100) seed(9)

* ---- estat gof, plot
capture erase `"`out'/gof1.gph"'
capture erase `"`out'/gof1.png"'
estat gof, nsim(100) seed(3) plot name(tg1) saving(`"`out'/gof1.gph"') export(`"`out'/gof1.png"', width(900)) scheme(s2color) xtitle("t")
graph describe tg1
confirm file `"`out'/gof1.gph"'
confirm file `"`out'/gof1.png"'
di as text "  PASS: gof plot, name(), saving(), export() and twoway options"

* ---- estat gof, distplot (in-model statistics; all adds the extra statistics)
capture erase `"`out'/gd1.png"'
estat gof, nsim(100) seed(3) distplot name(tgd1) export(`"`out'/gd1.png"', width(1200)) title("simulated")
graph describe tgd1
confirm file `"`out'/gd1.png"'
estat gof, nsim(100) seed(3) distplot all name(tgd2) saving(`"`out'/gd2.gph"')
graph describe tgd2
confirm file `"`out'/gd2.gph"'
* default names
estat gof, nsim(100) seed(3) distplot
graph describe nwalaam_gofdist
estat gof, nsim(100) seed(3) plot
graph describe nwalaam_gof
* both plots
estat gof, nsim(100) seed(3) plot distplot
di as text "  PASS: gof distplot"

* ---- estat mcmcdiag, plot (stochastic approximation and EE)
capture erase `"`out'/mc1.png"'
estat mcmcdiag, plot name(tm1) saving(`"`out'/mc1.gph"') export(`"`out'/mc1.png"') scheme(s2color)
graph describe tm1
confirm file `"`out'/mc1.gph"'
confirm file `"`out'/mc1.png"'
estat mcmcdiag, plot
graph describe nwalaam_trace
nwalaam faction karate, effects(density contagion) method(ee) runs(3) iterations(4000) thin(50) burnin(500) seed(5)
estat mcmcdiag, plot name(tm2) export(`"`out'/mc2.png"') title("EE trace")
graph describe tm2
confirm file `"`out'/mc2.png"'
di as text "  PASS: mcmcdiag plot"

* ---- error paths
capture estat gof, nsim(100) name(x)
assert _rc == 198
capture estat gof, nsim(100) saving(`"`out'/x.gph"')
assert _rc == 198
capture estat gof, nsim(100) export(`"`out'/x.png"')
assert _rc == 198
capture estat gof, nsim(100) plot distplot name(x)
assert _rc == 198
capture estat mcmcdiag, name(x)
assert _rc == 198
capture estat mcmcdiag, export(`"`out'/x.png"')
assert _rc == 198
capture estat gof, nsim(100) plot nosuchoption(1)
assert _rc != 0
nwalaam faction karate, effects(density contagion) mple
capture estat mcmcdiag, plot
assert _rc == 498
di as text "  PASS: error paths"
di as text "ALL TESTS PASSED"
