* Plots of estat gof and estat mcmcdiag on the karate club (figures into ../fig)
clear all
set linesize 90
set seed 20261118
set scheme s1network
import delimited "data/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "data/karate_nodes.csv", clear
quietly nwalaam faction karate, effects(density activity contagion) seed(20261118)
log using out/plots_gof.txt, text replace
estat gof, nsim(1000) seed(7) plot title("") xsize(5) ysize(3) name(karate_gof) export(../fig/fig_gof.pdf)
estat gof, nsim(1000) seed(7) distplot xsize(6.5) ysize(2.4) name(karate_dist) export(../fig/fig_gofdist.pdf)
log close
log using out/plots_err.txt, text replace
capture noisily estat gof, nsim(100) seed(7) export(never.pdf)
display "return code: " _rc
log close
quietly nwalaam faction karate, effects(density activity contagion) method(ee) seed(20261118)
estat mcmcdiag, plot title("") xsize(6.5) ysize(4) name(karate_trace) export(../fig/fig_trace.pdf)
