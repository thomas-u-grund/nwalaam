* Zachary's karate club, faction as outcome (34 members); edges and attributes from the ALAAMEE repository (MIT licence)
clear all
set linesize 90
set seed 20261118
log using out/karate_data.txt, text replace
import delimited "data/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "data/karate_nodes.csv", clear
tabulate faction
nwsummarize karate
log close
log using out/karate_sa.txt, text replace
nwalaam faction karate, effects(density activity contagion) seed(20261118)
log close
log using out/karate_diag.txt, text replace
estat mcmcdiag
log close
log using out/karate_ee.txt, text replace
nwalaam faction karate, effects(density activity contagion) method(ee) seed(20261118)
estat mcmcdiag
log close
log using out/karate_mple.txt, text replace
capture noisily nwalaam faction karate, effects(density activity contagion) mple
nwalaam faction karate, effects(density contagion) mple
log close
log using out/karate_gof.txt, text replace
nwalaam faction karate, effects(density activity contagion) seed(20261118)
estat gof, nsim(1000) seed(7)
log close
log using out/karate_covs.txt, text replace
quietly summarize age
generate double age_c = age - r(mean)
nwalaam faction karate, effects(density contagion ego(age_c) same(gender)) seed(20261118)
log close
log using out/karate_boot.txt, text replace
nwalaam faction karate, effects(density activity contagion) seed(20261118)
timer clear 1
timer on 1
nwalaam bootstrap, reps(40) seed(11) nodots
timer off 1
timer list 1
log close
