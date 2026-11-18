* Corporate or litigation practice among the 71 lawyers of Lazega's law firm, given the friendship network (directed)
clear all
set linesize 90
log using out/lazega_data.txt, text replace
nwwebuse lazega, nwclear
nwset
generate byte corporate = (practicenum == 2)
quietly summarize age
generate double age_c = age - r(mean)
tabulate corporate
log close
log using out/lazega_sa.txt, text replace
nwalaam corporate lazega_fr, effects(density contagion sender receiver ego(age_c)) seed(20261118)
log close
log using out/lazega_diag.txt, text replace
estat mcmcdiag
estat gof, nsim(1000) seed(7)
log close
log using out/lazega_ee.txt, text replace
nwalaam corporate lazega_fr, effects(density contagion) method(ee) runs(30) seed(20261118)
log close
log using out/lazega_threads.txt, text replace
timer clear
timer on 1
nwalaam corporate lazega_fr, effects(density contagion) method(ee) runs(30) seed(1) threads(1)
timer off 1
matrix B1 = e(b)
timer on 2
nwalaam corporate lazega_fr, effects(density contagion) method(ee) runs(30) seed(1) threads(8)
timer off 2
matrix B8 = e(b)
display "largest difference between the two coefficient vectors: " mreldif(B1, B8)
timer list
log close
log using out/lazega_sim.txt, text replace
nwalaam corporate lazega_fr, effects(density contagion sender receiver ego(age_c)) seed(20261118)
nwalaam simulate, nsim(200) frame(simfr) seed(3)
frame simfr: collapse (sum) sim*
frame simfr: xpose, clear
frame simfr: summarize v1
display "observed number of corporate lawyers: " e(n_ones)
nwalaam simulate, nsim(3) generate(draw) seed(4)
summarize draw1 draw2 draw3
log close
