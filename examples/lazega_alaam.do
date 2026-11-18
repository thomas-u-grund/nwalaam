* Corporate or litigation practice among the 71 lawyers of Lazega's law firm, given the friendship network (undirected):
* are friends in the same kind of practice? Run in a fresh Stata with nwcommands and nwalaam installed:
* do examples/lazega_alaam.do
clear all
set linesize 80
nwwebuse lazega, nwclear
generate byte corporate = (practicenum == 2)
quietly summarize age
generate double age_c = age - r(mean)
tabulate corporate
nwalaam corporate lazega_fr, effects(density contagion activity ego(age_c)) seed(20261118)
estat mcmcdiag
estat gof, nsim(1000) seed(7)
* a smaller model and the equilibrium-expectation estimate of ALAAMEE for comparison
nwalaam corporate lazega_fr, effects(density contagion) method(ee) runs(30) seed(20261118)
