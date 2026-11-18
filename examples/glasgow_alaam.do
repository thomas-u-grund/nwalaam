* Smoking in the Glasgow friendship network (Teenage Friends and Lifestyle Study, wave 1; 50 pupils, directed friendship ties):
* who smokes, given the network? Run in a fresh Stata with nwcommands and nwalaam installed:  do examples/glasgow_alaam.do
* A cross-sectional ALAAM describes the association between a pupil's smoking and the smoking of the friends; it cannot tell
* whether friends influence each other or whether pupils pick friends who smoke (selection).
clear all
set linesize 80
nwwebuse glasgow, nwclear
generate byte smoker1 = (smoke1 > 1)          // occasional or regular smoker
tabulate smoker1

* 1. The pseudo-likelihood estimate: a logit of the outcome on the change statistics (fast, but standard errors too small)
nwalaam smoker1 glasgow1, effects(density contagion) mple

* 2. The default estimator (stochastic approximation): the maximum likelihood estimate
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) seed(20261118)
estat mcmcdiag

* 3. Covariates: alcohol use and sport; the covariates are not centred by the command, so centre them
quietly summarize alcohol1
generate double alc_c = alcohol1 - r(mean)
quietly summarize sport1
generate double sport_c = sport1 - r(mean)
nwalaam smoker1 glasgow1, effects(density contagion sender receiver ego(alc_c) ego(sport_c)) seed(20261118)

* 4. Goodness of fit: does the model reproduce statistics that are not in it (reciprocity, triangles, two-stars)?
estat gof, nsim(1000) seed(7)

* 5. The method of ALAAMEE (equilibrium expectation) for comparison
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) method(ee) runs(30) seed(20261118)

* 6. Relation to the exposure regression of the book (Chapter 11)
* The contagion term is the number of ties to smokers (a mutual tie counts twice). The pseudo-likelihood with density and
* contagion alone is a logit of smoking on that count; the exposure regression uses the share of friends who smoke instead.
nwaltergen frsmoke = proportion(alter.smoker1==1), net(glasgow1)
regress smoker1 frsmoke
nwaltergen nsmokefr = sum(alter.smoker1), net(glasgow1)
logit smoker1 nsmokefr
