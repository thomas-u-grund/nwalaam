* Smoking in the Glasgow friendship network (Teenage Friends and Lifestyle Study, wave 1; 50 pupils, directed friendship ties)
clear all
set linesize 90
log using out/glasgow_data.txt, text replace
nwwebuse glasgow, nwclear
generate byte smoker1 = (smoke1 > 1)
tabulate smoker1
log close
log using out/glasgow_mple.txt, text replace
nwalaam smoker1 glasgow1, effects(density contagion) mple
log close
log using out/glasgow_sa.txt, text replace
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) seed(20261118)
log close
log using out/glasgow_cov.txt, text replace
quietly summarize alcohol1
generate double alc_c = alcohol1 - r(mean)
quietly summarize sport1
generate double sport_c = sport1 - r(mean)
nwalaam smoker1 glasgow1, effects(density contagion sender receiver ego(alc_c) ego(sport_c)) seed(20261118)
log close
log using out/glasgow_gof.txt, text replace
estat gof, nsim(1000) seed(7)
log close
log using out/glasgow_ee.txt, text replace
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) method(ee) runs(30) seed(20261118)
log close
log using out/glasgow_replay.txt, text replace
nwalaam, ci
nwalaam, noci level(90)
log close
