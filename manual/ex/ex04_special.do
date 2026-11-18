* Two-mode networks, setting network, distance matrix, zero-or-one sampler, snowball zones
clear all
set linesize 90
log using out/special_twomode.txt, text replace
nwwebuse davis, nwclear
nwset
generate byte mrs = (title == "Mrs") if _nwmode == "1"
tabulate mrs
nwalaam mrs davis, effects(density activity fourcycle1 fourcycle2) mode(1) seed(20261118)
log close
log using out/special_setting.txt, text replace
nwwebuse lazega, nwclear
nwset
generate byte corporate = (practicenum == 2)
quietly summarize age
generate double age_c = age - r(mean)
nwalaam corporate lazega_fr, effects(density contagion ego(age_c) settinghomophily(lazega_cow)) seed(20261118)
log close
log using out/special_distance.txt, text replace
mata: ag = st_data(., "age"); st_matrix("DIST", abs(ag * J(1, rows(ag), 1) - J(rows(ag), 1, 1) * ag'))
nwalaam corporate lazega_fr, effects(density contagion geographichomophily(DIST)) seed(20261118)
log close
log using out/special_zoo.txt, text replace
nwalaam corporate lazega_fr, effects(density contagion) sampler(zoo) seed(20261118)
log close
log using out/special_zone.txt, text replace
generate zone = (_n > 10) + (_n > 40)
tabulate zone
nwalaam corporate lazega_fr, effects(density contagion) zone(zone) seed(20261118)
log close
