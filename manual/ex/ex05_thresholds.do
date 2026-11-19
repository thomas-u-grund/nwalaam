* Exposure thresholds on the karate club: contagion against thresholdcontagion(2); the fractions and the share
clear all
set linesize 90
set seed 20261118
import delimited "data/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "data/karate_nodes.csv", clear
log using out/thresholds_karate.txt, text replace
nwalaam faction karate, effects(density activity thresholdcontagion(2)) seed(20261118)
nwalaam faction karate, effects(density activity sharecontagion) seed(20261118)
log close
