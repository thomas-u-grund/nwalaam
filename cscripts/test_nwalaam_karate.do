cscript
* Zachary's karate club with the faction as outcome (the example of ALAAMEE; edges and attributes in dev/karate, from the ALAAMEE
* repository, MIT licence): the observed statistics, and the EE estimates against the MPNet estimates published with ALAAMEE
* (Density 0.3809 (se 0.906), Activity -0.7975 (0.292), Contagion 1.4523 (0.400)).
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
_nwsyntax karate
mata: AlaamChangeStats(*(`netobj'->get_matrix()), 0, st_data(., "faction"), (1, 2, 6, 3, 11, 12, 13, 7, 8), J(1, 9, 0), J(1, 9, 0), J(0, 34, .), D = ., z = .)
mata: st_matrix("OBS", z)
matrix list OBS
// Density 18, Activity 80, Contagion 35, Ego2Star 271, Alter2Star1 566, Alter2Star2 219, T1 57, T2 55, T3 18
matrix REF = (18, 80, 35, 271, 57, 55, 18, 566, 219)
assert mreldif(OBS, REF) == 0
di "PASS: the observed statistics of the karate club equal the ALAAMEE values"
nwalaam faction karate, effects(density activity contagion) runs(40) iterations(100000) seed(20261118)
matrix B = e(b)
matrix SE = vecdiag(e(V))
matrix MPNET = (0.3809, -0.7975, 1.4523)
matrix MPSE = (0.906, 0.292, 0.400)
forvalues j = 1/3 {
	di "effect `j': nwalaam " %8.4f B[1, `j'] " (se " %6.4f sqrt(SE[1, `j']) ")   MPNet " %8.4f MPNET[1, `j'] " (se " %6.4f MPSE[1, `j'] ")"
	* the MPNet estimates are themselves Monte-Carlo estimates: agreement within 0.7 standard errors, standard errors within 25 percent
	assert abs(B[1, `j'] - MPNET[1, `j']) < 0.7 * MPSE[1, `j']
	assert inrange(sqrt(SE[1, `j']) / MPSE[1, `j'], 0.75, 1.3)
}
assert e(runs_used) == 40
di "PASS: karate club estimates agree with MPNet"
