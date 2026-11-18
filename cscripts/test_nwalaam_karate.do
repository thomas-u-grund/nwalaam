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
// the exact maximum likelihood estimate, from the enumeration of all 2^34 outcome vectors (histogram of the three statistics, then
// Newton's method): Density 0.535 (se 0.941), Activity -0.832 (0.295), Contagion 1.461 (0.394)
matrix EXACT = (0.535, -0.832, 1.461)
matrix EXSE = (0.941, 0.295, 0.394)
// 1. the default (stochastic approximation) agrees with the exact MLE: mean over 4 fits within 0.1 exact se, each fit within 0.4
matrix ACC = J(1, 3, 0)
forvalues r = 1/4 {
	nwalaam faction karate, effects(density activity contagion) seed(`r'3)
	assert "`e(method)'" == "sa"
	matrix B = e(b)
	matrix SE = vecdiag(e(V))
	forvalues j = 1/3 {
		di "SA fit `r', effect `j': " %8.4f B[1, `j'] " (se " %6.4f sqrt(SE[1, `j']) ")   exact " %8.4f EXACT[1, `j'] " (se " %6.4f EXSE[1, `j'] ")"
		assert abs(B[1, `j'] - EXACT[1, `j']) < 0.4 * EXSE[1, `j']
		assert inrange(sqrt(SE[1, `j']) / EXSE[1, `j'], 0.6, 1.4)
		matrix ACC[1, `j'] = ACC[1, `j'] + (B[1, `j'] - EXACT[1, `j']) / EXSE[1, `j'] / 4
	}
}
forvalues j = 1/3 {
	assert abs(ACC[1, `j']) < 0.1
}
di "PASS: stochastic approximation agrees with the exact maximum likelihood of the karate club"
// 2. method(ee), the ALAAMEE update: agrees with MPNet (published with ALAAMEE) within 0.7 MPNet se, but is biased against the
// exact MLE (Density about 0.6 exact se lower): the test only records that it stays within 1 exact se
nwalaam faction karate, effects(density activity contagion) method(ee) runs(40) iterations(100000) seed(20261118)
matrix B = e(b)
matrix SE = vecdiag(e(V))
matrix MPNET = (0.3809, -0.7975, 1.4523)
matrix MPSE = (0.906, 0.292, 0.400)
forvalues j = 1/3 {
	di "effect `j': EE " %8.4f B[1, `j'] " (se " %6.4f sqrt(SE[1, `j']) ")   MPNet " %8.4f MPNET[1, `j'] " (se " %6.4f MPSE[1, `j'] ")   exact " %8.4f EXACT[1, `j']
	assert abs(B[1, `j'] - MPNET[1, `j']) < 0.7 * MPSE[1, `j']
	assert inrange(sqrt(SE[1, `j']) / MPSE[1, `j'], 0.75, 1.3)
	assert abs(B[1, `j'] - EXACT[1, `j']) < 1 * EXSE[1, `j']
}
assert e(runs_used) == 40
di "PASS: method(ee) estimates agree with MPNet"
