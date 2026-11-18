cscript
* The independent EE runs execute in parallel threads. The result must not depend on the number of threads (each run has its own
* generator and memory and its own place in the output), and a bad threads() is refused. Karate club, 9 runs of 4000 iterations.
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
_nwsyntax karate
local eff "density activity contagion"
nwalaam faction karate, effects(`eff') method(ee) runs(9) iterations(4000) thin(50) burnin(500) seed(8641) threads(1)
matrix B1 = e(b)
matrix V1 = e(V)
matrix R1 = e(run_diag)
matrix T1 = e(trace1)
scalar acc1 = e(acceptance)
local used1 = e(runs_used)
foreach th in 2 3 4 9 64 0 {
	nwalaam faction karate, effects(`eff') method(ee) runs(9) iterations(4000) thin(50) burnin(500) seed(8641) threads(`th')
	assert mreldif(e(b), B1) == 0
	assert mreldif(e(V), V1) == 0
	assert mreldif(e(run_diag), R1) == 0
	assert mreldif(e(trace1), T1) == 0
	assert e(acceptance) == acc1
	assert e(runs_used) == `used1'
	di "  threads(`th'): identical to threads(1)"
}
di "PASS: the EE result does not depend on the number of threads"
* with fewer runs than threads, and a single run
nwalaam faction karate, effects(`eff') method(ee) runs(1) iterations(3000) thin(50) burnin(500) seed(5) threads(1)
matrix B1 = e(b)
nwalaam faction karate, effects(`eff') method(ee) runs(1) iterations(3000) thin(50) burnin(500) seed(5) threads(8)
assert mreldif(e(b), B1) == 0
di "PASS: a single run"
* the other estimators ignore threads(), and their results do not change
nwalaam faction karate, effects(`eff') seed(77) threads(1)
matrix B1 = e(b)
nwalaam faction karate, effects(`eff') seed(77) threads(6)
assert mreldif(e(b), B1) == 0
di "PASS: stochastic approximation does not depend on threads()"
* bad values
foreach bad in -1 257 {
	capture nwalaam faction karate, effects(`eff') method(ee) runs(2) iterations(1000) thin(50) threads(`bad')
	assert _rc == 198
}
di "PASS: threads() out of range -> r(198)"
di "test_nwalaam_threads: ALL OK"
