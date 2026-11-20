* nwalaam 1.3.1: missing covariate values stop the fit (r(416)), self-loops give a note and change nothing, and a seed whose
* Algorithm S start is wild (seed 111 of the glasgow smoking model) converges from the limited start while the other seeds
* keep their estimates. Needs the glasgow data (the core data folder through NWCORE, else the web).
set seed 123
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
gen double cc = age - 30
gen double ccmiss = cc
replace ccmiss = . in 3
gen y2 = faction
replace y2 = . in 5/6

* ---- 1. missing covariate values: error 416 for every estimator and for every covariate effect, also when the actor has no outcome
capture nwalaam faction karate, effects(density contagion ego(ccmiss)) mple
assert _rc == 416
capture nwalaam faction karate, effects(density contagion ego(ccmiss)) seed(1)
assert _rc == 416
capture nwalaam faction karate, effects(density contagion alter(ccmiss)) method(ee) runs(2) iterations(2000) thin(50) burnin(500) seed(1)
assert _rc == 416
capture nwalaam simulate faction karate, effects(density contagion ego(ccmiss)) theta(nosuch) nsim(10) generate(zz) seed(1)
assert _rc != 0
* the message names the variable and the number
tempfile lg
log using `lg', text replace name(p131)
capture noisily nwalaam faction karate, effects(density contagion ego(ccmiss)) mple
log close p131
file open fh using `lg', read text
local found 0
file read fh line
while r(eof) == 0 {
	if strpos(`"`line'"', "ccmiss") & strpos(`"`line'"', "1 missing") local found 1
	file read fh line
}
file close fh
assert `found' == 1
* no missing values: runs
nwalaam faction karate, effects(density contagion ego(cc)) mple
di as text "PASS: missing covariate values stop the fit with r(416)"

* ---- 2. self-loops: a note, and the same estimate as without the loops
nwalaam faction karate, effects(density contagion ego(cc)) mple
matrix B0 = e(b)
preserve
import delimited "dev/karate/karate_edges.csv", clear
local n1 = _N + 1
local n2 = _N + 2
set obs `n2'
replace from = 1 in `n1'
replace to = 1 in `n1'
replace from = 5 in `n2'
replace to = 5 in `n2'
nwfromedge from to, undirected name(kloop) selfloop
restore
log using `lg', text replace name(p131b)
capture noisily nwalaam faction kloop, effects(density contagion ego(cc)) mple
log close p131b
assert _rc == 0
matrix B1 = e(b)
assert mreldif(B0, B1) < 1e-12
file open fh using `lg', read text
local found 0
local clean 1
file read fh line
while r(eof) == 0 {
	if strpos(`"`line'"', "2 self-loop") local found 1
	file read fh line
}
file close fh
assert `found' == 1
log using `lg', text replace name(p131c)
nwalaam faction karate, effects(density contagion ego(cc)) mple
log close p131c
file open fh using `lg', read text
file read fh line
while r(eof) == 0 {
	assert strpos(`"`line'"', "self-loop") == 0
	file read fh line
}
file close fh
di as text "PASS: self-loops give a note and are ignored"

* ---- 3. seed 111 of the book's smoking model
local core : environment NWCORE
if "`core'" != "" global nwwebpath "`core'/data"
nwclear
nwwebuse glasgow, nwclear
generate byte smoker1 = (smoke1 > 1)
quietly summarize alcohol1
generate double alc_c = alcohol1 - r(mean)
quietly summarize sport1
generate double sport_c = sport1 - r(mean)
local eff "density contagion sender receiver ego(alc_c) ego(sport_c)"
nwalaam smoker1 glasgow1, effects(`eff') seed(111)
assert e(sa_converged) == 1
matrix b111 = e(b)
matrix st111 = e(b_start)
assert max(abs(st111[1,1]), abs(st111[1,2]), abs(st111[1,3]), abs(st111[1,4]), abs(st111[1,5]), abs(st111[1,6])) <= 10
* estimates agree with the other seeds (all within a few Monte-Carlo errors)
nwalaam smoker1 glasgow1, effects(`eff') seed(1)
matrix b1 = e(b)
forvalues j = 1/6 {
	assert abs(b111[1,`j'] - b1[1,`j']) < 0.25
}
* the seeds that converged before give the estimates of 1.3.0 (seed 1, 4 digits)
assert abs(b1[1,1] - (-1.23734)) < 1e-4 & abs(b1[1,5] - 1.85765) < 1e-4
di as text "PASS: seed 111 converges; seed 1 unchanged"
di "test_nwalaam_patch131: ALL OK"
