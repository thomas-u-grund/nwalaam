cscript
* nwalaam.ado: option checks, error messages, e() results, replay and display options, reproducibility, missing values,
* directed networks, covariates. Short runs only (the statistical checks are in test_nwalaam_exact.do and test_nwalaam_karate.do).
set seed 123
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
* a small directed network
clear
set obs 120
gen from = ceil(runiform() * 15)
gen to = ceil(runiform() * 15)
drop if from == to
nwfromedge from to, directed name(dnet)
import delimited "dev/karate/karate_nodes.csv", clear
gen double cc = age - 30
gen double ccmiss = cc
replace ccmiss = . in 3
gen y2 = faction
replace y2 = . in 5/6
gen ybad = faction
replace ybad = 2 in 1
gen yconst = 1
gen ystr = "a"

local fast "runs(3) iterations(4000) thin(50) burnin(500)"

* ---- errors
capture nwalaam faction karate, `fast'
assert _rc == 198
capture nwalaam faction karate, effects(nosuch) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(sender) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(density density) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(ego) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(density(3)) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(gwactivity(-1)) `fast'
assert _rc == 198
capture nwalaam faction karate, effects(ego(nosuchvar)) `fast'
assert _rc != 0
capture nwalaam ystr karate, effects(density) `fast'
assert _rc == 198
capture nwalaam ybad karate, effects(density) `fast'
assert _rc == 450
capture nwalaam yconst karate, effects(density) `fast'
assert _rc == 2000
capture nwalaam faction karate, effects(density) runs(0)
assert _rc == 198
capture nwalaam faction karate, effects(density) iterations(10)
assert _rc == 198
capture nwalaam faction karate, effects(density) thin(60000)
assert _rc == 198
capture nwalaam faction karate, effects(density) seed(-1)
assert _rc == 198
capture nwalaam faction karate, effects(density) ci noci `fast'
assert _rc == 198
capture nwalaam faction nonet, effects(density) `fast'
assert _rc != 0
capture nwalaam faction karate, effects(density contagion) start(nomatrix) `fast'
assert _rc == 111
matrix S1 = (0, 0)
capture nwalaam faction karate, effects(density contagion activity) start(S1) `fast'
assert _rc == 503
di as text "  PASS: error paths"

* ---- fit, e() results
nwalaam faction karate, effects(density activity contagion) `fast' seed(5)
assert "`e(cmd)'" == "nwalaam"
assert "`e(depvar)'" == "faction"
assert "`e(network)'" == "karate"
assert "`e(method)'" == "ee"
assert "`e(directed)'" == "undirected"
assert "`e(effects)'" == "density activity contagion"
assert e(N) == 34 & e(n_ones) == 18 & e(nodes) == 34
assert e(runs) == 3 & e(runs_used) >= 1 & e(seed) == 5
assert colsof(e(b)) == 3 & rowsof(e(V)) == 3
matrix OBS = e(obs_stats)
assert OBS[1, 1] == 18 & OBS[1, 2] == 80 & OBS[1, 3] == 35
assert "`: colnames e(b)'" == "density activity contagion"
assert !missing(e(acceptance)) & e(acceptance) > 0 & e(acceptance) < 1
assert "`e(version)'" == "0.1.0"
matrix B1 = e(b)
di as text "  PASS: fit and stored results"

* ---- replay and display options (batch linesize 80: names in full, no intervals)
nwalaam
nwalaam, ci
nwalaam, noci
nwalaam, level(90)
nwalaam, nopvalues
capture nwalaam, ci noci
assert _rc == 198
capture nwalaam, level(5)
assert _rc == 198
estimates store m1
estimates replay m1
matrix B2 = e(b)
assert mreldif(B1, B2) == 0
di as text "  PASS: replay and display options"

* ---- reproducibility
nwalaam faction karate, effects(density contagion) `fast' seed(77)
matrix A1 = e(b)
nwalaam faction karate, effects(density contagion) `fast' seed(77)
matrix A2 = e(b)
assert mreldif(A1, A2) == 0
set seed 4
nwalaam faction karate, effects(density contagion) `fast'
matrix A3 = e(b)
set seed 4
nwalaam faction karate, effects(density contagion) `fast'
matrix A4 = e(b)
assert mreldif(A3, A4) == 0
nwalaam faction karate, effects(density contagion) `fast' seed(78)
matrix A5 = e(b)
assert mreldif(A1, A5) > 0
di as text "  PASS: seed()"

* ---- missing outcome, covariates, long names, parameters in names
nwalaam y2 karate, effects(density contagion ego(cc) gwactivity(0.5) same(class)) `fast' seed(1)
assert e(N) == 32
assert "`e(effects)'" == "density contagion ego_cc gwactivity_0p5 same_class"
nwalaam y2 karate, effects(density contagion ego(ccmiss) alter(cc) diff(class) samecontagion(class)) `fast' seed(1)
nwalaam faction karate, effects(density gwactivity gwcontagion logcontagion powercontagion(3)) `fast' seed(1)
assert "`: colnames e(b)'" == "density gwactivity gwcontagion logcontagion powercontagion_3"
foreach set in "density triangle1 alter2star1" "density twostar partneractivity" "density triangle3 alter2star2" "density threestar triangle2" "density partnerpartner" {
	capture nwalaam faction karate, effects(`set') runs(3) iterations(4000) thin(50) burnin(500) seed(2)
	assert inlist(_rc, 0, 430)
}
di as text "  PASS: missing outcome, covariates, parameters, all undirected effects run"

* ---- mple (karate is separated: r(430); a random outcome is not)
capture nwalaam faction karate, effects(density activity contagion) mple
assert _rc == 430
set seed 99
gen byte yrand = runiform() < 0.4
nwalaam yrand karate, effects(density contagion ego(cc)) mple
assert "`e(method)'" == "mple" & colsof(e(V)) == 3
di as text "  PASS: mple"

* ---- directed network
gen byte yd = runiform() < 0.45 in 1/15
replace yd = . in 16/34
keep in 1/15
nwalaam yd dnet, effects(density reciprocity contagion ego(cc)) runs(2) iterations(3000) thin(50) burnin(500) seed(3)
foreach set in "density sender receiver" "density transitive1 cyclic3" "density gwsender alter(cc)" "density instar2 outstar2 mixedstar" {
	capture nwalaam yd dnet, effects(`set') runs(2) iterations(3000) thin(50) burnin(500) seed(3)
	assert inlist(_rc, 0, 430)
}
assert "`e(directed)'" == "directed"
capture nwalaam yd dnet, effects(density activity) runs(2) iterations(2000) thin(50) burnin(500)
assert _rc == 198
di as text "  PASS: directed network"

* ---- version
nwalaam version
assert "`r(version)'" == "0.1.0"
di "test_nwalaam_ado: ALL OK"
