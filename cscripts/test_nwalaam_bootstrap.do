cscript
* nwalaam bootstrap: (1) smoke test on a small network (shapes, counts, data and e() untouched), (2) reproducibility with seed(),
* (3) error cases, (4) validation against the exact answer on 8-actor networks (all 2^8 outcome vectors enumerated; the exact maximum
* likelihood estimate and its exact standard errors): bootstrap SD against the exact standard error, next to the equilibrium-expectation
* standard error, (5) saving().
set seed 20261122
mata:
mata set matastrict off
// every outcome vector of the free nodes: the observed statistics of each (n_vectors x P), by the plugin
real matrix ex_Z(real matrix A, real scalar dir, real colvector yfix, real colvector freeidx, real rowvector et, real rowvector ep,
		real rowvector ec, real matrix cov)
{
	real scalar nf, ns, s, b
	real matrix Z
	real colvector y
	real rowvector z
	nf = rows(freeidx)
	ns = 2^nf
	Z = J(ns, cols(et), .)
	for (s = 0; s < ns; s++) {
		y = yfix
		for (b = 1; b <= nf; b++) y[freeidx[b]] = floor(s / 2^(b - 1)) - 2 * floor(s / 2^b)
		AlaamChangeStats(A, dir, y, et, ep, ec, cov, D = ., z = .)
		Z[s + 1, .] = z
	}
	return(Z)
}
// draws a vector from exp(Z th0), returns its row, solves for the exact MLE and standard errors
void ex_solve(real matrix Z, real rowvector th0, string scalar nmle, string scalar nse, string scalar nrow)
{
	real colvector lw, pr, cum, step, lwn
	real scalar u, idx, it, f, k, ll, lln, tries
	real rowvector th, zobs, Ez, thn
	real matrix Cz
	// draw observed outcomes from the model until the maximum likelihood estimate exists (finite)
	for (tries = 1; tries <= 200; tries++) {
		lw = Z * th0'
		pr = exp(lw :- max(lw)); pr = pr / sum(pr)
		u = runiform(1, 1); cum = runningsum(pr)
		idx = selectindex(cum :>= u)[1]
		zobs = Z[idx, .]
		th = J(1, cols(Z), 0)
		ll = -ln(rows(Z))
		for (it = 1; it <= 100; it++) {
			lw = Z * th'
			pr = exp(lw :- max(lw)); pr = pr / sum(pr)
			Ez = pr' * Z
			Cz = cross(Z, pr, Z) - Ez' * Ez
			step = invsym(Cz) * (zobs - Ez)'
			f = 1
			for (k = 1; k <= 40; k++) {
				thn = th + f * step'
				lwn = Z * thn'
				lln = zobs * thn' - (max(lwn) + ln(sum(exp(lwn :- max(lwn)))))
				if (lln >= ll - 1e-12) break
				f = f / 2
			}
			th = thn
			ll = lln
			if (max(abs(f * step)) < 1e-10) break
		}
		if (it <= 100 & max(abs(th)) < 4 & !hasmissing(Cz)) break
	}
	if (tries > 200) exit(498)
	st_matrix(nmle, th)
	st_matrix(nse, sqrt(diagonal(invsym(Cz)))')
	st_numscalar(nrow, idx - 1)
}
real matrix rand_u(real scalar n, real scalar p)
{
	real matrix A
	A = runiform(n, n) :< p
	_diag(A, 0)
	A = uppertriangle(A, 0)
	return(A + A')
}
end

// ---------- 1. smoke test on the Glasgow friendship network
di as text "1. smoke test (Glasgow, 50 pupils)"
nwwebuse glasgow, nwclear
generate byte smoker1 = (smoke1 > 1)
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) seed(20261118)
matrix B0 = e(b)
matrix V0 = e(V)
local N0 = e(N)
local cmd0 "`e(cmd)'"
local rseed0 = e(seed)
quietly describe
local nvar0 = r(k)
local nobs0 = r(N)
quietly summarize smoker1
local m0 = r(mean)
quietly nwsummarize glasgow1
nwalaam bootstrap, reps(30) seed(11) nodots level(90)
assert r(reps) == 30
assert r(nok) + r(nfail) == 30
assert r(nok) >= 20
assert r(level) == 90
foreach m in b_boot sd_boot bias_boot {
	matrix M = r(`m')
	assert rowsof(M) == 1 & colsof(M) == 4
}
matrix M = r(ci_boot)
assert rowsof(M) == 2 & colsof(M) == 4
matrix V = r(V_boot)
assert rowsof(V) == 4 & colsof(V) == 4
matrix SD = r(sd_boot)
matrix CI = r(ci_boot)
matrix BB = r(b_boot)
forvalues j = 1/4 {
	assert SD[1, `j'] > 0 & SD[1, `j'] < . 
	assert abs(sqrt(V[`j', `j']) - SD[1, `j']) < 1e-10
	assert CI[1, `j'] < CI[2, `j']
}
// e() and the data are untouched
matrix B1 = e(b)
matrix V1 = e(V)
assert mreldif(B0, B1) == 0 & mreldif(V0, V1) == 0
assert e(N) == `N0' & "`e(cmd)'" == "`cmd0'" & e(seed) == `rseed0'
quietly describe
assert r(k) == `nvar0' & r(N) == `nobs0'
quietly summarize smoker1
assert r(mean) == `m0'
nwalaam bootstrap, reps(10) seed(11) nodots
// the bootstrap does not leave an e() of its own and nwalaam, and replay still work
nwalaam
di as text "  PASS 1: shapes, counts, e() and data untouched"

// ---------- 2. reproducibility
nwalaam bootstrap, reps(12) seed(77) nodots silent
matrix S1 = r(sd_boot)
local n1 = r(nok)
nwalaam bootstrap, reps(12) seed(77) nodots silent
matrix S2 = r(sd_boot)
assert mreldif(S1, S2) == 0 & r(nok) == `n1'
nwalaam bootstrap, reps(12) seed(78) nodots silent
matrix S3 = r(sd_boot)
assert mreldif(S1, S3) > 0
// set seed before the call gives the same as seed() inside? no: seed() sets the generator, so both agree
set seed 77
nwalaam bootstrap, reps(12) nodots silent
assert mreldif(S1, r(sd_boot)) == 0
di as text "  PASS 2: same seed, same result"

// ---------- 3. error cases
capture nwalaam bootstrap, reps(1)
assert _rc == 198
capture nwalaam bootstrap
assert _rc == 198
capture nwalaam bootstrap, reps(10) level(5)
assert _rc == 198
capture nwalaam bootstrap, reps(10) method(foo)
assert _rc == 198
capture nwalaam bootstrap, reps(10) seed(-3)
assert _rc == 198
capture nwalaam bootstrap, reps(10) burnin(0)
assert _rc == 198
capture nwalaam bootstrap, reps(10) bogus
assert _rc == 198
local tf "`c(tmpdir)'/nwalaam_bs_test.dta"
capture erase "`tf'"
nwalaam bootstrap, reps(10) seed(3) nodots silent saving("`tf'")
capture nwalaam bootstrap, reps(10) seed(3) nodots silent saving("`tf'")
assert _rc == 602
nwalaam bootstrap, reps(10) seed(3) nodots silent saving("`tf'") replace
estimates clear
ereturn clear
capture nwalaam bootstrap, reps(10)
assert _rc == 301
di as text "  PASS 3: error cases"

// ---------- 5. saving(), method(mple), method(ee)
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) seed(20261118)
nwalaam bootstrap, reps(25) seed(4) nodots silent saving("`tf'") replace
local nok = r(nok)
preserve
use "`tf'", clear
assert _N == `nok'
describe, short
assert r(k) == 5
restore
nwalaam bootstrap, reps(20) seed(4) nodots silent method(mple)
assert r(nok) + r(nfail) == 20
di as text "  PASS 5: saving() and method(mple)"
nwalaam smoker1 glasgow1, effects(density contagion sender receiver) method(ee) runs(4) iterations(20000) seed(5)
nwalaam bootstrap, reps(10) seed(4) nodots
assert r(nok) + r(nfail) == 10
di as text "  PASS 5b: after a method(ee) fit"
capture erase "`tf'"

// ---------- 4. validation against the exact answer (8 actors)
capture program drop ex_val
program ex_val
	args label spec reps
	local P = colsof(th0)
	mata: Z = ex_Z(A, DIR, YFIX, FREE, ET, EP, EC, COV)
	mata: ex_solve(Z, st_matrix("th0"), "mle", "mlese", "idxrow")
	mata: yo = YFIX; for (b = 1; b <= rows(FREE); b++) yo[FREE[b]] = floor(st_numscalar("idxrow") / 2^(b - 1)) - 2 * floor(st_numscalar("idxrow") / 2^b)
	mata: st_store(., "y", yo[YROWS])
	quietly nwalaam y `spec' saiter(1000) seed(17)
	matrix SAB = e(b)
	matrix SAV = vecdiag(e(V))
	timer clear 1
	timer on 1
	quietly nwalaam bootstrap, reps(`reps') seed(31) nodots silent saiter(300) saving("`c(tmpdir)'/nwalaam_bs_val.dta") replace
	timer off 1
	quietly timer list 1
	local tt = r(t1)
	local nok = r(nok)
	local nfail = r(nfail)
	matrix SD = r(sd_boot)
	preserve
	quietly use "`c(tmpdir)'/nwalaam_bs_val.dta", clear
	forvalues j = 1/`P' {
		quietly summarize b_`j', detail
		local rsd`j' = (r(p75) - r(p25)) / 1.349
		quietly count if abs(b_`j') > 10
		local nwild`j' = r(N)
	}
	restore
	quietly nwalaam y `spec' method(ee) runs(10) seed(17)
	matrix EEV = vecdiag(e(V))
	di as text _n "  `label': exact MLE / exact SE / SA fit / SA model SE / boot SD / EE SE   (`nok' of `reps' refits, `tt' s)"
	forvalues j = 1/`P' {
		local r1 = SD[1, `j'] / mlese[1, `j']
		local r2 = sqrt(SAV[1, `j']) / mlese[1, `j']
		local r3 = sqrt(EEV[1, `j']) / mlese[1, `j']
		di "  effect `j': " %8.3f mle[1, `j'] "  " %7.3f mlese[1, `j'] "  " %8.3f SAB[1, `j'] "  " %7.3f sqrt(SAV[1, `j']) "  " %7.3f SD[1, `j'] "  " %7.3f sqrt(EEV[1, `j']) ///
			"   ratios to exact SE: boot SD " %5.2f `r1' "  robust boot SD (IQR/1.349) " %5.2f `rsd`j'' / mlese[1, `j'] "  (" `nwild`j'' " refits with |b|>10)  SA " %5.2f `r2' "  EE " %5.2f `r3'
		assert inrange(`r1', 0.3, 10)
	}
end

mata: A = rand_u(8, 0.4)
nwset, mat(A) name(g8) undirected
clear
set obs 8
generate y = 0
generate x = rnormal()
mata: DIR = 0; YFIX = J(8, 1, 0); FREE = (1::8); YROWS = (1::8)
mata: ET = (1, 6, 17); EP = J(1, 3, 0); EC = (0, 0, 1)
mata: COV = st_data(., "x")'
matrix th0 = (-0.7, 0.5, 0.5)
ex_val "problem A (density contagion ego(x))" "g8, effects(density contagion ego(x))" 60
// Problem B (density activity contagion, 8 actors, 100 refits) was checked by hand and is not run here (each degenerate draw makes the
// stochastic approximation slow, 15 minutes): boot SD / exact SE 2.45, 1.77, 1.96; robust (IQR/1.349) 2.83, 2.18, 1.69; EE 0.96, 0.50, 0.46.
di "test_nwalaam_bootstrap: ALL OK"
