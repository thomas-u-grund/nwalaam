* nwalaam against the exact model for the extensions: the zero-or-one sampler with setting-network, geographic-homophily and contagion-distance
* effects (12 nodes), a two-mode network with the outcome on one level (four-cycle effects; 5 + 9 nodes), and the snowball (zone) sampler
* (13 nodes, the outermost zone held fixed). For each, all outcome vectors of the nodes that can change are enumerated, the observed outcome is drawn
* from the model, the exact maximum likelihood estimate is computed, and (1) the sampler is stationary at it (equilibrium-expectation t-ratios
* of about 0), (2) the default method (stochastic approximation) agrees with it within 0.3 exact standard errors, with standard errors within 40 percent.
set seed 20261118
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
capture program drop ex_run
// ex_run "<outcome network and options of nwalaam>" : Mata objects A, DIR, YFIX, FREE, ET, EP, EC, COV, YROWS and the Stata matrix th0 are set
program ex_run
	args spec tol
	mata: Z = ex_Z(A, DIR, YFIX, FREE, ET, EP, EC, COV)
	mata: ex_solve(Z, st_matrix("th0"), "mle", "mlese", "idxrow")
	mata: yo = YFIX; for (b = 1; b <= rows(FREE); b++) yo[FREE[b]] = floor(st_numscalar("idxrow") / 2^(b - 1)) - 2 * floor(st_numscalar("idxrow") / 2^b)
	mata: st_store(., "y", yo[YROWS])
	local P = colsof(mle)
	matrix list mle
	matrix list mlese
	// 1. the sampler is stationary at the exact MLE
	quietly nwalaam y `spec' method(ee) runs(4) iterations(100000) thin(100) burnin(1000) lrate(1e-12) start(mle) seed(99)
	matrix T = e(tratio)
	forvalues j = 1/`P' {
		di "  stationarity, effect `j': t-ratio " %7.4f T[1, `j']
		assert abs(T[1, `j']) < 0.08
	}
	di as text "  PASS: the sampler is stationary at the exact MLE"
	// 2. stochastic approximation against the exact MLE
	quietly nwalaam y `spec' saiter(1000) seed(17)
	matrix B = e(b)
	matrix SE = vecdiag(e(V))
	forvalues j = 1/`P' {
		local sa = B[1, `j']
		local se = sqrt(SE[1, `j'])
		local ml = mle[1, `j']
		local mlse = mlese[1, `j']
		di "  effect `j': SA " %8.4f `sa' " (se " %6.4f `se' ")   exact " %8.4f `ml' " (se " %6.4f `mlse' ")   diff in exact se " %6.3f (`sa' - `ml') / `mlse' "   se ratio " %5.2f `se' / `mlse'
		assert inrange(`se' / `mlse', 0.6, 1.4)
		assert abs(`sa' - `ml') < `tol' * `mlse'
	}
	di as text "  PASS: stochastic approximation agrees with the exact maximum likelihood"
end

* ---- 1. zero-or-one sampler; setting network, geographic homophily, contagion distance (12 nodes)
di as text "1. sampler(zoo): setting network, geographic homophily, contagion distance"
mata: A = rand_u(12, 0.3); SM = rand_u(12, 0.25)
mata: DM = runiform(12, 12); DM = (DM + DM') / 2; DM = DM :* (1 :- I(12))
nwset, mat(A) name(g) undirected
nwset, mat(SM) name(sn) undirected
mata: st_matrix("dm", DM)
clear
set obs 12
generate y = 0
generate x = rnormal()
mata: DIR = 0; YFIX = J(12, 1, 0); FREE = (1::12); YROWS = (1::12)
mata: ET = (1, 6, 17, 40, 41, 42); EP = J(1, 6, 0); EC = (0, 0, 1, 1, 1, 1)
mata: COV = st_data(., "x")'
mata: __nwalaam_set1 = SM; __nwalaam_dist1 = DM
matrix th0 = (-0.7, 0.4, 0.5, 0.3, -0.2, 0.3)
ex_run "g, effects(density contagion ego(x) settinghomophily(sn) geographichomophily(dm) contagiondist(dm)) sampler(zoo)" 0.3
capture mata: mata drop __nwalaam_set1 __nwalaam_dist1

* ---- 2. two-mode network, outcome on one level (nine actors on level 1, five events on level 2), four-cycle effects
di as text "2. two-mode network, mode(1): four-cycle effects"
mata: BM = runiform(5, 9) :< 0.4
nwset, bipartite mat(BM) name(tm)
_nwsyntax tm
mata: assert(`netobj'->get_nodes_mode1() == 9 & `netobj'->get_nodes_mode2() == 5)
mata: AM = (editmissing(*(`netobj'->get_matrix()), 0) :!= 0); AM = AM :| AM'; AM = AM :+ 0
clear
set obs 9
generate y = 0
mata: DIR = 0; YFIX = (J(9, 1, 0) \ J(5, 1, .)); FREE = (1::9); YROWS = (1::9)
mata: A = AM
mata: ET = (1, 2, 8, 28, 29); EP = J(1, 5, 0); EC = J(1, 5, 0); COV = J(0, 14, .)
matrix th0 = (-0.5, 0.3, 0.1, 0.1, -0.1)
ex_run "tm, effects(density activity alter2star2 fourcycle1 fourcycle2) mode(1)" 0.3

* ---- 3. snowball (conditional) sampler: the outermost zone is held at its observed values (13 nodes, 5 of them fixed)
di as text "3. zone(): the snowball conditional sampler"
mata: A = rand_u(13, 0.3)
nwset, mat(A) name(sg) undirected
clear
set obs 13
generate y = 0
generate x = rnormal()
generate zone = (_n > 3) + (_n > 8)
mata: yfx = (J(8, 1, 0) \ (1 \ 0 \ 1 \ 1 \ 0))
mata: DIR = 0; YFIX = yfx; FREE = (1::8); YROWS = (1::13)
mata: ET = (1, 6, 17); EP = J(1, 3, 0); EC = (0, 0, 1); COV = st_data(., "x")'
matrix th0 = (-0.6, 0.4, 0.5)
mata: st_store((9::13), "y", yfx[(9::13)])
ex_run "sg, effects(density contagion ego(x)) zone(zone)" 0.3
// the outermost zone keeps its observed outcomes in a simulation
nwalaam simulate y sg, effects(density contagion ego(x)) theta(th0) zone(zone) nsim(50) frame(simfr) replace
frame simfr: mata: st_numscalar("fixedok", all(st_data(9::13, "sim1") :== (1 \ 0 \ 1 \ 1 \ 0)))
assert fixedok == 1
di as text "  PASS: simulated outcomes keep the outermost zone fixed"
di "test_nwalaam_exact2: ALL OK"
