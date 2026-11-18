cscript
* nwalaam against the exact model: a network of 16 nodes, so that all 2^16 outcome vectors can be enumerated.
*   1. the observed statistics (plugin) equal the statistics computed by definition on every outcome vector (density, activity,
*      contagion, two-star, triangle1-3: counted directly from the adjacency matrix);
*   2. the exact maximum likelihood estimate and its exact standard errors (inverse Fisher information) are computed;
*   3. the maximum pseudo-likelihood estimate of nwalaam equals Stata's logit without constant on the change statistics;
*   4. the EE estimate (40 runs) agrees with the exact MLE within Monte-Carlo error, and its standard errors with the exact ones.
set seed 20261118
local nn = 16
mata:
mata set matastrict off
real matrix rand_graph(real scalar n, real scalar p)
{
	real matrix A
	A = runiform(n, n) :< p
	_diag(A, 0)
	A = uppertriangle(A, 0)
	A = A + A'
	return(A)
}
// statistics by definition, for one outcome vector: density activity contagion twostar triangle1 triangle2 triangle3
real rowvector z_def(real matrix A, real colvector y)
{
	real scalar n, i, j, k, c, t1, t2, t3, ts
	real colvector d
	n = rows(A)
	d = rowsum(A)
	c = 0
	for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) c = c + A[i,j] * y[i] * y[j]
	ts = 0
	for (i = 1; i <= n; i++) ts = ts + y[i] * d[i] * (d[i] - 1) / 2
	t1 = 0; t2 = 0; t3 = 0
	for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) for (k = j + 1; k <= n; k++) {
		if (A[i,j] & A[i,k] & A[j,k]) {
			t1 = t1 + y[i] + y[j] + y[k]
			t2 = t2 + comb_(y[i] + y[j] + y[k])
			t3 = t3 + (y[i] & y[j] & y[k])
		}
	}
	return((sum(y), sum(y :* d), c, ts, t1, t2, t3))
}
// density, contagion and the ego effect of a covariate, by definition
real rowvector z_def3(real matrix A, real colvector y, real colvector x)
{
	real scalar n, i, j, c
	n = rows(A)
	c = 0
	for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) c = c + A[i,j] * y[i] * y[j]
	return((sum(y), c, sum(y :* x)))
}
real scalar comb_(real scalar m) return(m * (m - 1) / 2)
end
mata: A = rand_graph(`nn', 0.3)
// edge list for nwfromedge
mata: E = selectindex(vec(A')); E = E[selectindex(floor((E :- 1) / `nn') :< mod(E :- 1, `nn'))]; st_local("ne", strofreal(rows(E)))
clear
set obs `ne'
generate from = .
generate to = .
mata: st_store(., "from", floor((E :- 1) / `nn') :+ 1); st_store(., "to", mod(E :- 1, `nn') :+ 1)
nwfromedge from to, undirected name(g12)
clear
set obs `nn'
generate y = .
generate x = rnormal()
_nwsyntax g12
mata: A = *(`netobj'->get_matrix())
mata: st_numscalar("nedge", sum(A) / 2)
assert nedge == `ne'
// 1. statistics by definition vs plugin on random outcome vectors
mata:
et = (1, 2, 6, 3, 11, 12, 13)
ep = J(1, 7, 0); ec = J(1, 7, 0)
nbad = 0
for (r = 1; r <= 60; r++) {
	yy = (runiform(16, 1) :< runiform(1, 1))
	AlaamChangeStats(A, 0, yy, et, ep, ec, J(0, 16, .), D=., z=.)
	if (max(abs(z - z_def(A, yy))) > 1e-9) nbad++
}
st_numscalar("nbad", nbad)
end
assert nbad == 0
di "PASS 1: observed statistics equal the statistics by definition (7 effects, 60 random outcome vectors)"
// 2. exact ML by enumeration
mata:
xx = st_data(., "x")
ns = 2^16
Z = J(ns, 3, .)
for (s = 0; s < ns; s++) {
	yy = J(16, 1, 0)
	for (b = 1; b <= 16; b++) yy[b] = floor(s / 2^(b - 1)) - 2 * floor(s / 2^b)
	Z[s + 1, .] = z_def3(A, yy, xx)
}
th0 = (-0.8, 0.5, 0.6)
lw = Z * th0'
pr = exp(lw :- max(lw)); pr = pr / sum(pr)
u = runiform(1, 1); cum = runningsum(pr)
idx = selectindex(cum :>= u)[1]
zobs = Z[idx, .]
yo = J(16, 1, 0)
s = idx - 1
for (b = 1; b <= 16; b++) yo[b] = floor(s / 2^(b - 1)) - 2 * floor(s / 2^b)
// exact MLE: Newton on  zobs - E_theta[z] = 0
th = J(1, 3, 0)
for (it = 1; it <= 200; it++) {
	lw = Z * th'
	pr = exp(lw :- max(lw)); pr = pr / sum(pr)
	Ez = pr' * Z
	Cz = cross(Z, pr, Z) - Ez' * Ez
	step = invsym(Cz) * (zobs - Ez)'
	th = th + step'
	if (max(abs(step)) < 1e-10) break
}
st_matrix("mle", th)
st_matrix("mlese", sqrt(diagonal(invsym(Cz)))')
st_matrix("yo", yo)
st_numscalar("mleit", it)
end
matrix list mle
matrix list mlese
getmata yy = yo, force
replace y = yy
// 3. pseudo-likelihood = logit
nwalaam y g12, effects(density contagion ego(x)) mple
matrix BM = e(b)
mata: AlaamChangeStats(A, 0, st_data(., "y"), (1, 6, 17), J(1, 3, 0), (0, 0, 1), st_data(., "x")', D=., z=.)
mata: st_matrix("DD", D)
matrix colnames DD = d1 d2 d3
svmat DD, names(col)
logit y d1 d2 d3, noconstant
matrix BL = e(b)
assert mreldif(BM, BL) < 1e-6
di "PASS 3: the pseudo-likelihood estimate equals logit (max rel. diff " mreldif(BM, BL) ")"
// 4. the sampler is stationary at the exact MLE: with the learning rate essentially zero and theta fixed at the exact MLE, the
// accumulated differences between simulated and observed statistics have mean zero (t-ratio |t| < 0.05)
nwalaam y g12, effects(density contagion ego(x)) runs(4) iterations(100000) thin(100) burnin(1000) lrate(1e-12) start(mle) seed(99)
matrix T = e(tratio)
forvalues j = 1/3 {
	assert abs(T[1, `j']) < 0.05
}
di as text "  PASS 4: the sampler is stationary at the exact MLE (t-ratios " %7.4f T[1,1] " " %7.4f T[1,2] " " %7.4f T[1,3] ")"
// 5. EE against the exact MLE: standard errors (information part) within 15 percent of the exact ones; estimates within 0.35 exact
// standard errors (the EE trajectory mean differs from the maximum likelihood estimate by about 0.2 standard errors on this network)
nwalaam y g12, effects(density contagion ego(x)) runs(40) iterations(50000) seed(777)
matrix B = e(b)
matrix SE = vecdiag(e(V))
forvalues j = 1/3 {
	local ee = B[1, `j']
	local se = sqrt(SE[1, `j'])
	local ml = mle[1, `j']
	local mlse = mlese[1, `j']
	di "effect `j': EE " %8.4f `ee' " (se " %6.4f `se' ")   exact " %8.4f `ml' " (se " %6.4f `mlse' ")   diff in exact se " %6.3f (`ee' - `ml') / `mlse' "   se ratio " %5.2f `se' / `mlse'
	assert inrange(`se' / `mlse', 0.85, 1.15)
	assert abs(`ee' - `ml') < 0.35 * `mlse'
}
di as text "  PASS 5: EE estimates and standard errors agree with the exact maximum likelihood"
di "test_nwalaam_exact: ALL OK"
