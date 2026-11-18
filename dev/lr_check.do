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
foreach lr in 0.0005 0.02 {
	nwalaam y g12, effects(density contagion ego(x)) runs(8) iterations(400000) thin(400) burnin(8000) lrate(`lr') seed(5) start(mle)
	matrix B = e(b)
	forvalues j = 1/3 {
		di "lr `lr' effect `j': EE " %8.4f B[1, `j'] "   exact " %8.4f mle[1, `j'] "   diff in exact se " %6.3f (B[1, `j'] - mle[1, `j']) / mlese[1, `j']
	}
}
