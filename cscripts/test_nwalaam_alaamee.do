* The change statistics of the plugin against the Python ALAAMEE software: dev/ref_alaamee holds, for six random networks (undirected,
* directed, two-mode with the outcome on either level; outcomes, covariates and ties partly missing) and every node, the change statistic of
* every effect as ALAAMEE computes it (dev/gen_ref_alaamee.py; only the numbers are stored). Two-mode: only the nodes of the outcome level.
mata:
real matrix alaam_adj(real matrix E, real scalar n, real scalar sym)
{
	real matrix A
	real scalar r
	A = J(n, n, 0)
	for (r = 1; r <= rows(E); r++) {
		A[E[r,1], E[r,2]] = 1
		if (sym) A[E[r,2], E[r,1]] = 1
	}
	return(A)
}
end
local tot 0
foreach c in x_und_a x_und_b x_dir_a x_dir_b x_bip_b x_bip_a {
	quietly import delimited "dev/ref_alaamee/`c'_meta.csv", clear asdouble
	local n = n[1]
	local dir = directed[1]
	local bip = (substr("`c'", 3, 3) == "bip")
	quietly import delimited "dev/ref_alaamee/`c'_nodes.csv", clear asdouble
	mata: Y = st_data(., "y"); CC = st_data(., "c")'; KK = st_data(., "k")'; BB = st_data(., "b")'
	quietly import delimited "dev/ref_alaamee/`c'_edges.csv", clear asdouble
	mata: A = alaam_adj(st_data(., ("from", "to")), `n', !`dir')
	quietly import delimited "dev/ref_alaamee/`c'_set.csv", clear asdouble
	capture drop __nwalaam_set1
	mata: SETM = alaam_adj(st_data(., ("from", "to")), `n', 1)
	quietly import delimited "dev/ref_alaamee/`c'_dist.csv", clear asdouble
	mata: DISTM = st_data(., .)
	quietly import delimited "dev/ref_alaamee/`c'_spec.csv", clear asdouble stringcols(1)
	mata: ET = st_data(., "code")'; EP = st_data(., "param")'; KD = st_data(., "cov")'
	mata: EC = (KD :== 1) :+ 2 * (KD :== 2) :+ 3 * (KD :== 3) :+ (KD :== 4) :+ (KD :== 5)
	mata: NAMES = st_sdata(., "name")
	quietly import delimited "dev/ref_alaamee/`c'_delta.csv", clear asdouble
	mata: REF = st_data(., .)
	capture mata: mata drop __nwalaam_set1 __nwalaam_dist1
	mata: __nwalaam_set1 = SETM; __nwalaam_dist1 = DISTM
	mata: AlaamChangeStats(A, `dir', Y, ET, EP, EC, (CC \ KK \ BB), D = ., z = .)
	mata: ok = selectindex(rowmissing(REF) :== 0)
	mata: st_numscalar("maxdiff", max(abs(D[ok, .] - REF[ok, .])))
	mata: st_numscalar("maxcol", select(1::cols(REF), colmax(abs(D[ok, .] - REF[ok, .]))' :== max(abs(D[ok, .] - REF[ok, .])))[1])
	mata: st_local("worst", NAMES[st_numscalar("maxcol")])
	mata: st_local("ne", strofreal(cols(REF)))
	di "`c': `ne' effects, `n' nodes: max abs difference plugin - ALAAMEE = " maxdiff
	if maxdiff > 1e-9 {
		di as error "mismatch in effect `worst'"
		exit 9
	}
	local ++tot
}
di "PASS: every change statistic equals ALAAMEE's (`tot' networks: undirected, directed, two-mode)"
capture mata: mata drop __nwalaam_set1 __nwalaam_dist1
