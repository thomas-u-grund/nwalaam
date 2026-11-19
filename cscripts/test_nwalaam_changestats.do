cscript
* The change statistics of the plugin against the fixtures of dev/ref (made by dev/gen_ref_definitions.py): for four random
* networks (two undirected, two directed) and every node, the change statistic of every effect is computed from the DEFINITION of
* the statistic (a count over the whole outcome vector, y_i set to 1 minus y_i set to 0) and compared with the plugin.
mata:
real matrix build_adj(real matrix E, real scalar n, real scalar dir)
{
	real matrix A
	real scalar r
	A = J(n, n, 0)
	for (r = 1; r <= rows(E); r++) {
		A[E[r,1], E[r,2]] = 1
		if (!dir) A[E[r,2], E[r,1]] = 1
	}
	return(A)
}
end
foreach c in und_a und_b dir_a dir_b {
	local dir = (substr("`c'", 1, 3) == "dir")
	quietly import delimited "dev/ref/`c'_nodes.csv", clear asdouble
	local n = _N
	mata: Y = st_data(., "y"); CC = st_data(., "c")'; KK = st_data(., "k")'; BB = st_data(., "b")'
	quietly import delimited "dev/ref/`c'_edges.csv", clear asdouble
	mata: A = build_adj(st_data(., ("from", "to")), `n', `dir')
	quietly import delimited "dev/ref/`c'_delta.csv", clear asdouble
	unab names : _all
	mata: REF = st_data(., .)
	// effect codes (nwalaam.ado / alaam_sim.c): the names of the fixture columns in order
	local und "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 43 43 43 44 44 44 45 46 46 46 47 47 47 48"
	local dire "1 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 17 124 125 126 127 128 129 130 131 132 133 134 135 136 137 137 137 138 138 138 139 140 140 140 141 141 141 142"
	mata: st_local("nmcols", strofreal(cols(REF)))
	if `dir' mata: ET = strtoreal(tokens("`dire'"))
	else mata: ET = strtoreal(tokens("`und'"))
	// parameter (alpha / beta) and covariate number of every effect
	mata: EP = J(1, cols(ET), 0); EC = J(1, cols(ET), 0)
	local k = 0
	foreach nm of local names {
		local ++k
		if inlist("`nm'", "gwactivity_0p5", "gwcontagion_0p5", "gwsender_0p5", "gwreceiver_0p5") | regexm("`nm'", "^gw.*_.*5$") mata: EP[`k'] = 0.5
		if regexm("`nm'", "^powercontagion") mata: EP[`k'] = 2
		if regexm("`nm'", "^threshold(contagion|exposure)_([0-9]+)$") {
			local pv = regexs(2)
			mata: EP[`k'] = `pv'
		}
		if regexm("`nm'", "^fraction(contagion|exposure)_([0-9])([0-9])$") {
			local pv = real(regexs(2)) + 0.1 * real(regexs(3))
			mata: EP[`k'] = `pv'
		}
		if regexm("`nm'", "_c$") mata: EC[`k'] = 1
		if regexm("`nm'", "_k$") mata: EC[`k'] = 2
		if regexm("`nm'", "_b$") mata: EC[`k'] = 3
	}
	mata: AlaamChangeStats(A, `dir', Y, ET, EP, EC, (CC \ KK \ BB), D = ., z = .)
	mata: st_numscalar("maxdiff", max(abs(D - REF)))
	mata: st_numscalar("maxrow", select(1::cols(REF), colmax(abs(D - REF))' :== max(abs(D - REF)))[1])
	di "`c': `nmcols' effects, `n' nodes: max abs difference plugin - definition = " maxdiff ""
	if maxdiff > 1e-9 {
		di as error "mismatch in effect " maxrow ": " word("`names'", maxrow)
		exit 9
	}
	assert maxdiff < 1e-9
}
di "PASS: every change statistic equals the definition (4 networks, 2 undirected and 2 directed, all effects, including the threshold, fraction and share contagion effects)"
