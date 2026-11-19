*! version 1.2.0  19nov2026 nwalaam_estat
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).

capture program drop nwalaam_estat
program nwalaam_estat, rclass
	version 14
	_nwalaam_libcheck
	gettoken subcmd rest : 0, parse(" ,")
	local subcmd = lower(trim("`subcmd'"))
	if "`e(cmd)'" != "nwalaam" error 301
	if "`subcmd'" == "gof" {
		_nwalaam_gof `rest'
		return add
	}
	else if "`subcmd'" == "mcmcdiag" {
		_nwalaam_mcmcdiag `rest'
		return add
	}
	else {
		di as err "estat `subcmd' not allowed after nwalaam"
		exit 321
	}
end

* estat gof: simulates outcome vectors from the fitted model and compares the statistics with the observed ones
capture program drop _nwalaam_gof
program _nwalaam_gof, rclass
	version 14
	syntax [, EXTRA(string) NOEXtra NSIM(integer 1000) BURNin(string) INTERval(string) SEED(string) PLOT NAME(string) ALL]
	if `nsim' < 20 {
		di as error "nsim() must be at least 20"
		exit 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			exit 198
		}
		local __seed = `seed'
	}
	else local __seed = floor(runiform() * 2147483646) + 1
	local depvar `"`e(depvar)'"'
	local netname `"`e(network)'"'
	local effects `"`e(effects_spec)'"'
	if "`noextra'" == "" & `"`extra'"' == "" {
		if "`e(directed)'" == "directed" local extra "sender receiver reciprocity outstar2 instar2 mixedstar mixedsource mixedsink contagion contagionrecip transitive1 transitive3 transitived1 transitiveu1 cyclic1 cyclic3 alterinstar2 alteroutstar2"
		else if "`e(twomode)'" == "twomode" local extra "twostar threestar alter2star1 alter2star2 fourcycle1 fourcycle2"
		else local extra "twostar threestar alter2star1 triangle1 contagion alter2star2 partneractivity partnerpartner triangle2 triangle3"
	}
	if "`noextra'" != "" local extra
	tempname TH
	matrix `TH' = e(b)
	local __P = colsof(`TH')
	local __mo = cond("`e(modeopt)'" == "", "1", "`e(modeopt)'")
	_nwalaam_prep `depvar' `netname', effects(`effects') extra(`extra') mode(`__mo') zone(`e(zone)') sampler(`e(samplerarg)')
	local __Pall : word count `__allnames'
	mata: __nwalaam_n = rows(__nwalaam_adj)
	if `"`burnin'"' == "" local burnin = max(10000, 100 * `nodes')
	if `"`interval'"' == "" local interval = max(1000, 10 * `nodes')
	tempname OBS SIM
	mata: st_matrix("`OBS'", AlaamGofStats(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov))
	mata: st_matrix("`SIM'", AlaamSim(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov, ///
		(st_matrix("`TH'"), J(1, `__Pall' - `__P', 0)), `nsim', `burnin', `interval', `__seed', 0))
	capture mata: mata drop __nwalaam_*
	tempname RES MAH
	mata: AlaamGofTable("`OBS'", "`SIM'", "`RES'", "`MAH'")
	matrix rownames `RES' = `__allnames'
	matrix colnames `RES' = observed mean sd tratio
	local w = 12
	foreach nm of local __allnames {
		local w = max(`w', min(28, udstrlen("`nm'")))
	}
	local w1 = `w' + 1
	di
	di as txt "Goodness of fit: " as res trim(string(`nsim',"%12.0fc")) as txt " simulated outcome vectors from the fitted model (burn-in " as res trim(string(`burnin',"%12.0fc")) as txt ", interval " as res trim(string(`interval',"%12.0fc")) as txt ")"
	di as txt "{hline `w1'}{c TT}{hline 44}"
	di as txt %`w's "Statistic" " {c |}" %11s "Observed" %11s "Sim. mean" %11s "Sim. SD" %11s "t-ratio"
	di as txt "{hline `w1'}{c +}{hline 44}"
	local k = 0
	local nbad = 0
	foreach nm of local __allnames {
		local ++k
		if `k' == `__P' + 1 {
			di as txt "{hline `w1'}{c +}{hline 44}"
			di as txt %`w's "not in the model:" " {c |}"
		}
		local t = `RES'[`k', 4]
		local flag
		if missing(`t') local flag " ."
		else if abs(`t') > 0.3 & `k' <= `__P' local flag " *"
		else if abs(`t') > 0.3 local flag " *"
		if abs(`t') > 0.3 & !missing(`t') local ++nbad
		di as txt %`w's abbrev("`nm'", `w') " {c |}" as res %11.3f `RES'[`k', 1] %11.3f `RES'[`k', 2] %11.3f `RES'[`k', 3] %11.3f `t' as txt "`flag'"
	}
	di as txt "{hline `w1'}{c BT}{hline 44}"
	di as txt "Mahalanobis distance of the observed statistics from the simulated ones: " as res %9.3f `MAH'[1, 1]
	di as txt "{p 0 0 2}t-ratio = (simulated mean - observed) / simulated SD; * marks |t| > 0.3 (t-ratios of the statistics in the model should be near 0 after a converged fit; for the other statistics, values within about 0.3 indicate a good fit).{p_end}"
	if "`plot'" != "" {
		_nwalaam_gofplot, res(`RES') p(`__P') name(`name')
	}
	return scalar nsim = `nsim'
	return scalar seed = `__seed'
	return scalar mahalanobis = `MAH'[1, 1]
	return scalar n_extreme = `nbad'
	return matrix gof = `RES'
end

capture program drop _nwalaam_gofplot
program _nwalaam_gofplot
	syntax, RES(name) P(integer) [NAME(string)]
	preserve
	tempname R
	matrix `R' = `res'
	local nm : rownames `R'
	local K = rowsof(`R')
	clear
	quietly set obs `K'
	quietly gen str32 stat = ""
	quietly gen double t = .
	quietly gen byte inmodel = .
	local k = 0
	foreach s of local nm {
		local ++k
		quietly replace stat = "`s'" in `k'
		quietly replace t = `R'[`k', 4] in `k'
		quietly replace inmodel = (`k' <= `p') in `k'
	}
	quietly gen int pos = _n
	quietly levelsof pos, local(ps)
	local lab
	foreach i of local ps {
		local lab `"`lab' `i' "`=stat[`i']'""'
	}
	local nopt
	if "`name'" != "" local nopt name(`name', replace)
	twoway (bar t pos if inmodel == 1, horizontal barwidth(0.7)) (bar t pos if inmodel == 0, horizontal barwidth(0.7)), ///
		xline(-0.3 0.3, lpattern(dash)) ylabel(`lab', angle(0) valuelabel nogrid labsize(small)) yscale(reverse) ///
		xtitle("Goodness-of-fit t-ratio") ytitle("") legend(order(1 "in the model" 2 "not in the model") rows(1)) ///
		title("Goodness of fit") `nopt'
	restore
end

* estat mcmcdiag: convergence diagnostics of the equilibrium-expectation runs
capture program drop _nwalaam_mcmcdiag
program _nwalaam_mcmcdiag, rclass
	version 14
	syntax [, PLOT NAME(string) RUN(integer 1) THRESHold(real 0.3)]
	if !inlist("`e(method)'", "ee", "sa") {
		di as err "{p 0 0 2}{bf:estat mcmcdiag} is only available after an equilibrium-expectation or stochastic-approximation fit; a {bf:mple} fit runs no simulation.{p_end}"
		exit 498
	}
	tempname RD TR1 b
	matrix `RD' = e(run_diag)
	matrix `TR1' = e(trace1)
	matrix `b' = e(b)
	local P = colsof(`b')
	local names : colnames `b'
	local nr = rowsof(`RD')
	tempname SUM
	matrix `SUM' = J(`P', 7, .)
	mata: AlaamDiagTable("`RD'", `P', `threshold', "`SUM'")
	matrix rownames `SUM' = `names'
	matrix colnames `SUM' = est_mean est_sd t_mean t_maxabs nrun_bad acf_theta acf_dza
	local w = 12
	foreach nm of local names {
		local w = max(`w', min(28, udstrlen("`nm'")))
	}
	local w1 = `w' + 1
	di
	if "`e(method)'" == "sa" {
		di as txt "Stochastic-approximation diagnostics: " as res e(repetitions) as txt " repetition(s) of the three phases, " as res trim(string(e(proposals),"%12.0fc")) as txt " proposals per iteration, " as res trim(string(4*e(saiter),"%12.0fc")) as txt " iterations in the final phase"
		di as txt "{hline `w1'}{c TT}{hline 40}"
		di as txt %`w's "Parameter" " {c |}" %10s "Estimate" %10s "t" %10s "|t| max" %10s "AC dz"
		di as txt "{hline `w1'}{c +}{hline 40}"
		forvalues k = 1/`P' {
			local nm : word `k' of `names'
			di as txt %`w's abbrev("`nm'", `w') " {c |}" as res %10.4f `SUM'[`k', 1] %10.4f `SUM'[`k', 3] %10.4f `SUM'[`k', 4] %10.3f `SUM'[`k', 7]
		}
		di as txt "{hline `w1'}{c BT}{hline 40}"
		di as txt "{p 0 0 2}t is the convergence t-ratio, mean(dz)/sd(dz) over the iterations of the final phase after a burn-in of one tenth, where dz is the difference between the simulated and the observed statistics; |t| below 0.1 counts as converged (the estimation is repeated, up to 20 times, until it is reached). AC is the lag-1 autocorrelation of dz.{p_end}"
	}
	else {
	di as txt "Equilibrium-expectation diagnostics: " as res e(runs_used) as txt " runs used of " as res e(runs) as txt ", " as res trim(string(e(iterations),"%12.0fc")) as txt " iterations each, thinned to every " as res e(thin) as txt ", burn-in " as res trim(string(e(burnin),"%12.0fc"))
	di as txt "Acceptance rate " as res %5.3f e(acceptance)
	di as txt "{hline `w1'}{c TT}{hline 66}"
	di as txt %`w's "Parameter" " {c |}" %10s "Mean" %10s "SD runs" %10s "Mean t" %10s "Max |t|" %9s "Runs>`threshold'" %9s "AC theta" %9s "AC dzA"
	di as txt "{hline `w1'}{c +}{hline 66}"
	forvalues k = 1/`P' {
		local nm : word `k' of `names'
		di as txt %`w's abbrev("`nm'", `w') " {c |}" as res %10.4f `SUM'[`k', 1] %10.4f `SUM'[`k', 2] %10.4f `SUM'[`k', 3] %10.4f `SUM'[`k', 4] %9.0f `SUM'[`k', 5] %9.3f `SUM'[`k', 6] %9.3f `SUM'[`k', 7]
	}
	di as txt "{hline `w1'}{c BT}{hline 66}"
	di as txt "{p 0 0 2}t is the convergence t-ratio of a run, mean(dzA)/sd(dzA) over the saved iterations after the burn-in; |t| <= `threshold' counts as converged. AC is the lag-1 autocorrelation of the saved (thinned) values of the coefficient and of dzA; values near 1 mean the thinning interval is too short for the number of runs to give independent information.{p_end}"
	}
	if "`plot'" != "" {
		if `run' != 1 {
			di as txt "(only the trace of the first run is stored; run() ignored)"
		}
		if "`e(method)'" == "sa" _nwalaam_traceplot, trace(`TR1') p(`P') names(`names') thin(1) burnin(`=ceil(0.1*4*e(saiter))') name(`name') dzonly
		else _nwalaam_traceplot, trace(`TR1') p(`P') names(`names') thin(`=e(thin)') burnin(`=e(burnin)') name(`name')
	}
	return scalar runs_used = e(runs_used)
	return scalar threshold = `threshold'
	return matrix diag = `SUM'
end

capture program drop _nwalaam_traceplot
program _nwalaam_traceplot
	syntax, TRACE(name) P(integer) NAMES(string) THIN(integer) BURNin(integer) [NAME(string) DZONLY]
	preserve
	tempname T
	matrix `T' = `trace'
	local ns = rowsof(`T')
	clear
	quietly set obs `ns'
	quietly gen double it = (_n - 1) * `thin'
	local k = 0
	local graphs
	foreach nm of local names {
		local ++k
		quietly gen double th`k' = .
		quietly gen double dz`k' = .
		forvalues i = 1/`ns' {
			quietly replace th`k' = `T'[`i', `k'] in `i'
			quietly replace dz`k' = `T'[`i', `P' + `k'] in `i'
		}
		if "`dzonly'" == "" quietly twoway line th`k' it, title("`nm': coefficient", size(medsmall)) xtitle("iteration") ytitle("") xline(`burnin', lpattern(dash)) name(__nwalaam_g`k'a, replace) nodraw
		quietly twoway line dz`k' it, title("`nm': dzA", size(medsmall)) xtitle("iteration") ytitle("") xline(`burnin', lpattern(dash)) yline(0) name(__nwalaam_g`k'b, replace) nodraw
		if "`dzonly'" == "" local graphs `graphs' __nwalaam_g`k'a __nwalaam_g`k'b
		else local graphs `graphs' __nwalaam_g`k'b
	}
	local nopt
	if "`name'" != "" local nopt name(`name', replace)
	graph combine `graphs', cols(`=cond("`dzonly'" == "", 2, min(2, `k'))') title("`=cond("`dzonly'" == "", "Equilibrium-expectation trace of the first run", "Stochastic approximation: statistics of the final phase")'") `nopt'
	forvalues j = 1/`k' {
		capture graph drop __nwalaam_g`j'a
		capture graph drop __nwalaam_g`j'b
	}
	restore
end
