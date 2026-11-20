*! version 1.3.1  20nov2026 nwalaam
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwalaam, a Stata package for autologistic actor attribute models (ALAAMs). It follows the
* method of the ALAAMEE software (Stivala, Wang and Lomi; MIT licence) and is written independently.
* nwalaam is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

capture program drop nwalaam
program nwalaam, eclass
	version 14
	gettoken __w1 __rest : 0, parse(" ,")
	if `"`__w1'"' == "version" {
		nwalaam_version
		exit
	}
	if `"`__w1'"' == "simulate" {
		_nwalaam_simulate `__rest'
		exit
	}
	if replay() {
		capture _nwalaam_tabopts __t check `0'
		if !_rc {
			if "`e(cmd)'" != "nwalaam" error 301
			_nwalaam_display `0'
			exit
		}
	}
	syntax [anything(name=args)] [, EFFects(string) MPLE METHod(string) SAITer(integer 500) RUNs(integer 30) ITERations(integer 50000) THIN(integer 100) ///
		BURNin(integer 1000) PROPosals(integer 1000) LRate(real 0.01) SEED(string) START(name) MODE(integer 1) ZONE(string) SAMPler(string) THReads(integer 0) *]
	local __tabopts `"`options'"'
	local method = lower(trim(`"`method'"'))
	if "`mple'" != "" & "`method'" != "" & "`method'" != "mple" {
		di as error "mple and method(`method') cannot be combined"
		error 198
	}
	if "`method'" == "" local method "sa"
	if "`mple'" != "" local method "mple"
	if "`method'" == "mple" local mple "mple"
	if !inlist("`method'", "ee", "sa", "mple") {
		di as error "method() must be sa (stochastic approximation, the default) or ee (equilibrium expectation)"
		error 198
	}
	if `saiter' < 20 {
		di as error "saiter() must be at least 20"
		error 198
	}
	_nwalaam_tabopts __chk check, `__tabopts'
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
	}
	if `runs' < 1 | `iterations' < 100 | `thin' < 1 | `burnin' < 0 | `proposals' < 1 | `lrate' <= 0 {
		di as error "{p 0 0 2}runs() and proposals() must be at least 1, iterations() at least 100, thin() at least 1, burnin() not negative, lrate() positive.{p_end}"
		error 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256 (default: one thread per core)"
		error 198
	}
	if `thin' > `iterations' {
		di as error "thin() cannot exceed iterations()"
		error 198
	}
	if `"`args'"' == "" {
		di as error "{p 0 0 2}syntax: {bf:nwalaam} {it:outcome} [{it:netname}]{bf:, effects(}{it:list}{bf:)}; see {help nwalaam}.{p_end}"
		error 198
	}
	set more off
	_nwalaam_libcheck
	gettoken depvar netname : args
	local netname = trim("`netname'")
	_nwalaam_prep `depvar' `netname', effects(`effects') mode(`mode') zone(`zone') sampler(`sampler')
	local __P : word count `__names'
	if `__nloops' > 0 di as text "{p 0 6 2}Note: the network `netname' has `__nloops' self-loop(s); they are ignored.{p_end}"
	// start values: given as a 1 x P matrix, or (default) from Algorithm S
	if "`start'" != "" {
		capture confirm matrix `start'
		if _rc {
			di as error "start(): matrix `start' not found"
			error 111
		}
		if rowsof(`start') != 1 | colsof(`start') != `__P' {
			di as error "start() must be a 1 x `__P' matrix, one value per effect"
			error 503
		}
		tempname __S0
		matrix `__S0' = `start'
	}
	else {
		tempname __S0
		matrix `__S0' = J(1, `__P', .)
	}
	// seed: given, or drawn from Stata's own generator (so that set seed makes the fit reproducible)
	if `"`seed'"' == "" local __seed = floor(runiform() * 2147483646) + 1
	else local __seed = `seed'
	tempname B V DG RUNS TRC
	capture noisily mata: AlaamFit(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov, ///
		`=("`mple'" != "")', `runs', `iterations', `thin', `burnin', `proposals', `lrate', `__seed', `__seed' + 17, "`__S0'", "`B'", "`V'", "`DG'", "`RUNS'", "`TRC'", `=("`method'" == "sa")', `saiter', `threads')
	local __rc = _rc
	capture mata: mata drop __nwalaam_*
	if `__rc' exit `__rc'
	matrix colnames `B' = `__names'
	matrix rownames `V' = `__names'
	matrix colnames `V' = `__names'
	matrix rownames `DG' = `__names'
	matrix colnames `DG' = observed pseudolik start ee se tratio mcse
	tempname OBS MP ST TR
	matrix `OBS' = `DG'[1..., 1]'
	matrix `MP' = `DG'[1..., 2]'
	matrix `ST' = `DG'[1..., 3]'
	matrix `TR' = `DG'[1..., 6]'
	matrix colnames `OBS' = `__names'
	matrix colnames `MP' = `__names'
	matrix colnames `ST' = `__names'
	matrix colnames `TR' = `__names'
	if "`mple'" != "" {
		mata: st_local("__hm", strofreal(hasmissing(st_matrix("`B'"))))
		if `__hm' {
			di as error "{p 0 0 2}the pseudo-likelihood estimate does not exist (perfect prediction or collinear effects).{p_end}"
			error 430
		}
	}
	ereturn post `B' `V', obs(`__ny')
	ereturn local cmd "nwalaam"
	ereturn local estat_cmd "nwalaam_estat"
	ereturn local version "1.3.1"
	ereturn local depvar "`depvar'"
	ereturn local network "`netname'"
	ereturn local effects "`__names'"
	ereturn local effects_spec `"`effects'"'
	ereturn local directed = cond(`__dir', "directed", "undirected")
	ereturn local modeopt = cond(`__two', "`mode'", "")
	ereturn local zone "`zone'"
	ereturn local sampler = cond("`zone'" != "", "conditional (snowball)", cond("`sampler'" == "", "toggle", "`sampler'"))
	ereturn local samplerarg "`sampler'"
	ereturn local twomode = cond(`__two', "twomode", "onemode")
	ereturn local method = cond("`mple'" != "", "mple", "`method'")
	ereturn local title = cond("`mple'" != "", "Autologistic actor attribute model (maximum pseudo-likelihood)", cond("`method'" == "sa", "Autologistic actor attribute model (stochastic approximation)", "Autologistic actor attribute model (equilibrium expectation)"))
	ereturn scalar N = `__ny'
	ereturn scalar n_ones = `__ny1'
	ereturn scalar nodes = `nodes'
	ereturn matrix diag = `DG'
	ereturn matrix obs_stats = `OBS'
	ereturn matrix b_mple = `MP'
	if "`mple'" == "" & "`method'" == "sa" {
		ereturn matrix b_start = `ST'
		ereturn matrix tratio = `TR'
		matrix colnames `RUNS' = `__names' `__names' `__names' `__names'
		ereturn matrix run_diag = `RUNS'
		ereturn matrix trace1 = `TRC'
		ereturn scalar runs = 1
		ereturn scalar repetitions = `runs_used'
		ereturn scalar sa_converged = `sa_converged'
		ereturn scalar proposals = `proposals'
		ereturn scalar saiter = `saiter'
		ereturn scalar seed = `__seed'
	}
	else if "`mple'" == "" {
		ereturn matrix b_start = `ST'
		ereturn matrix tratio = `TR'
		matrix colnames `RUNS' = `__names' `__names' `__names' `__names'
		ereturn matrix run_diag = `RUNS'
		ereturn matrix trace1 = `TRC'
		ereturn scalar runs = `runs'
		ereturn scalar runs_used = `runs_used'
		ereturn scalar iterations = `iterations'
		ereturn scalar thin = `thin'
		ereturn scalar burnin = `burnin'
		ereturn scalar proposals = `proposals'
		ereturn scalar acceptance = `accrate'
		ereturn scalar seed = `__seed'
	}
	_nwalaam_display, `__tabopts'
end

capture program drop _nwalaam_display
program _nwalaam_display
	syntax [, *]
	di as text "{hline 70}"
	di as text "`e(title)'"
	di as text "Outcome: " as result "`e(depvar)'" as text " (" as result e(n_ones) as text " of " as result e(N) as text ")" _col(48) as text "Network: " as result "`e(network)'"
	if "`e(method)'" == "sa" di as text "SA repetitions: " as result e(repetitions) as text " (converged: " as result cond(e(sa_converged), "yes", "no") as text ")" _col(48) as text "Proposals per iteration: " as result e(proposals)
	if "`e(method)'" == "ee" di as text "EE runs used: " as result e(runs_used) as text " of " as result e(runs) _col(48) as text "Acceptance rate: " as result %5.3f e(acceptance)
	di as text "{hline 70}"
	_nwalaam_tabopts __t, `options'
	ereturn display, `__t'
	if inlist("`e(method)'", "ee", "sa") {
		local __thr = cond("`e(method)'" == "sa", 0.1, 0.3)
		tempname tr
		matrix `tr' = e(tratio)
		local __k = colsof(`tr')
		local __bad
		local __names : colnames `tr'
		forvalues j = 1/`__k' {
			if abs(`tr'[1, `j']) > `__thr' | missing(`tr'[1, `j']) local __bad `__bad' `: word `j' of `__names''
		}
		if "`__bad'" != "" {
			di as text "{p 0 6 2}Note: the convergence t-ratio exceeds `__thr' in absolute value for " as result "`__bad'" as text ". The estimation may not have converged, or the model may be degenerate; see {bf:matrix list e(tratio)}.{p_end}"
		}
		else di as text "All convergence t-ratios are within `__thr' in absolute value."
	}
end

capture program drop _nwalaam_tabopts
program _nwalaam_tabopts
	version 14
	gettoken __mac 0 : 0, parse(" ,")
	gettoken __chk : 0, parse(" ,")
	if `"`__chk'"' == "check" gettoken __chk 0 : 0, parse(" ,")
	else local __chk
	syntax [, Level(string) *]
	local __lev `"`level'"'
	local ci
	local __rest
	while `"`options'"' != "" {
		gettoken __o options : options, bind
		if strlower(`"`__o'"') == "ci" local ci ci
		else local __rest `"`__rest' `__o'"'
	}
	_get_diopts __di, `__rest'
	local __noci : list posof "noci" in __di
	local __nopv : list posof "nopvalues" in __di
	if "`ci'" != "" & `__noci' {
		di as err "only one of options ci or noci is allowed"
		exit 198
	}
	local __show = ("`ci'" != "" | `"`__lev'"' != "")
	if `"`__lev'"' != "" {
		local 0 , level(`__lev')
		syntax [, Level(cilevel)]
		local __lev `level'
	}
	else local __lev = c(level)
	if "`__chk'" != "" exit
	local __add
	if !`__show' & !`__noci' & !`__nopv' {
		if c(linesize) < 100 | c(lstretch) == "off" | `: list posof "nolstretch" in __di' local __add noci
	}
	c_local `__mac' `"level(`__lev') `__di' `__add'"'
end
