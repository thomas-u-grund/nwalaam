*! version 0.3.0  18nov2026 _nwalaam_simulate
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
* nwalaam simulate: draws outcome vectors from the model (after a fit, or from given coefficients)
capture program drop _nwalaam_simulate
program _nwalaam_simulate, rclass
	version 14
	syntax [anything(name=args)] [, EFFects(string) THETA(name) NSIM(integer 1000) BURNin(string) INTERval(string) SEED(string) ///
		GENerate(name) FRame(name) REPlace MODE(integer 1) ZONE(string) SAMPler(string)]
	_nwalaam_libcheck
	local usefit = 0
	if `"`args'"' == "" {
		if "`e(cmd)'" != "nwalaam" {
			di as error "{p 0 0 2}{bf:nwalaam simulate} without an outcome and a network needs a preceding {bf:nwalaam} fit; otherwise: {bf:nwalaam simulate} {it:outcome netname}{bf:, effects(}{it:list}{bf:) theta(}{it:matrix}{bf:)}.{p_end}"
			error 301
		}
		local usefit = 1
		gettoken depvar netname : args
		local depvar `"`e(depvar)'"'
		local netname `"`e(network)'"'
		local effects `"`e(effects_spec)'"'
		local mode = cond("`e(modeopt)'" == "", 1, `e(modeopt)')
		local zone `"`e(zone)'"'
		local sampler `"`e(samplerarg)'"'
	}
	else {
		gettoken depvar netname : args
		local netname = trim("`netname'")
		if `"`effects'"' == "" | "`theta'" == "" {
			di as error "{p 0 0 2}give effects() and theta() (a 1 x P coefficient matrix, one value per effect) when no model has been fitted.{p_end}"
			error 198
		}
	}
	if `nsim' < 1 {
		di as error "nsim() must be at least 1"
		error 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
		local __seed = `seed'
	}
	else local __seed = floor(runiform() * 2147483646) + 1
	if "`generate'" == "" & "`frame'" == "" {
		di as error "{p 0 0 2}give generate(stub) to store the draws as variables stub1, stub2, ... or frame(name) to store them in a frame.{p_end}"
		error 198
	}
	tempname TH
	if "`theta'" != "" {
		capture confirm matrix `theta'
		if _rc {
			di as error "theta(): matrix `theta' not found"
			error 111
		}
		matrix `TH' = `theta'
	}
	else matrix `TH' = e(b)
	_nwalaam_prep `depvar' `netname', effects(`effects') mode(`mode') zone(`zone') sampler(`sampler')
	local __P : word count `__names'
	if colsof(`TH') != `__P' {
		di as error "the coefficient matrix has `=colsof(`TH')' values but the model has `__P' effects"
		error 503
	}
	if `"`burnin'"' == "" local burnin = max(10000, 100 * `nodes')
	if `"`interval'"' == "" local interval = max(1000, 10 * `nodes')
	if "`generate'" != "" & "`replace'" == "" {
		forvalues j = 1/`=min(`nsim', 3)' {
			capture confirm new variable `generate'`j'
			if _rc {
				di as error "variable `generate'`j' already exists; use replace or another stub"
				error 110
			}
		}
	}
	tempname SIM
	mata: st_matrix("`SIM'", AlaamSim(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov, ///
		st_matrix("`TH'"), `nsim', `burnin', `interval', `__seed', 1))
	local __rc = _rc
	if `__rc' {
		capture mata: mata drop __nwalaam_*
		exit `__rc'
	}
	// SIM: nsim x n; write the draws as variables (rows of the dataset = actors, as the outcome) or into a frame
	if "`generate'" != "" {
		forvalues j = 1/`nsim' {
			if "`replace'" != "" capture drop `generate'`j'
			quietly gen byte `generate'`j' = .
		}
		mata: __nwalaam_S = st_matrix("`SIM'")
		forvalues j = 1/`nsim' {
			mata: st_store(__nwalaam_rows, "`generate'`j'", __nwalaam_S[`j', .]')
		}
		mata: mata drop __nwalaam_S
	}
	if "`frame'" != "" {
		capture frame drop `frame'
		frame create `frame'
		mata: __nwalaam_S = st_matrix("`SIM'")
		frame `frame' {
			quietly set obs `nodes'
			quietly gen long actor = _n
			forvalues j = 1/`nsim' {
				quietly gen byte sim`j' = .
				mata: st_store(., "sim`j'", __nwalaam_S[`j', .]')
			}
		}
		mata: mata drop __nwalaam_S
	}
	// statistics of the draws for the model's effects (mean), compared with the observed ones
	tempname OBS MN
	mata: st_matrix("`OBS'", AlaamGofStats(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov))
	matrix colnames `OBS' = `__names'
	capture mata: mata drop __nwalaam_*
	di as text "Simulated " as result `nsim' as text " outcome vectors on " as result `nodes' as text " actors" _c
	if "`generate'" != "" di as text ", stored in variables " as result "`generate'1" as text " to " as result "`generate'`nsim'" _c
	if "`frame'" != "" di as text ", stored in frame " as result "`frame'" _c
	di as text " (seed " as result `__seed' as text ", burn-in " as result `burnin' as text " and interval " as result `interval' as text " proposals)."
	return scalar nsim = `nsim'
	return scalar seed = `__seed'
	return scalar burnin = `burnin'
	return scalar interval = `interval'
	return matrix observed = `OBS'
end

