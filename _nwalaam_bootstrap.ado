*! version 1.3.1  20nov2026 _nwalaam_bootstrap
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
* nwalaam bootstrap: parametric bootstrap after an nwalaam fit. Draws outcome vectors from the fitted model (the machinery of nwalaam simulate),
* refits each with the same method and effects, and summarises the refitted coefficients. e() is not touched.
capture program drop _nwalaam_bootstrap
program _nwalaam_bootstrap, rclass
	version 14
	syntax [, REPS(integer 0) SEED(string) Level(cilevel) METHod(string) SAVing(string asis) REPLACE NODOTS SILENT THReads(integer 0) ///
		BURNin(string) INTERval(string) SAITer(integer -1) RUNs(integer -1) ITERations(integer -1) THIN(integer -1) FITBURNin(integer -1) ///
		PROPosals(integer -1) LRate(real 0.01)]
	if "`e(cmd)'" != "nwalaam" {
		di as error "{p 0 0 2}{bf:nwalaam bootstrap} needs a preceding {bf:nwalaam} fit.{p_end}"
		error 301
	}
	if `reps' < 2 {
		di as error "reps() must be given and be at least 2"
		error 198
	}
	if `"`seed'"' != "" {
		capture confirm integer number `seed'
		if _rc | real(`"`seed'"') < 0 | real(`"`seed'"') > 2147483647 {
			di as error "seed() must be an integer between 0 and 2147483647"
			error 198
		}
	}
	local method = lower(trim(`"`method'"'))
	if "`method'" == "" local method "`e(method)'"
	if !inlist("`method'", "sa", "ee", "mple") {
		di as error "method() must be sa, ee or mple"
		error 198
	}
	if `threads' < 0 | `threads' > 256 {
		di as error "threads() must be between 1 and 256 (default: one thread per core)"
		error 198
	}
	if `lrate' <= 0 {
		di as error "lrate() must be positive"
		error 198
	}
	foreach o in burnin interval {
		if `"``o''"' != "" {
			capture confirm integer number ``o''
			if _rc | real(`"``o''"') < 1 {
				di as error "`o'() must be a positive integer"
				error 198
			}
		}
	}
	if `"`saving'"' != "" {
		local saving : subinstr local saving `"""' "", all
		if !strpos(`"`saving'"', ".") local saving `"`saving'.dta"'
		if "`replace'" == "" {
			capture confirm new file `"`saving'"'
			if _rc {
				di as error `"file `saving' already exists; use replace"'
				error 602
			}
		}
	}
	// settings of the refits: those of the fit when the method is the fit's, else the defaults of nwalaam
	local same = ("`method'" == "`e(method)'")
	if `saiter' < 0 local saiter = cond(`same' & "`method'" == "sa", e(saiter), 500)
	if `proposals' < 0 local proposals = cond(`same' & "`method'" != "mple", e(proposals), 1000)
	if `runs' < 0 local runs = cond(`same' & "`method'" == "ee", e(runs), 30)
	if `iterations' < 0 local iterations = cond(`same' & "`method'" == "ee", e(iterations), 50000)
	if `thin' < 0 local thin = cond(`same' & "`method'" == "ee", e(thin), 100)
	if `fitburnin' < 0 local fitburnin = cond(`same' & "`method'" == "ee", e(burnin), 1000)
	if `saiter' < 20 | `runs' < 1 | `iterations' < 100 | `thin' < 1 | `fitburnin' < 0 | `proposals' < 1 | `thin' > `iterations' {
		di as error "the settings of the refits are out of range (saiter() >= 20, runs() >= 1, iterations() >= 100, thin() >= 1, proposals() >= 1)"
		error 198
	}
	_nwalaam_libcheck
	tempname TH SE0 S0 B V DG RUNS TRC EST
	matrix `TH' = e(b)
	matrix `SE0' = vecdiag(e(V))
	local P = colsof(`TH')
	local effnames : colnames `TH'
	local depvar `"`e(depvar)'"'
	local netname `"`e(network)'"'
	local effects `"`e(effects_spec)'"'
	local mode = cond("`e(modeopt)'" == "", "1", "`e(modeopt)'")
	local zone `"`e(zone)'"'
	local sampler `"`e(samplerarg)'"'
	// the work is done in a subprogram so that the Mata objects are dropped on every exit
	if `"`seed'"' != "" set seed `seed'
	capture noisily _nwalaam_bs_run, reps(`reps') method(`method') threads(`threads') burnin(`burnin') interval(`interval') saiter(`saiter') ///
		runs(`runs') iterations(`iterations') thin(`thin') fitburnin(`fitburnin') proposals(`proposals') lrate(`lrate') `nodots' `silent' ///
		th(`TH') est(`EST') depvar(`depvar') netname(`netname') effects(`effects') mode(`mode') zone(`zone') sampler(`sampler')
	local __rc = _rc
	local nok = r(nok)
	local nfail = r(nfail)
	local nnoconv = r(nnoconv)
	local simseed = r(simseed)
	local burn = r(burn)
	local intv = r(interval)
	capture mata: mata drop __nwalaam_*
	if `__rc' exit `__rc'
	// the results of the run
	if `nok' < 2 {
		di as error "{p 0 0 2}only `nok' of `reps' refits succeeded; no summary is possible (the fitted model may be degenerate: try other effects or more proposals).{p_end}"
		error 430
	}
	tempname BB SD VB CI BIAS
	mata: __nwalaam_E = st_matrix("`EST'")
	mata: st_matrix("`BB'", mean(__nwalaam_E))
	mata: st_matrix("`VB'", variance(__nwalaam_E))
	mata: st_matrix("`SD'", sqrt(diagonal(variance(__nwalaam_E)))')
	mata: __nwalaam_a = (100 - `level') / 200; __nwalaam_C = J(2, `P', .)
	mata: for (__nwalaam_j = 1; __nwalaam_j <= `P'; __nwalaam_j++) __nwalaam_C[., __nwalaam_j] = (_nwalaam_pct(__nwalaam_E[., __nwalaam_j], __nwalaam_a) \ _nwalaam_pct(__nwalaam_E[., __nwalaam_j], 1 - __nwalaam_a))
	mata: st_matrix("`CI'", __nwalaam_C)
	if `"`saving'"' != "" {
		mata: st_matrix("`EST'", __nwalaam_E)
	}
	capture mata: mata drop __nwalaam_*
	matrix `BIAS' = `BB' - `TH'
	foreach m in BB SD BIAS {
		matrix colnames ``m'' = `effnames'
	}
	matrix colnames `VB' = `effnames'
	matrix rownames `VB' = `effnames'
	matrix colnames `CI' = `effnames'
	matrix rownames `CI' = lower upper
	if `"`saving'"' != "" {
		preserve
		quietly {
			clear
			svmat double `EST', names(b_)
			local j 0
			foreach en of local effnames {
				local ++j
				label variable b_`j' `"`en'"'
			}
			generate long rep = _n
			order rep
			save `"`saving'"', `replace'
		}
		restore
	}
	if "`silent'" == "" {
		di as text _n "{hline 78}"
		di as text "Parametric bootstrap of the " as result "`e(method)'" as text " fit (refits by " as result "`method'" as text "), " as result `reps' as text " outcomes drawn from the fitted model"
		di as text "Successful refits: " as result `nok' as text "   failed: " as result `nfail' as text "   (seed " as result %10.0f `simseed' as text "; simulation burn-in " as result `burn' as text ", interval " as result `intv' as text ")"
		if `nnoconv' > 0 & "`method'" != "mple" di as text "Refits that did not meet the convergence criterion (kept): " as result `nnoconv'
		di as text "{hline 78}"
		di as text %-18s "Coefficient" _col(20) %10s "Estimate" %10s "Model SE" %10s "Boot SD" %10s "Boot bias" %20s "`level'% percentile int."
		di as text "{hline 78}"
		local j 0
		foreach en of local effnames {
			local ++j
			di as text %-18s abbrev("`en'", 18) as result _col(20) %10.4f `TH'[1, `j'] %10.4f sqrt(`SE0'[1, `j']) %10.4f `SD'[1, `j'] %10.4f `BIAS'[1, `j'] ///
				%10.4f `CI'[1, `j'] %10.4f `CI'[2, `j']
		}
		di as text "{hline 78}"
		di as text "{p 0 6 2}Statistics of the successful refits only. Boot bias = mean of the refits minus the fitted estimate. Model-based: the draws come from the fitted model.{p_end}"
	}
	if `nfail' > 0.1 * `reps' di as error "{p 0 6 2}Warning: `nfail' of `reps' refits failed (more than 10 percent); the summary may be biased towards well-behaved draws.{p_end}"
	if `nok' < 50 di as error "{p 0 6 2}Warning: only `nok' refits succeeded (fewer than 50); percentile intervals are unreliable.{p_end}"
	return scalar reps = `reps'
	return scalar nok = `nok'
	return scalar nfail = `nfail'
	return scalar nnoconv = `nnoconv'
	return scalar level = `level'
	return scalar seed = `simseed'
	return scalar burnin = `burn'
	return scalar interval = `intv'
	return matrix b_boot = `BB'
	return matrix sd_boot = `SD'
	return matrix bias_boot = `BIAS'
	return matrix V_boot = `VB'
	return matrix ci_boot = `CI'
	return local method "`method'"
end

// draws the outcome vectors, refits each one; leaves r(nok) r(nfail) r(nnoconv) and the matrix of successful estimates in est()
capture program drop _nwalaam_bs_run
program _nwalaam_bs_run, rclass
	syntax, REPS(integer) METHod(string) THREADS(integer) [BURNin(string) INTERval(string)] SAITER(integer) RUNS(integer) ITERATIONS(integer) ///
		THIN(integer) FITBURNin(integer) PROPOSALS(integer) LRATE(real) TH(name) EST(name) DEPVAR(string) NETNAME(string) EFFECTS(string) ///
		MODE(string) [ZONE(string) SAMPLER(string) NODOTS SILENT]
	_nwalaam_prep `depvar' `netname', effects(`effects') mode(`mode') zone(`zone') sampler(`sampler')
	local P : word count `__names'
	if colsof(`th') != `P' {
		di as error "the effects of the fit (`P') do not match the coefficients; refit the model"
		error 503
	}
	if `"`burnin'"' == "" local burnin = max(10000, 100 * `nodes')
	if `"`interval'"' == "" local interval = max(1000, 10 * `nodes')
	tempname S0 B V DG RUNS TRC
	matrix `S0' = `th'
	local simseed = floor(runiform() * 2147483646) + 1
	mata: __nwalaam_SIM = AlaamSim(__nwalaam_adj, `__dir', __nwalaam_y, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov, ///
		st_matrix("`th'"), `reps', `burnin', `interval', `simseed', 1)
	mata: __nwalaam_E = J(0, `P', .)
	local mple = ("`method'" == "mple")
	local sa = ("`method'" == "sa")
	local nok 0
	local nfail 0
	local nnoconv 0
	if "`nodots'" == "" & "`silent'" == "" di as text "Refitting (. = refit, x = failed):"
	forvalues b = 1/`reps' {
		local __good 0
		// the draw: missing where the observed outcome is missing (the other level of a two-mode network)
		mata: __nwalaam_yb = __nwalaam_SIM[`b', .]'; __nwalaam_yb[selectindex(__nwalaam_y :>= .)] = J(sum(__nwalaam_y :>= .), 1, .)
		mata: st_local("__var", strofreal(sum(__nwalaam_yb :== 1) > 0 & sum(__nwalaam_yb :== 0) > 0))
		if `__var' {
			local sa_converged ""
			local __s = floor(runiform() * 2147483646) + 1
			capture mata: AlaamFit(__nwalaam_adj, `__dir', __nwalaam_yb, __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_cov, ///
				`mple', `runs', `iterations', `thin', `fitburnin', `proposals', `lrate', `__s', `__s' + 17, "`S0'", "`B'", "`V'", "`DG'", ///
				"`RUNS'", "`TRC'", `sa', `saiter', `threads')
			if !_rc {
				mata: st_local("__good", strofreal(!hasmissing(st_matrix("`B'")) & max(abs(st_matrix("`B'"))) < 1e10))
				if `__good' {
					mata: __nwalaam_E = __nwalaam_E \ st_matrix("`B'")
					if "`sa_converged'" == "0" local ++nnoconv
				}
			}
		}
		if `__good' local ++nok
		else local ++nfail
		if "`nodots'" == "" & "`silent'" == "" {
			if `__good' di as text "." _c
			else di as error "x" _c
			if mod(`b', 50) == 0 | `b' == `reps' di as text " `b'"
		}
	}
	mata: st_matrix("`est'", __nwalaam_E)
	return scalar nok = `nok'
	return scalar nfail = `nfail'
	return scalar nnoconv = `nnoconv'
	return scalar simseed = `simseed'
	return scalar burn = `burnin'
	return scalar interval = `interval'
end

version 14
mata:
mata set matastrict off
// percentile of the values x at probability p (linear interpolation of the order statistics, type 7)
real scalar _nwalaam_pct(real colvector x, real scalar p)
{
	real colvector s
	real scalar n, h, lo
	s = sort(x, 1)
	n = rows(s)
	h = (n - 1) * p + 1
	lo = floor(h)
	if (lo >= n) return(s[n])
	return(s[lo] + (h - lo) * (s[lo + 1] - s[lo]))
}
end
