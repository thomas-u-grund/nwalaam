*! version 0.2.0  18nov2026 _nwalaam_prep
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
* reads the data of a fit into Mata objects (__nwalaam_*): outcome, covariates, the effects (and extra effects appended, for goodness of fit)
capture program drop _nwalaam_prep
program _nwalaam_prep
	version 14
	syntax anything(name=args) , EFFects(string) [EXTRA(string)]
	gettoken depvar netname : args
	local netname = trim("`netname'")
	capture confirm numeric variable `depvar'
	if _rc {
		di as error "{p 0 0 2}the first argument must be the outcome variable (numeric, 0/1) in the current dataset, then the network: {bf:nwalaam outcome netname, effects(...)}.{p_end}"
		error 198
	}
	_nwsyntax `netname', max(1)
	if "`is2mode'" == "true" {
		di as error "nwalaam does not support two-mode (bipartite) networks yet; {bf:`netname'} is two-mode."
		error 198
	}
	local __dir = ("`directed'" == "true")
	if `"`effects'"' == "" {
		di as error "{p 0 0 2}effects() must name at least one effect, e.g. {bf:effects(density contagion activity)}.{p_end}"
		error 198
	}
	capture mata: mata drop __nwalaam_*
	mata: __nwalaam_et = J(1, 0, .); __nwalaam_ep = __nwalaam_et; __nwalaam_ec = __nwalaam_et; __nwalaam_nm = J(1, 0, ""); __nwalaam_cn = J(1, 0, "")
	mata: alaam_effects(`"`effects'"', `__dir', __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_nm, __nwalaam_cn)
	mata: st_local("__names", invtokens(__nwalaam_nm))
	local __Pm : word count `__names'
	if `"`extra'"' != "" {
		mata: __nwalaam_xt = J(1, 0, .); __nwalaam_xp = __nwalaam_xt; __nwalaam_xc = __nwalaam_xt; __nwalaam_xm = J(1, 0, ""); __nwalaam_xcn = J(1, 0, "")
		mata: alaam_effects(`"`extra'"', `__dir', __nwalaam_xt, __nwalaam_xp, __nwalaam_xc, __nwalaam_xm, __nwalaam_xcn)
		mata: (void) alaam_append(__nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_nm, __nwalaam_cn, __nwalaam_xt, __nwalaam_xp, __nwalaam_xc, __nwalaam_xm, __nwalaam_xcn)
		mata: st_local("__allnames", invtokens(__nwalaam_nm))
	}
	else local __allnames `"`__names'"'
	mata: st_local("__covs", invtokens(__nwalaam_cn))
	// per-actor variables: matched to the actors by label when the dataset carries them, else one row per actor in order
	unw_defs
	capture confirm string variable `nw_nodename'
	if !_rc {
		mata: __nwalaam_rows = alaam_node_rows(`netobj'->get_nodenames(), st_sdata(., "`nw_nodename'"))
		mata: st_local("__nf", strofreal(sum(__nwalaam_rows :== 0)))
		mata: st_local("__nd", strofreal(sum(__nwalaam_rows :< 0)))
		if `__nf' > 0 | `__nd' > 0 {
			capture mata: mata drop __nwalaam_rows
			di as error "{p 0 0 2}the actors of `netname' do not each have exactly one row in the current dataset (variable `nw_nodename'); the outcome and the covariates need one row per actor (see {help nwload}'s {bf:xvars} option).{p_end}"
			error 198
		}
	}
	else {
		if _N != `nodes' {
			di as error "{p 0 0 2}the current dataset has `=_N' observations but `netname' has `nodes' actors and no node labels (variable `nw_nodename'); the outcome and the covariates need one row per actor, in actor order.{p_end}"
			error 198
		}
		mata: __nwalaam_rows = (1::`nodes')
	}
	mata: __nwalaam_y = st_data(__nwalaam_rows, "`depvar'")
	mata: st_local("__bad", strofreal(sum((__nwalaam_y :< .) :& (__nwalaam_y :!= 0) :& (__nwalaam_y :!= 1))))
	mata: st_local("__ny1", strofreal(sum(__nwalaam_y :== 1)))
	mata: st_local("__ny", strofreal(sum(__nwalaam_y :< .)))
	if `__bad' > 0 {
		capture mata: mata drop __nwalaam_*
		di as error "{p 0 0 2}the outcome `depvar' must be 0 or 1 (or missing).{p_end}"
		error 450
	}
	if `__ny1' == 0 | `__ny1' == `__ny' {
		capture mata: mata drop __nwalaam_*
		di as error "{p 0 0 2}the outcome `depvar' has no variation among the actors with a non-missing value.{p_end}"
		error 2000
	}
	mata: __nwalaam_cov = J(0, `nodes', .)
	foreach __v of local __covs {
		capture confirm numeric variable `__v'
		if _rc {
			capture mata: mata drop __nwalaam_*
			di as error "{p 0 0 2}covariate {bf:`__v'} is not a numeric variable in the current dataset.{p_end}"
			error 198
		}
		mata: __nwalaam_cov = (__nwalaam_cov \ st_data(__nwalaam_rows, "`__v'")')
	}
	c_local __names `"`__names'"'
	c_local __allnames `"`__allnames'"'
	c_local __Pm `__Pm'
	c_local __covs `"`__covs'"'
	c_local netobj `netobj'
	c_local nodes `nodes'
	c_local __dir `__dir'
	c_local __ny `__ny'
	c_local __ny1 `__ny1'
end

