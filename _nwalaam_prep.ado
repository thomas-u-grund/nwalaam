*! version 1.2.0  19nov2026 _nwalaam_prep
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
* reads the data of a fit into Mata objects (__nwalaam_*): outcome, covariates, the effects (and extra effects appended, for goodness of fit)
capture program drop _nwalaam_prep
program _nwalaam_prep
	version 14
	syntax anything(name=args) , EFFects(string) [EXTRA(string) MODE(integer 1) ZONE(string) SAMPler(string)]
	gettoken depvar netname : args
	local netname = trim("`netname'")
	capture confirm numeric variable `depvar'
	if _rc {
		di as error "{p 0 0 2}the first argument must be the outcome variable (numeric, 0/1) in the current dataset, then the network: {bf:nwalaam outcome netname, effects(...)}.{p_end}"
		error 198
	}
	_nwsyntax `netname', max(1)
	local __two = ("`is2mode'" == "true")
	local __dir = ("`directed'" == "true")
	if `__two' {
		if `__dir' {
			di as error "nwalaam does not support directed two-mode networks."
			error 198
		}
		if !inlist(`mode', 1, 2) {
			di as error "mode() must be 1 or 2: the level of the two-mode network on which the outcome is defined."
			error 198
		}
	}
	else if `mode' != 1 {
		di as error "mode() applies to two-mode networks only."
		error 198
	}
	local sampler = lower(trim("`sampler'"))
	if !inlist("`sampler'", "", "toggle", "zoo") {
		di as error "sampler() must be toggle (the default) or zoo."
		error 198
	}
	if "`sampler'" == "" local sampler "toggle"
	if "`zone'" != "" {
		if `__two' {
			di as error "zone() (snowball sampling) cannot be combined with a two-mode network."
			error 198
		}
		if "`sampler'" == "zoo" {
			di as error "the zero-or-one sampler (sampler(zoo)) cannot be combined with zone(); the snowball (conditional) sampler is used with zone()."
			error 198
		}
		capture confirm numeric variable `zone'
		if _rc {
			di as error "zone() must name a numeric variable: the snowball sampling zone (0, 1, 2, ...) of each actor."
			error 198
		}
	}
	if `"`effects'"' == "" {
		di as error "{p 0 0 2}effects() must name at least one effect, e.g. {bf:effects(density contagion activity)}.{p_end}"
		error 198
	}
	capture mata: mata drop __nwalaam_*
	mata: __nwalaam_et = J(1, 0, .); __nwalaam_ep = __nwalaam_et; __nwalaam_ec = __nwalaam_et; __nwalaam_nm = J(1, 0, ""); __nwalaam_cn = J(1, 0, "")
	mata: __nwalaam_sn = J(1, 0, ""); __nwalaam_dn = J(1, 0, "")
	mata: alaam_effects(`"`effects'"', `__dir', `__two', __nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_nm, __nwalaam_cn, __nwalaam_sn, __nwalaam_dn)
	mata: st_local("__names", invtokens(__nwalaam_nm))
	local __Pm : word count `__names'
	if `"`extra'"' != "" {
		mata: __nwalaam_xt = J(1, 0, .); __nwalaam_xp = __nwalaam_xt; __nwalaam_xc = __nwalaam_xt; __nwalaam_xm = J(1, 0, ""); __nwalaam_xcn = J(1, 0, "")
		mata: __nwalaam_xsn = J(1, 0, ""); __nwalaam_xdn = J(1, 0, "")
		mata: alaam_effects(`"`extra'"', `__dir', `__two', __nwalaam_xt, __nwalaam_xp, __nwalaam_xc, __nwalaam_xm, __nwalaam_xcn, __nwalaam_xsn, __nwalaam_xdn)
		mata: st_local("__xaux", strofreal(cols(__nwalaam_xsn) + cols(__nwalaam_xdn)))
		if `__xaux' > 0 {
			capture mata: mata drop __nwalaam_*
			di as error "extra(): effects that take a setting network or a distance matrix are not allowed here."
			error 198
		}
		mata: (void) alaam_append(__nwalaam_et, __nwalaam_ep, __nwalaam_ec, __nwalaam_nm, __nwalaam_cn, __nwalaam_xt, __nwalaam_xp, __nwalaam_xc, __nwalaam_xm, __nwalaam_xcn)
		mata: st_local("__allnames", invtokens(__nwalaam_nm))
	}
	else local __allnames `"`__names'"'
	mata: st_local("__covs", invtokens(__nwalaam_cn))
	// the network as a symmetric 0/1 tie matrix (for a two-mode network the stored matrix has missing values within a level)
	mata: __nwalaam_adj = *(`netobj'->get_matrix())
	if `__two' {
		mata: __nwalaam_adj = (editmissing(__nwalaam_adj, 0) :!= 0); __nwalaam_adj = (__nwalaam_adj :| __nwalaam_adj'); __nwalaam_adj = __nwalaam_adj :+ 0
		mata: __nwalaam_mo = (strtoreal(`netobj'->get_modes()') :== `mode')
		mata: st_local("__nm", strofreal(sum(__nwalaam_mo)))
	}
	else {
		mata: __nwalaam_mo = J(rows(__nwalaam_adj), 1, 1)
		local __nm `nodes'
	}
	// per-actor variables: matched to the actors by label when the dataset carries them, else one row per actor in order (for a two-mode
	// network: one row per node of both levels, or one row per node of the level of the outcome)
	unw_defs
	capture confirm string variable `nw_nodename'
	if !_rc {
		mata: __nwalaam_rows = alaam_node_rows(`netobj'->get_nodenames(), st_sdata(., "`nw_nodename'"))
		mata: st_local("__nf", strofreal(sum((__nwalaam_rows :== 0) :& __nwalaam_mo)))
		mata: st_local("__nd", strofreal(sum((__nwalaam_rows :< 0) :& __nwalaam_mo)))
		if `__nf' > 0 | `__nd' > 0 {
			capture mata: mata drop __nwalaam_*
			di as error "{p 0 0 2}the actors of `netname' do not each have exactly one row in the current dataset (variable `nw_nodename'); the outcome and the covariates need one row per actor (see {help nwload}'s {bf:xvars} option).{p_end}"
			error 198
		}
		mata: __nwalaam_rows = __nwalaam_rows :* (__nwalaam_rows :> 0)
	}
	else {
		if _N == `nodes' mata: __nwalaam_rows = (1::`nodes')
		else if `__two' & _N == `__nm' {
			mata: __nwalaam_rows = J(`nodes', 1, 0); __nwalaam_rows[selectindex(__nwalaam_mo)] = (1::`__nm')
		}
		else {
			di as error "{p 0 0 2}the current dataset has `=_N' observations but `netname' has `nodes' actors and no node labels (variable `nw_nodename'); the outcome and the covariates need one row per actor, in actor order.{p_end}"
			capture mata: mata drop __nwalaam_*
			error 198
		}
	}
	mata: __nwalaam_y = alaam_pull(__nwalaam_rows, "`depvar'")
	if `__two' mata: __nwalaam_y[selectindex(!__nwalaam_mo)] = J(sum(!__nwalaam_mo), 1, .)
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
		mata: __nwalaam_cov = (__nwalaam_cov \ alaam_pull(__nwalaam_rows, "`__v'")')
	}
	// setting networks and distance matrices of the effects
	mata: st_local("__sets", invtokens(__nwalaam_sn))
	local __k 0
	foreach __s of local __sets {
		local ++__k
		capture _nwsyntax `__s', max(1) other(__st)
		if _rc {
			di as error "setting network {bf:`__s'} not found."
			capture mata: mata drop __nwalaam_*
			error 198
		}
		if "`__stis2mode'" == "true" | `__stnodes' != `nodes' {
			di as error "{p 0 0 2}the setting network {bf:`__s'} must be a one-mode network of the same `nodes' actors as {bf:`netname'}.{p_end}"
			capture mata: mata drop __nwalaam_*
			error 198
		}
		mata: __nwalaam_set`__k' = (editmissing(*(`__stnetobj'->get_matrix()), 0) :!= 0) :+ 0
	}
	mata: st_local("__dists", invtokens(__nwalaam_dn))
	local __k 0
	foreach __s of local __dists {
		local ++__k
		capture mata: __nwalaam_dist`__k' = st_matrix("`__s'")
		mata: st_local("__dok", strofreal(rows(__nwalaam_dist`__k') == `nodes' & cols(__nwalaam_dist`__k') == `nodes'))
		if !`__dok' {
			// not a Stata matrix of the right size: try a Mata matrix of that name
			capture mata: __nwalaam_dist`__k' = `__s'
			capture mata: st_local("__dok", strofreal(rows(__nwalaam_dist`__k') == `nodes' & cols(__nwalaam_dist`__k') == `nodes'))
			if _rc local __dok 0
		}
		if !`__dok' {
			di as error "{p 0 0 2}the distance matrix {bf:`__s'} must be a Stata matrix (or a Mata matrix) with `nodes' rows and columns, one for each actor of {bf:`netname'} in the order of the network.{p_end}"
			capture mata: mata drop __nwalaam_*
			error 198
		}
	}
	// sampler: toggle (default), zero-or-one, or conditional on the snowball sampling structure
	mata: __nwalaam_smode = `=("`sampler'" == "zoo")'
	if "`zone'" != "" {
		mata: __nwalaam_z = alaam_pull(__nwalaam_rows, "`zone'")
		mata: st_local("__zbad", strofreal(hasmissing(__nwalaam_z) | any(__nwalaam_z :!= floor(__nwalaam_z)) | min(__nwalaam_z) != 0))
		mata: st_local("__zmax", strofreal(max(__nwalaam_z)))
		mata: st_local("__zgap", strofreal(rows(uniqrows(__nwalaam_z)) != max(__nwalaam_z) + 1))
		if `__zbad' | `__zgap' {
			capture mata: mata drop __nwalaam_*
			di as error "{p 0 0 2}zone() must hold the snowball zone of every actor: integers 0 (seeds), 1, ..., up to the outermost zone, every zone occurring.{p_end}"
			error 198
		}
		if `__zmax' < 1 {
			capture mata: mata drop __nwalaam_*
			di as error "zone() has only one zone; a snowball sample needs at least the seeds (0) and one more zone."
			error 198
		}
		mata: __nwalaam_fixed = (__nwalaam_z :== `__zmax')
	}
	c_local __names `"`__names'"'
	c_local __allnames `"`__allnames'"'
	c_local __two `__two'
	c_local __Pm `__Pm'
	c_local __covs `"`__covs'"'
	c_local netobj `netobj'
	c_local nodes `nodes'
	c_local __dir `__dir'
	c_local __ny `__ny'
	c_local __ny1 `__ny1'
end

