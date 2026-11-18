/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwalaam, a Stata package for autologistic actor attribute models (ALAAMs). It follows the
 * method of the ALAAMEE software (Stivala, Wang and Lomi; MIT licence) and is written independently.
 * nwalaam is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	unw_alaam.do -- Mata glue of nwalaam (built into lib/lnwalaam.mlib). Every computation on the network runs in the compiled
	plugin native/alaam_sim.c (change statistics, observed statistics, the equilibrium-expectation sampler); Mata sorts out
	the effects, writes the model into the plugin's frame, fits the pseudo-likelihood start (a logit), runs Algorithm S and
	the EE runs and computes the standard errors from the runs (batch means). There is no Mata fallback: without a plugin
	the command stops (r(601), r(3495)).
*/
capture mata: mata drop AlaamPluginFilename()
capture mata: mata drop AlaamPluginSubdir()
capture mata: mata drop AlaamPluginPath()
capture mata: mata drop AlaamNativeAvailable()
capture mata: mata drop ALAAM_PROTOCOL_NEEDED()
capture mata: mata drop alaam_node_rows()
capture mata: mata drop alaam_pull()
capture mata: mata drop AlaamNat()
capture mata: mata drop alaam_effects()
capture mata: mata drop alaam_logit()
capture mata: mata drop alaam_fixed()
capture mata: mata drop alaam_bmcov()
capture mata: mata drop AlaamFit()
capture mata: mata drop AlaamDeltas()
capture mata: mata drop AlaamStats()
capture mata: mata drop AlaamChangeStats()
capture mata: mata drop alaam_acf1()
capture mata: mata drop alaam_append()
capture mata: mata drop AlaamSim()
capture mata: mata drop AlaamSA()
capture mata: mata drop AlaamGofTable()
capture mata: mata drop AlaamDiagTable()
capture mata: mata drop AlaamGofStats()

mata:
mata set matastrict off

real colvector alaam_node_rows(string rowvector names, string colvector lab)
{
	transmorphic A
	real colvector pos
	real scalar k

	A = asarray_create("string", 1)
	for (k = 1; k <= rows(lab); k++) {
		if (asarray_contains(A, lab[k])) asarray(A, lab[k], -1)
		else asarray(A, lab[k], k)
	}
	pos = J(cols(names), 1, 0)
	for (k = 1; k <= cols(names); k++) {
		if (asarray_contains(A, names[k])) pos[k] = asarray(A, names[k])
	}
	return(pos)
}

// one value per node from the dataset rows of the nodes (0 = the node has no row: missing)
real colvector alaam_pull(real colvector rows, string scalar var)
{
	real colvector v, out
	real scalar k

	out = J(rows(rows), 1, .)
	v = st_data(., var)
	for (k = 1; k <= rows(rows); k++) if (rows[k] > 0) out[k] = v[rows[k]]
	return(out)
}

string scalar AlaamPluginFilename()
{
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("alaam_sim_windows.plugin")
	if (os == "Unix") return("alaam_sim_unix.plugin")
	return("alaam_sim_macos.plugin")
}

string scalar AlaamPluginSubdir()
{
	string scalar os

	os = st_global("c(os)")
	if (os == "Windows") return("windows")
	if (os == "Unix") return("unix")
	return("macos")
}

// (1) findfile() on the file name alone (after a net install), (2) lib/plugins/<os>/ next to nwalaam.ado (a checkout)
string scalar AlaamPluginPath()
{
	string scalar fname, found, full, dir, fn

	fname = AlaamPluginFilename()
	found = findfile(fname)
	if (found != "") return(found)
	full = findfile("nwalaam.ado")
	if (full == "") return("")
	pathsplit(full, dir, fn)
	return(pathjoin(pathjoin(dir, "lib"), pathjoin("plugins", pathjoin(AlaamPluginSubdir(), fname))))
}

real scalar AlaamNativeAvailable()
{
	string scalar p

	p = AlaamPluginPath()
	if (p == "") return(0)
	return(fileexists(p))
}

real scalar ALAAM_PROTOCOL_NEEDED() return(5)

/*
	The plugin's frame. setup() writes the model (see the header of native/alaam_sim.c); call() runs one plugin code and
	returns the first nres rows of the result variable.
*/
class AlaamNat {
	real scalar isopen, L, nres
	real colvector blob
	string scalar origframe
	void open()
	void close()
	void setup()
	real colvector call()
}

void AlaamNat::open()
{
	string scalar p, fname, q
	real scalar rc, ver

	isopen = 0
	fname = AlaamPluginFilename()
	p = AlaamPluginPath()
	if (p == "" | !fileexists(p)) {
		errprintf("{p 0 0 2}nwalaam: the compiled plugin {bf:%s} was not found. nwalaam estimates every model with this plugin and has no Mata fallback; reinstall the package with {bf:net install nwalaam, replace}.{p_end}\n", fname)
		exit(601)
	}
	origframe = st_framecurrent()
	stata("capture frame drop __nwalaam_native")
	stata("frame create __nwalaam_native")
	st_framecurrent("__nwalaam_native")
	isopen = 1
	st_addobs(4)
	(void) st_addvar("double", "v1")
	(void) st_addvar("double", "v2")
	q = char(34)
	stata("capture program drop alaamnative")
	stata("capture program alaamnative, plugin using(" + q + p + q + ")")
	rc = _stata("plugin call alaamnative v1 v2, " + q + "99 0" + q, 1)
	ver = .
	if (rc == 0) ver = st_data(1, 2)
	if (rc != 0 | ver != ALAAM_PROTOCOL_NEEDED()) {
		close()
		errprintf("{p 0 0 2}nwalaam: the installed plugin {bf:%s} is too old for this version of nwalaam (or it cannot be loaded). Reinstall the package with {bf:net install nwalaam, replace}.{p_end}\n", fname)
		exit(3495)
	}
}

void AlaamNat::close()
{
	if (isopen) {
		st_framecurrent(origframe)
		stata("capture frame drop __nwalaam_native")
		isopen = 0
	}
}

// adj: n x n tie matrix, directed: 0/1, y: outcome (. = missing), et/ep/ec: effect type, parameter, covariate number (1 x P),
// cov: one covariate per row (n columns)
void AlaamNat::setup(real matrix adj, real scalar directed, real colvector y, real rowvector et, real rowvector ep,
		real rowvector ec, real matrix cov)
{
	real scalar n, P, nt, k, nset, ndist, smode, hasfixed
	real matrix E, Sk
	real colvector yy, aux, fixed
	pointer scalar pp

	n = rows(adj)
	P = cols(et)
	E = adj :!= 0
	_diag(E, 0)
	E = E :& (adj :< .)
	if (!directed) E = uppertriangle(E, 0)
	E = selectindex(vec(E')) 
	nt = rows(E)
	E = (floor((E :- 1) / n) :+ 1), (mod(E :- 1, n) :+ 1)
	yy = editmissing(y, -1)
	// auxiliary objects of the fit, set by the caller as Mata externals: setting networks (n x n tie matrices), distance
	// matrices (n x n), the sampler (0 toggle, 1 zero-or-one) and the nodes whose outcome is held fixed (snowball)
	nset = 0
	while (findexternal("__nwalaam_set" + strofreal(nset + 1)) != NULL) nset++
	ndist = 0
	while (findexternal("__nwalaam_dist" + strofreal(ndist + 1)) != NULL) ndist++
	smode = 0
	pp = findexternal("__nwalaam_smode")
	if (pp != NULL) smode = *pp
	hasfixed = 0
	fixed = J(n, 1, 0)
	pp = findexternal("__nwalaam_fixed")
	if (pp != NULL) {
		fixed = *pp
		hasfixed = 1
	}
	blob = (ALAAM_PROTOCOL_NEEDED() \ n \ directed \ P \ nt \ rows(cov) \ nset \ ndist \ smode \ hasfixed \ J(6, 1, 0)) \ vec((et \ ep \ ec)) \ vec(E') \ yy
	if (rows(cov) > 0) blob = blob \ vec(cov')
	for (k = 1; k <= nset; k++) {
		pp = findexternal("__nwalaam_set" + strofreal(k))
		Sk = (*pp) :!= 0
		_diag(Sk, 0)
		Sk = Sk :& ((*pp) :< .)
		Sk = Sk :| Sk'
		Sk = uppertriangle(Sk, 0)
		E = selectindex(vec(Sk'))
		E = (floor((E :- 1) / n) :+ 1), (mod(E :- 1, n) :+ 1)
		blob = blob \ rows(E) \ vec(E')
	}
	for (k = 1; k <= ndist; k++) {
		pp = findexternal("__nwalaam_dist" + strofreal(k))
		blob = blob \ vec((*pp)')
	}
	if (hasfixed) blob = blob \ fixed
	L = rows(blob)
	nres = 0
	if (L > st_nobs()) st_addobs(L - st_nobs())
	st_store((1::L), 1, blob)
}

real colvector AlaamNat::call(string scalar arg, real scalar nrows)
{
	real scalar need

	need = max((L, nrows, 4))
	if (need > st_nobs()) st_addobs(need - st_nobs())
	stata("plugin call alaamnative v1 v2, " + char(34) + arg + char(34))
	return(st_data((1::nrows), 2))
}

/*
	Effects. spec: the effects as the user typed them ("density contagion gwactivity(0.69) ego(age)"), directed: 0/1,
	vars: the names of the covariates (rows of the covariate matrix are in the order of first use, returned in covnames).
	Returns et, ep, ec (1 x P) and the coefficient names; errors with r(198).
*/
void alaam_effects(string scalar spec, real scalar directed, real scalar twomode, real rowvector et, real rowvector ep, real rowvector ec,
		string rowvector names, string rowvector covnames, string rowvector setnames, string rowvector distnames)
{
	string rowvector bip
	string rowvector toks, un, dn
	real rowvector uc, dc
	string scalar tok, nm, arg, c, ps
	real scalar k, p, j, code, hasarg, par
	string matrix tbl

	// name, undirected code, directed code, kind (0 none, 1 parameter (alpha/beta), 2 covariate), default parameter
	tbl = J(0, 5, "")
	tbl = tbl \ ("density", "1", "1", "0", "" \
		"activity", "2", "0", "0", "" \
		"twostar", "3", "0", "0", "" \
		"threestar", "4", "0", "0", "" \
		"gwactivity", "5", "0", "1", "0.693147180559945" \
		"contagion", "6", "112", "0", "" \
		"alter2star1", "7", "0", "0", "" \
		"alter2star2", "8", "0", "0", "" \
		"partneractivity", "9", "0", "0", "" \
		"partnerpartner", "10", "0", "0", "")
	tbl = tbl \ ("triangle1", "11", "0", "0", "" \
		"triangle2", "12", "0", "0", "" \
		"triangle3", "13", "0", "0", "" \
		"gwcontagion", "14", "125", "1", "0.693147180559945" \
		"logcontagion", "15", "126", "0", "" \
		"powercontagion", "16", "127", "1", "2" \
		"ego", "17", "17", "2", "" \
		"alter", "18", "124", "2", "" \
		"same", "19", "128", "2", "" \
		"diff", "20", "129", "2", "")
	tbl = tbl \ ("samecontagion", "21", "130", "2", "" \
		"samepartneractivity", "22", "0", "2", "" \
		"diffpartneractivity", "23", "0", "2", "" \
		"sameindirect", "24", "0", "2", "" \
		"diffindirect", "25", "0", "2", "" \
		"alterbinarytwostar1", "26", "0", "2", "" \
		"alterbinarytwostar2", "27", "0", "2", "" \
		"fourcycle1", "28", "0", "0", "" \
		"fourcycle2", "29", "0", "0", "" \
		"settinghomophily", "40", "40", "3", "" \
		"geographichomophily", "41", "41", "4", "" \
		"contagiondist", "42", "0", "4", "" \
		"sendermatch", "0", "131", "2", "" \
		"sendermismatch", "0", "132", "2", "" \
		"receivermatch", "0", "133", "2", "")
	tbl = tbl \ ("receivermismatch", "0", "134", "2", "" \
		"reciprocitymatch", "0", "135", "2", "" \
		"reciprocitymismatch", "0", "136", "2", "" \
		"sender", "0", "102", "0", "" \
		"receiver", "0", "103", "0", "" \
		"reciprocity", "0", "104", "0", "" \
		"instar2", "0", "105", "0", "" \
		"instar3", "0", "106", "0", "" \
		"outstar2", "0", "107", "0", "" \
		"outstar3", "0", "108", "0", "")
	tbl = tbl \ ("mixedstar", "0", "109", "0", "" \
		"mixedsource", "0", "110", "0", "" \
		"mixedsink", "0", "111", "0", "" \
		"contagionrecip", "0", "113", "0", "" \
		"transitive1", "0", "114", "0", "" \
		"transitive3", "0", "115", "0", "" \
		"transitived1", "0", "116", "0", "" \
		"transitiveu1", "0", "117", "0", "" \
		"cyclic1", "0", "118", "0", "" \
		"cyclic3", "0", "119", "0", "")
	tbl = tbl \ ("alterinstar2", "0", "120", "0", "" \
		"alteroutstar2", "0", "121", "0", "" \
		"gwsender", "0", "122", "1", "0.693147180559945" \
		"gwreceiver", "0", "123", "1", "0.693147180559945")
	toks = tokens(subinstr(subinstr(spec, "(", " ( "), ")", " ) "))
	// re-join name and (argument): tokens() with "()" parse chars gives name, "(", arg, ")"
	et = J(1, 0, .); ep = et; ec = et; names = J(1, 0, ""); covnames = J(1, 0, ""); setnames = J(1, 0, ""); distnames = J(1, 0, "")
	// the effects defined for two-mode (bipartite) networks: those of the one-mode undirected list that do not need ties within a mode
	bip = ("density", "activity", "twostar", "threestar", "gwactivity", "alter2star1", "alter2star2", "fourcycle1", "fourcycle2",
		"ego", "alter", "samepartneractivity", "diffpartneractivity", "sameindirect", "diffindirect", "alterbinarytwostar1", "alterbinarytwostar2",
		"settinghomophily", "geographichomophily")
	k = 1
	while (k <= cols(toks)) {
		nm = strlower(strtrim(toks[k])); k++
		arg = ""; hasarg = 0
		if (k <= cols(toks)) {
			if (toks[k] == "(") {
				k++
				while (k <= cols(toks)) {
					if (toks[k] == ")") break
					arg = arg + toks[k]
					k++
				}
				k++
				hasarg = 1
			}
		}
		if (nm == "") continue
		j = selectindex(tbl[., 1] :== nm)
		if (rows(j) == 0) {
			errprintf("{p 0 0 2}nwalaam: unknown effect {bf:%s}.{p_end}\n", nm)
			exit(198)
		}
		if (twomode & !any(bip :== nm)) {
			errprintf("{p 0 0 2}nwalaam: effect {bf:%s} is not available for two-mode networks (the effects for two-mode networks are: %s).{p_end}\n", nm, invtokens(bip))
			exit(198)
		}
		code = strtoreal(tbl[j, (directed ? 3 : 2)])
		if (code == 0) {
			errprintf("{p 0 0 2}nwalaam: effect {bf:%s} is not available for %s networks.{p_end}\n", nm, (directed ? "directed" : "undirected"))
			exit(198)
		}
		p = strtoreal(tbl[j, 4])
		par = .
		c = ""
		if (p == 0 & hasarg) {
			errprintf("{p 0 0 2}nwalaam: effect {bf:%s} takes no argument.{p_end}\n", nm)
			exit(198)
		}
		if (p == 1) {
			par = strtoreal(tbl[j, 5])
			if (hasarg) {
				par = strtoreal(strtrim(arg))
				if (par >= . | par <= 0) {
					errprintf("{p 0 0 2}nwalaam: the argument of {bf:%s} must be a positive number.{p_end}\n", nm)
					exit(198)
				}
			}
		}
		if (p >= 3) {
			c = strtrim(arg)
			if (c == "") {
				errprintf("{p 0 0 2}nwalaam: effect {bf:%s} needs %s: %s(%s).{p_end}\n", nm, (p == 3 ? "a setting network" : "a distance matrix"), nm, (p == 3 ? "netname" : "matname"))
				exit(198)
			}
		}
		if (p == 2) {
			c = strtrim(arg)
			if (c == "") {
				errprintf("{p 0 0 2}nwalaam: effect {bf:%s} needs a covariate: %s(varname).{p_end}\n", nm, nm)
				exit(198)
			}
		}
		et = et, code
		ep = ep, (par >= . ? 0 : par)
		if (p == 2) {
			if (!any(covnames :== c)) covnames = covnames, c
			ec = ec, selectindex(covnames :== c)
			names = names, abbrev(nm + "_" + c, 32)
		}
		else if (p == 3) {
			if (!any(setnames :== c)) setnames = setnames, c
			ec = ec, selectindex(setnames :== c)
			names = names, abbrev(nm + "_" + c, 32)
		}
		else if (p == 4) {
			if (!any(distnames :== c)) distnames = distnames, c
			ec = ec, selectindex(distnames :== c)
			names = names, abbrev(nm + "_" + c, 32)
		}
		else if (p == 1 & hasarg) {
			ps = strtrim(strofreal(par, "%9.4g"))
			if (substr(ps, 1, 1) == ".") ps = "0" + ps
			ps = subinstr(ps, ".", "p")
			ec = ec, 0
			names = names, abbrev(nm + "_" + strtrim(ps), 32)
		}
		else {
			ec = ec, 0
			names = names, nm
		}
	}
	if (cols(et) == 0) {
		errprintf("{p 0 0 2}nwalaam: effects() must name at least one effect.{p_end}\n")
		exit(198)
	}
	if (cols(uniqrows(names')') != cols(names)) {
		errprintf("{p 0 0 2}nwalaam: an effect is listed twice.{p_end}\n")
		exit(198)
	}
}

// the nodes whose outcome is held at the observed value (outermost snowball zone), n x 1 of 0/1
real colvector alaam_fixed(real scalar n)
{
	pointer scalar pp

	pp = findexternal("__nwalaam_fixed")
	if (pp == NULL) return(J(n, 1, 0))
	return(*pp)
}

// logit without constant: Newton-Raphson; returns the coefficients (1 x P), missing when it does not converge, and the
// number of iterations in conv
real rowvector alaam_logit(real matrix X, real colvector y, real scalar conv)
{
	real colvector b, p, g, d, bn
	real matrix H
	real scalar it, ll, lln, step, k

	b = J(cols(X), 1, 0)
	ll = -rows(X) * ln(2)
	conv = 0
	for (it = 1; it <= 100; it++) {
		p = invlogit(X * b)
		g = X' * (y - p)
		H = cross(X, p :* (1 :- p), X)
		if (rank(H) < cols(X)) return(J(1, cols(X), .))
		d = invsym(H) * g
		step = 1
		for (k = 1; k <= 40; k++) {
			bn = b + step * d
			lln = sum(y :* ln(invlogit(X * bn)) + (1 :- y) :* ln(invlogit(-X * bn)))
			if (lln >= ll - 1e-12) break
			step = step / 2
		}
		b = bn
		ll = lln
		if (max(abs(step * d)) < 1e-9) {
			conv = it
			break
		}
		if (max(abs(b)) > 30) break
	}
	if (conv == 0 | max(abs(b)) > 30) return(J(1, cols(X), .))
	return(b')
}

// batch-means estimate of the asymptotic covariance matrix of the mean of the columns of X (as mcmcse's "bm" with "sqroot")
real matrix alaam_bmcov(real matrix X)
{
	real scalar n, b, a, k
	real rowvector mu
	real matrix Y, D

	n = rows(X)
	b = floor(sqrt(n))
	a = floor(n / b)
	mu = mean(X)
	Y = J(a, cols(X), .)
	for (k = 1; k <= a; k++) Y[k, .] = mean(X[((k - 1) * b + 1)..(k * b), .])
	D = Y :- mu
	return(b * cross(D, D) / (a - 1))
}

// change statistics of every node at the observed outcome (n x P)
real matrix AlaamDeltas(class AlaamNat scalar N, real scalar n, real scalar P)
{
	return(rowshape(N.call("10 " + strofreal(N.L, "%18.0f"), n * P), n))
}

real rowvector AlaamStats(class AlaamNat scalar N, real scalar P)
{
	return(N.call("11 " + strofreal(N.L, "%18.0f"), P)')
}

// the change statistics of every node (n x P matrix D) and the observed statistics (z) at the outcome y: used by the tests
void AlaamChangeStats(real matrix adj, real scalar directed, real colvector y, real rowvector et, real rowvector ep,
		real rowvector ec, real matrix cov, real matrix D, real rowvector z)
{
	class AlaamNat scalar N

	N.open()
	N.setup(adj, directed, y, et, ep, ec, cov)
	z = AlaamStats(N, cols(et))
	D = AlaamDeltas(N, rows(adj), cols(et))
	N.close()
}

/*
	AlaamFit(): the whole estimation. Posts the results to the Stata matrices mb (1 x P), mV (P x P) and mdiag (P x 6: observed
	statistic, pseudo-likelihood estimate, Algorithm S start, EE estimate, standard error, t-ratio) and returns through
	st_local: runs used, runs asked, mean acceptance rate. With mple = 1 only the pseudo-likelihood is computed.
*/
// appends the effects (xt, xp, xc, xm, xcn) that are not yet in the list (et, ep, ec, nm, cn); returns the number of effects added.
// Covariates are matched by name; xc refers to xcn.
real scalar alaam_append(real rowvector et, real rowvector ep, real rowvector ec, string rowvector nm, string rowvector cn,
		real rowvector xt, real rowvector xp, real rowvector xc, string rowvector xm, string rowvector xcn)
{
	real scalar k, add, c
	string scalar cname

	add = 0
	for (k = 1; k <= cols(xt); k++) {
		if (any(nm :== xm[k])) continue
		c = 0
		if (xc[k] > 0) {
			cname = xcn[xc[k]]
			if (!any(cn :== cname)) cn = cn, cname
			c = selectindex(cn :== cname)
		}
		et = et, xt[k]; ep = ep, xp[k]; ec = ec, c; nm = nm, xm[k]
		add++
	}
	return(add)
}

// lag-1 autocorrelation of every column
real rowvector alaam_acf1(real matrix X)
{
	real rowvector out
	real scalar j, m
	real colvector x

	out = J(1, cols(X), .)
	m = rows(X)
	for (j = 1; j <= cols(X); j++) {
		x = X[., j] :- mean(X[., j])
		if (cross(x, x) > 0) out[j] = cross(x[1..(m - 1)], x[2..m]) / cross(x, x)
	}
	return(out)
}

// nsim simulated outcome draws: statistics (vec = 0, nsim x P) or outcome vectors (vec = 1, nsim x n). The plugin model has
// the effects of the fit followed by the extra effects (theta = 0 for those)
real matrix AlaamSim(real matrix adj, real scalar directed, real colvector y, real rowvector et, real rowvector ep, real rowvector ec,
		real matrix cov, real rowvector theta, real scalar nsim, real scalar burn, real scalar interval, real scalar seed, real scalar vec)
{
	class AlaamNat scalar N
	real matrix out

	N.open()
	N.setup(adj, directed, y, et, ep, ec, cov)
	out = N.call((vec ? "32 " : "31 ") + strofreal(N.L, "%18.0f") + " " + strofreal(seed, "%18.0f") + " " + strofreal(nsim, "%18.0f") + " " +
		strofreal(burn, "%18.0f") + " " + strofreal(interval, "%18.0f") + " " + invtokens(strofreal(theta, "%25.17g")),
		nsim * (vec ? rows(adj) : cols(et)))
	N.close()
	return(rowshape(out', nsim))
}

// observed statistics of the effects et/ep/ec at the outcome y
real rowvector AlaamGofStats(real matrix adj, real scalar directed, real colvector y, real rowvector et, real rowvector ep, real rowvector ec,
		real matrix cov)
{
	class AlaamNat scalar N
	real rowvector z

	N.open()
	N.setup(adj, directed, y, et, ep, ec, cov)
	z = AlaamStats(N, cols(et))
	N.close()
	return(z)
}

// goodness of fit: obs (1 x K), sim (nsim x K) -> res (K x 4: observed, simulated mean, simulated sd, t-ratio) and the Mahalanobis distance
void AlaamGofTable(string scalar mobs, string scalar msim, string scalar mres, string scalar mmah)
{
	real rowvector o, mu, sd
	real matrix S, C
	real scalar d

	o = st_matrix(mobs)
	S = st_matrix(msim)
	mu = mean(S)
	sd = sqrt(diagonal(variance(S))')
	st_matrix(mres, (o', mu', sd', ((mu :- o) :/ sd)'))
	C = variance(S)
	d = sqrt(max(((o :- mu) * pinv(C) * (o :- mu)', 0)))
	st_matrix(mmah, d)
}

// per-run table: columns of RUNS are the estimate, convergence t-ratio, AC theta, AC dzA (P each)
void AlaamDiagTable(string scalar mrun, real scalar P, real scalar thr, string scalar mout)
{
	real matrix R, O
	real scalar k

	R = st_matrix(mrun)
	O = J(P, 7, .)
	for (k = 1; k <= P; k++) {
		O[k, 1] = mean(R[., k])
		O[k, 2] = rows(R) > 1 ? sqrt(variance(R[., k])) : .
		O[k, 3] = mean(R[., P + k])
		O[k, 4] = max(abs(R[., P + k]))
		O[k, 5] = sum(abs(R[., P + k]) :> thr)
		O[k, 6] = mean(R[., 2 * P + k])
		O[k, 7] = mean(R[., 3 * P + k])
	}
	st_matrix(mout, O)
}

// stochastic approximation (Snijders' Robbins-Monro): repeats the plugin's phases 1 to 3 from the previous estimate until every
// convergence t-ratio is below 0.1 in absolute value (at most 20 repetitions). The covariance of the estimate is the inverse of the
// batch-means covariance of the simulated statistics (the information) at the estimate.
void AlaamSA(class AlaamNat scalar N, real scalar P, real rowvector zobs, real rowvector theta0, real scalar mprop, real scalar gain,
		real scalar saiter, real scalar s0, real scalar seed, real matrix dg, string scalar mb, string scalar mV, string scalar mdiag,
		string scalar mrun, string scalar mtrace)
{
	real scalar rep, n3, ok, drop, nr
	real rowvector theta, tr
	real matrix R, DZ, C, TR
	real colvector sdv

	n3 = 4 * saiter
	theta = theta0
	ok = 0
	for (rep = 1; rep <= 20; rep++) {
		nr = (2 + n3) * P
		R = N.call("21 " + strofreal(N.L, "%18.0f") + " " + strofreal(seed + 7919 * rep + s0, "%18.0f") + " " + strofreal(mprop) + " " + strofreal(saiter) + " " +
			strofreal(saiter) + " 4 " + strofreal(n3) + " 0.25" + " " + invtokens(strofreal(theta, "%25.17g")), nr)
		R = rowshape(R', 2 + n3)
		theta = R[1, .]
		if (hasmissing(theta) | max(abs(theta)) > 1e10) {
			// diverged: restart from the pseudo-likelihood-free start with smaller gain
			theta = theta0
			gain = gain / 2
			continue
		}
		DZ = R[3::(2 + n3), .]
		drop = ceil(n3 / 10)
		DZ = DZ[(drop + 1)::n3, .]
		sdv = sqrt(diagonal(variance(DZ)))
		tr = (mean(DZ)' :/ sdv)'
		if (max(abs(tr)) < 0.1) {
			ok = 1
			break
		}
	}
	if (hasmissing(theta) | max(abs(theta)) > 1e10) {
		N.close()
		errprintf("{p 0 0 2}nwalaam: the stochastic approximation diverged (the model may be degenerate for this network; try fewer effects, or a smaller lrate()).{p_end}\n")
		exit(430)
	}
	C = alaam_bmcov(DZ)
	if (hasmissing(C) | rank(C) < P) {
		N.close()
		errprintf("{p 0 0 2}nwalaam: the simulated statistics are collinear at the estimate; the covariance matrix of the estimate cannot be computed.{p_end}\n")
		exit(430)
	}
	dg[., 4] = theta'
	dg[., 5] = sqrt(diagonal(invsym(C)))
	dg[., 6] = tr'
	dg[., 7] = J(P, 1, 0)
	st_matrix(mb, theta)
	st_matrix(mV, invsym(C))
	st_matrix(mdiag, dg)
	st_matrix(mrun, (theta, tr, alaam_acf1(theta \ theta), alaam_acf1(DZ)))
	TR = (J(rows(DZ), P, .), DZ, J(rows(DZ), 1, .))
	st_matrix(mtrace, TR)
	st_local("runs_used", strofreal(rep > 20 ? 20 : rep))
	st_local("accrate", ".")
	st_local("sa_converged", strofreal(ok))
}

void AlaamFit(real matrix adj, real scalar directed, real colvector y, real rowvector et, real rowvector ep, real rowvector ec,
		real matrix cov, real scalar mple, real scalar runs, real scalar M, real scalar thin, real scalar burn, real scalar mprop,
		real scalar lr, real scalar seed, real scalar s0, string scalar mstart, string scalar mb, string scalar mV, string scalar mdiag,
		string scalar mrun, string scalar mtrace, real scalar sa, real scalar saiter)
{
	class AlaamNat scalar N
	real scalar n, P, r, ns, used, tmax, conv, e, k, acc
	real matrix D, R, Th, Dz, Corr, allDz, dg
	real colvector ok, keep, tt, pp
	real rowvector zobs, bmple, theta0, est, se, tr
	real matrix EST, CV, MC, Vm, TC, Vmc, RUNS, TR1

	n = rows(adj)
	P = cols(et)
	N.open()
	N.setup(adj, directed, y, et, ep, ec, cov)
	zobs = AlaamStats(N, P)
	D = AlaamDeltas(N, n, P)
	ok = selectindex((y :< .) :& !alaam_fixed(n))
	bmple = alaam_logit(D[ok, .], y[ok], conv)
	dg = J(P, 7, .)
	dg[., 1] = zobs'
	dg[., 2] = bmple'
	if (mple) {
		st_matrix(mb, bmple)
		if (hasmissing(bmple)) st_matrix(mV, J(P, P, .))
		else {
			pp = invlogit(D[ok, .] * bmple')
			st_matrix(mV, invsym(cross(D[ok, .], pp :* (1 :- pp), D[ok, .])))
		}
		st_matrix(mdiag, dg)
		st_matrix(mrun, J(0, 4 * P, .))
		st_matrix(mtrace, J(0, 0, .))
		st_local("runs_used", "0")
		N.close()
		return
	}
	// Algorithm S: starting values
	theta0 = st_matrix(mstart)
	if (hasmissing(theta0)) theta0 = N.call("12 " + strofreal(N.L, "%18.0f") + " " + strofreal(s0, "%18.0f") + " 100 " + strofreal(mprop), P)'
	dg[., 3] = theta0'
	if (sa) {
		AlaamSA(N, P, zobs, theta0, mprop, lr, saiter, s0, seed, dg, mb, mV, mdiag, mrun, mtrace)
		N.close()
		return
	}
	ns = ceil(M / thin)
	EST = J(0, P, .); CV = J(0, P * P, .); MC = J(0, P * P, .)
	allDz = J(0, P, .)
	RUNS = J(0, 4 * P, .); TR1 = J(0, 0, .)
	used = 0; acc = 0
	for (r = 1; r <= runs; r++) {
		R = N.call("20 " + strofreal(N.L, "%18.0f") + " " + strofreal(seed + 1000003 * r, "%18.0f") + " " + strofreal(M) + " " +
			strofreal(mprop) + " " + strofreal(thin) + " " + strofreal(lr, "%18.12g") + " " + invtokens(strofreal(theta0, "%25.17g")),
			ns * (2 * P + 1))
		R = rowshape(R', ns)
		tt = (0::(ns - 1)) * thin
		keep = selectindex(tt :> burn)
		if (rows(keep) < 20) continue
		Th = R[keep, 1..P]; Dz = R[keep, (P + 1)..(2 * P)]
		acc = acc + mean(R[keep, 2 * P + 1])
		if (hasmissing(Th) | max(abs(Th)) > 1e10) continue
		TC = alaam_bmcov(Dz)
		if (hasmissing(TC) | rank(TC) < P) continue
		EST = EST \ mean(Th)
		CV = CV \ rowshape(invsym(TC), 1)
		MC = MC \ rowshape(alaam_bmcov(Th) / rows(Th), 1)
		allDz = allDz \ Dz
		// per run: estimate, convergence t-ratio, lag-1 autocorrelation of theta and of dzA
		RUNS = RUNS \ (mean(Th), mean(Dz) :/ sqrt(diagonal(variance(Dz))'), alaam_acf1(Th), alaam_acf1(Dz))
		if (used == 0) TR1 = R[., .]
		used++
	}
	if (used == 0) {
		N.close()
		errprintf("{p 0 0 2}nwalaam: none of the %f EE runs gave usable output (the model may be degenerate, or the covariate scales too large; try fewer or different effects, centred covariates, or more iterations).{p_end}\n", runs)
		exit(430)
	}
	// the estimate is the mean over the runs. Variance: the information part, inverse of the batch-means covariance of dzA (the
	// same in every run: averaged, not divided by the number of runs), plus the Monte-Carlo variance of the mean over the runs.
	est = mean(EST)
	Vmc = rowshape(mean(MC), P) / used
	Vm = rowshape(mean(CV), P) + Vmc
	se = sqrt(diagonal(Vm))'
	dg[., 7] = sqrt(diagonal(Vmc))
	tr = mean(allDz) :/ sqrt(diagonal(variance(allDz))')
	dg[., 4] = est'; dg[., 5] = se'; dg[., 6] = tr'
	st_matrix(mb, est)
	st_matrix(mV, Vm)
	st_matrix(mdiag, dg)
	st_matrix(mrun, RUNS)
	st_matrix(mtrace, TR1)
	st_local("runs_used", strofreal(used))
	st_local("accrate", strofreal(acc / used))
	N.close()
}

end
