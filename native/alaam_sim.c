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
	alaam_sim.c -- native (C) engine of nwalaam: change statistics of an autologistic actor attribute model (a fixed
	network, a binary outcome per node), the observed statistics, and the equilibrium-expectation (EE) sampler.

	The model is P(A = a) proportional to exp(sum_k theta_k z_k(a)). The change statistic delta_k(i) is the change of z_k
	when the outcome of node i goes from 0 to 1; it is evaluated with the outcome of i set to 0 and the CURRENT outcomes of
	all other nodes. The log odds of A_i = 1 given the rest is theta . delta(i).

	Calls: Mata writes the model once per fit as a block of numbers ("blob") into variable 1 of a frame; each call passes
	the code and the length of the blob in the argument string and gets its result in variable 2.

	Argument string: "<code> <L> ..."
	  code 99  protocol number (row 1 of variable 2)
	  code 10  change statistics of every node at the observed outcome: n x P, by rows (node 1 effect 1.. P, node 2 ..)
	  code 11  observed statistics: P rows (built by adding the outcome nodes one at a time, as the sum of the change
	           statistics)
	  code 12  Algorithm S (starting values): "12 L seed M1 mprop" -> P rows with theta
	  code 31  simulation: "31 L seed nsim burn interval theta_1 ... theta_P" -> nsim rows of P statistics (the statistics of the
	           simulated outcome vectors, updated incrementally from the observed ones; burn = proposals before the first
	           draw, interval = proposals between draws); the draw starts at the observed outcomes
	  code 32  like 31, but returns the nsim simulated outcome vectors (n values each, missing outcomes stay missing)
	  code 21  stochastic approximation: "21 L seed mprop n1 n2 nsub n3 gain theta_1 ... theta_P". Phase 1: n1 iterations at theta give
	           the derivative matrix (covariance of the statistics); phase 2: nsub subphases of n2*2^(k-1) iterations with the
	           Robbins-Monro step theta -= gain/2^(k-1) * D^-1 * dz (dz = statistics of the current outcome minus the observed ones);
	           the estimate is the mean of theta over the last subphase; phase 3: n3 iterations at the estimate. Result: P values of the
	           estimate, P values of the mean of dz in phase 3 (not used), then n3 rows of P values of dz in phase 3
	  code 20  EE run: "20 L seed M mprop thin lr theta_1 ... theta_P" -> for every saved iteration (t = 0, thin, 2 thin, ...):
	           theta (P values), dzA (P values), acceptance rate (1 value), written one iteration after the other
	Blob (doubles): [0] protocol [1] n [2] directed [3] P [4] number of ties listed [5] number of covariates (HDR = 16),
	then P x (effect type, parameter, covariate number (1-based, 0 = none)), then the ties (from, to), 1-based (an undirected
	tie once), then the outcome of every node (0, 1, or -1 = missing), then the covariates (n values each; missing = Stata missing).
	Effect types: see the table in change_stat() below (this table is the one of nwalaam.ado).
*/
#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define ALAAM_PROTOCOL 3
#define HDR 16

typedef struct {
	int n, directed, P, ncov;
	int *etype; double *epar; int *ecov;
	unsigned char *bits;        /* n x n: bit (i,j) = tie i -> j (symmetric when undirected) */
	int *oofs, *oadj;           /* out (or all) neighbours, CSR */
	int *iofs, *iadj;           /* in neighbours, CSR (directed only) */
	double **cov;
	signed char *y;
} Model;

static int isarc(const Model *g, int i, int j) {
	size_t k = (size_t)i * (size_t)g->n + (size_t)j;
	return (g->bits[k >> 3] >> (k & 7)) & 1;
}
static void setarc(Model *g, int i, int j) {
	size_t k = (size_t)i * (size_t)g->n + (size_t)j;
	g->bits[k >> 3] |= (unsigned char)(1 << (k & 7));
}
#define OUTDEG(g, i) ((g)->oofs[(i) + 1] - (g)->oofs[(i)])
#define INDEG(g, i) ((g)->iofs[(i) + 1] - (g)->iofs[(i)])
#define ONB(g, i, k) ((g)->oadj[(g)->oofs[(i)] + (k)])
#define INB(g, i, k) ((g)->iadj[(g)->iofs[(i)] + (k)])
#define MISS(x) ((x) > 1e300 || (x) != (x))

static int twopaths(const Model *g, int i, int v) {
	int k, c = 0;
	if (i == v) return 0;
	for (k = 0; k < OUTDEG(g, i); k++) if (isarc(g, v, ONB(g, i, k))) c++;
	return c;
}

/* contagion-type counts for the non-linear (gw, log, power) contagion statistics */
static double gw_like(const Model *g, const signed char *A, int i, int kind, double par) {
	int k, l, diplus = 0;
	double delta = 0.0;
	for (k = 0; k < OUTDEG(g, i); k++) {
		int j = ONB(g, i, k), djplus = 0;
		if (A[j] == 1) {
			diplus++;
			for (l = 0; l < OUTDEG(g, j); l++) if (A[ONB(g, j, l)] == 1) djplus++;
			if (kind == 0) delta += exp(-par * (djplus + 1)) - exp(-par * djplus);
			else if (kind == 1) delta += log((djplus + 2.0) / (djplus + 1.0));
			else delta += pow(djplus + 1.0, 1.0 / par) - pow((double)djplus, 1.0 / par);
		}
	}
	if (kind == 0) delta += exp(-par * diplus);
	else if (kind == 1) delta += log(diplus + 1.0);
	else delta += pow((double)diplus, 1.0 / par);
	return delta;
}

/*
	Change statistic of one effect for node i; A[i] must be 0 on entry.
	Undirected: 1 density 2 activity 3 twostar 4 threestar 5 gwactivity(a) 6 contagion 7 alter2star1 (partner activity two-path)
	  8 alter2star2 (indirect partner attribute) 9 partneractivity 10 partnerpartner 11 triangle1 12 triangle2 13 triangle3
	  14 gwcontagion(a) 15 logcontagion 16 powercontagion(b) 17 ego(x) 18 alter(x) 19 same(c) 20 diff(c) 21 samecontagion(c)
	Directed (ties i -> j): 1 density 17 ego(x) and
	  102 sender 103 receiver 104 reciprocity 105 instar2 106 instar3 107 outstar2 108 outstar3 109 mixedstar 110 mixedsource
	  111 mixedsink 112 contagion 113 contagionrecip 114 transitive1 115 transitive3 116 transitived1 117 transitiveu1
	  118 cyclic1 119 cyclic3 120 alterin2star2 121 alterout2star2 122 gwsender(a) 123 gwreceiver(a) 124 alter(x)
*/
static double change_stat(const Model *g, const signed char *A, int i, int e) {
	int type = g->etype[e], k, l, d = 0, di;
	double par = g->epar[e], delta = 0.0;
	const double *x = g->ecov[e] > 0 ? g->cov[g->ecov[e] - 1] : NULL;
	if (type == 1) return 1.0;
	if (type == 17) return MISS(x[i]) ? 0.0 : x[i];
	if (!g->directed) {
		di = OUTDEG(g, i);
		switch (type) {
		case 2: return di;
		case 3: return di > 1 ? di * (di - 1) / 2.0 : 0.0;
		case 4: return di > 2 ? di * (di - 1.0) * (di - 2.0) / 6.0 : 0.0;
		case 5: return exp(-par * di);
		case 6: for (k = 0; k < di; k++) if (A[ONB(g, i, k)] == 1) d++; return d;
		case 7: for (k = 0; k < di; k++) delta += OUTDEG(g, ONB(g, i, k)) - 1; return delta;
		case 8:
			for (k = 0; k < di; k++) {
				int u = ONB(g, i, k);
				for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && A[v] == 1) d++; }
			}
			return d;
		case 9:
			for (k = 0; k < di; k++) { int u = ONB(g, i, k); if (A[u] == 1) delta += di + OUTDEG(g, u) - 2; }
			return delta;
		case 10:
			for (k = 0; k < di; k++) {
				int u = ONB(g, i, k);
				if (A[u] == 1) {
					for (l = 0; l < OUTDEG(g, u); l++) if (A[ONB(g, u, l)] == 1) d += 2;
					for (l = 0; l < di; l++) { int v = ONB(g, i, l); if (A[v] == 1 && v != u) d++; }
				}
			}
			return d;
		case 11:
			if (di < 2) return 0.0;
			for (k = 0; k < di; k++) d += twopaths(g, i, ONB(g, i, k));
			return d / 2;
		case 12:
			if (di < 2) return 0.0;
			for (k = 0; k < di; k++) { int v = ONB(g, i, k); if (A[v] == 1) d += twopaths(g, i, v); }
			return d;
		case 13:
			if (di < 2) return 0.0;
			for (k = 0; k < di; k++) {
				int u = ONB(g, i, k);
				if (A[u] == 1) for (l = 0; l < k; l++) { int v = ONB(g, i, l); if (A[v] == 1 && isarc(g, u, v)) d++; }
			}
			return d;
		case 14: return gw_like(g, A, i, 0, par);
		case 15: return gw_like(g, A, i, 1, par);
		case 16: return gw_like(g, A, i, 2, par);
		case 18: for (k = 0; k < di; k++) { double v = x[ONB(g, i, k)]; if (!MISS(v)) delta += v; } return delta;
		case 19: case 20: case 21:
			if (MISS(x[i])) return 0.0;
			for (k = 0; k < di; k++) {
				int u = ONB(g, i, k);
				if (MISS(x[u])) continue;
				if (type == 19 && x[u] == x[i]) delta += 1;
				else if (type == 20 && x[u] != x[i]) delta += 1;
				else if (type == 21 && x[u] == x[i] && A[u] == 1) delta += 1;
			}
			return delta;
		}
		return 0.0;
	}
	/* directed */
	{
		int o = OUTDEG(g, i), in = INDEG(g, i);
		switch (type) {
		case 102: return o;
		case 103: return in;
		case 104: for (k = 0; k < o; k++) if (isarc(g, ONB(g, i, k), i)) d++; return d;
		case 105: return in > 1 ? in * (in - 1) / 2.0 : 0.0;
		case 106: return in > 2 ? in * (in - 1.0) * (in - 2.0) / 6.0 : 0.0;
		case 107: return o > 1 ? o * (o - 1) / 2.0 : 0.0;
		case 108: return o > 2 ? o * (o - 1.0) * (o - 2.0) / 6.0 : 0.0;
		case 109: for (k = 0; k < o; k++) if (isarc(g, ONB(g, i, k), i)) d++; return (double)in * o - d;
		case 110: for (k = 0; k < o; k++) { int u = ONB(g, i, k); delta += OUTDEG(g, u) - (isarc(g, u, i) ? 1 : 0); } return delta;
		case 111: for (k = 0; k < in; k++) { int u = INB(g, i, k); delta += INDEG(g, u) - (isarc(g, i, u) ? 1 : 0); } return delta;
		case 112:
			for (k = 0; k < o; k++) if (A[ONB(g, i, k)] == 1) d++;
			for (k = 0; k < in; k++) if (A[INB(g, i, k)] == 1) d++;
			return d;
		case 113: for (k = 0; k < o; k++) { int u = ONB(g, i, k); if (A[u] == 1 && isarc(g, u, i)) d++; } return d;
		case 114:
			for (k = 0; k < o; k++) { int u = ONB(g, i, k); for (l = 0; l < INDEG(g, u); l++) { int v = INB(g, u, l); if (v != i && isarc(g, v, i)) d++; } }
			return d;
		case 115:
			for (k = 0; k < o; k++) {
				int u = ONB(g, i, k);
				if (A[u] != 1) continue;
				for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, i, v) && A[v] == 1) d++; }
				for (l = 0; l < INDEG(g, u); l++) { int v = INB(g, u, l); if (v != i && isarc(g, v, i) && A[v] == 1) d++; }
			}
			for (k = 0; k < in; k++) {
				int u = INB(g, i, k);
				if (A[u] != 1) continue;
				for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, v, i) && A[v] == 1) d++; }
			}
			return d;
		case 116:
			for (k = 0; k < o; k++) { int u = ONB(g, i, k); for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, i, v)) d++; } }
			return d;
		case 117:
			for (k = 0; k < in; k++) { int u = INB(g, i, k); for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, v, i)) d++; } }
			return d;
		case 118:
			for (k = 0; k < o; k++) { int u = ONB(g, i, k); for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, v, i)) d++; } }
			return d;
		case 119:
			for (k = 0; k < o; k++) {
				int u = ONB(g, i, k);
				if (A[u] != 1) continue;
				for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && isarc(g, v, i) && A[v] == 1) d++; }
			}
			return d;
		case 120:
			for (k = 0; k < in; k++) { int u = INB(g, i, k); for (l = 0; l < OUTDEG(g, u); l++) { int v = ONB(g, u, l); if (v != i && A[v] == 1) d++; } }
			return d;
		case 121:
			for (k = 0; k < o; k++) { int u = ONB(g, i, k); for (l = 0; l < INDEG(g, u); l++) { int v = INB(g, u, l); if (v != i && A[v] == 1) d++; } }
			return d;
		case 122: return exp(-par * o);
		case 123: return exp(-par * in);
		case 124:
			for (k = 0; k < o; k++) { double v = x[ONB(g, i, k)]; if (!MISS(v)) delta += v; }
			for (k = 0; k < in; k++) { double v = x[INB(g, i, k)]; if (!MISS(v)) delta += v; }
			return delta;
		}
	}
	return 0.0;
}

/* xoshiro256** */
typedef struct { uint64_t s[4]; } Rng;
static uint64_t splitmix(uint64_t *x) {
	uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}
static void rng_seed(Rng *r, uint64_t seed) { int k; for (k = 0; k < 4; k++) r->s[k] = splitmix(&seed); }
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rng_next(Rng *r) {
	uint64_t *s = r->s, res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
	s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
	return res;
}
static double rng_unif(Rng *r) { return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }

static Rng *cur_rng = NULL;   /* the generator of the running call */

/* the toggle sampler: mprop proposals; S1 and S0 are the sums of the change statistics of accepted moves to 1 and to 0 */
static double sampler(const Model *g, signed char *A, const int *elig, int ne, const double *theta, int perform,
		int mprop, double *S1, double *S0, double *cs) {
	int m, e, acc = 0;
	for (e = 0; e < g->P; e++) { S1[e] = 0.0; S0[e] = 0.0; }
	for (m = 0; m < mprop; m++) {
		int i = elig[(int)(rng_unif(cur_rng) * ne)];
		int to0 = (A[i] == 1);
		double tot = 0.0, u;
		if (to0) A[i] = 0;
		for (e = 0; e < g->P; e++) { cs[e] = change_stat(g, A, i, e); tot += theta[e] * cs[e]; }
		if (to0) tot = -tot;
		u = rng_unif(cur_rng);
		if (u < exp(tot)) {
			acc++;
			if (perform) { if (!to0) A[i] = 1; }
			else if (to0) A[i] = 1;
			for (e = 0; e < g->P; e++) { if (to0) S0[e] += cs[e]; else S1[e] += cs[e]; }
		}
		else if (to0) A[i] = 1;
	}
	return (double)acc / mprop;
}

/* inverse of the P x P matrix a (row-major) by Gauss-Jordan with partial pivoting; returns 0 when singular */
static int invert(const double *a, double *inv, int P) {
	int i, j, k, piv;
	double *m = malloc((size_t)P * 2 * P * sizeof(double));
	int W = 2 * P, ok = 1;
	for (i = 0; i < P; i++) for (j = 0; j < P; j++) { m[i * W + j] = a[i * P + j]; m[i * W + P + j] = (i == j) ? 1.0 : 0.0; }
	for (k = 0; k < P && ok; k++) {
		double best = fabs(m[k * W + k]), d;
		piv = k;
		for (i = k + 1; i < P; i++) if (fabs(m[i * W + k]) > best) { best = fabs(m[i * W + k]); piv = i; }
		if (best < 1e-12) { ok = 0; break; }
		if (piv != k) for (j = 0; j < W; j++) { double t = m[k * W + j]; m[k * W + j] = m[piv * W + j]; m[piv * W + j] = t; }
		d = m[k * W + k];
		for (j = 0; j < W; j++) m[k * W + j] /= d;
		for (i = 0; i < P; i++) if (i != k) {
			double f = m[i * W + k];
			if (f != 0.0) for (j = 0; j < W; j++) m[i * W + j] -= f * m[k * W + j];
		}
	}
	if (ok) for (i = 0; i < P; i++) for (j = 0; j < P; j++) inv[i * P + j] = m[i * W + P + j];
	free(m);
	return ok;
}
static void model_free(Model *g) {
	int k;
	free(g->etype); free(g->epar); free(g->ecov); free(g->bits); free(g->oofs); free(g->oadj); free(g->iofs); free(g->iadj); free(g->y);
	if (g->cov) { for (k = 0; k < g->ncov; k++) free(g->cov[k]); free(g->cov); }
}

static int model_read(Model *g, long L) {
	long pos = 0, k, nt;
	int n, P, i, e;
	double h[HDR];
	int *ea = NULL, *eb = NULL, *cnt = NULL;
	memset(g, 0, sizeof(*g));
	if (L < HDR) return 198;
	for (k = 0; k < HDR; k++) { ST_double v; SF_vdata(1, (ST_int)(k + 1), &v); h[k] = v; }
	if ((int)h[0] != ALAAM_PROTOCOL) return 198;
	n = (int)h[1]; g->n = n; g->directed = (int)h[2]; P = (int)h[3]; g->P = P; nt = (long)h[4]; g->ncov = (int)h[5];
	if (n < 1 || P < 1 || P > 1000 || nt < 0) return 198;
	pos = HDR;
	g->etype = calloc(P, sizeof(int)); g->epar = calloc(P, sizeof(double)); g->ecov = calloc(P, sizeof(int));
	for (e = 0; e < P; e++) {
		ST_double v;
		SF_vdata(1, (ST_int)(++pos), &v); g->etype[e] = (int)v;
		SF_vdata(1, (ST_int)(++pos), &v); g->epar[e] = v;
		SF_vdata(1, (ST_int)(++pos), &v); g->ecov[e] = (int)v;
	}
	g->bits = calloc(((size_t)n * (size_t)n + 7) / 8 + 1, 1);
	ea = malloc((nt + 1) * sizeof(int)); eb = malloc((nt + 1) * sizeof(int));
	if (!g->bits || !ea || !eb) return 909;
	for (k = 0; k < nt; k++) {
		ST_double a, b;
		SF_vdata(1, (ST_int)(++pos), &a); SF_vdata(1, (ST_int)(++pos), &b);
		ea[k] = (int)a - 1; eb[k] = (int)b - 1;
		if (ea[k] < 0 || eb[k] < 0 || ea[k] >= n || eb[k] >= n) { free(ea); free(eb); return 198; }
		if (ea[k] == eb[k]) continue;
		setarc(g, ea[k], eb[k]);
		if (!g->directed) setarc(g, eb[k], ea[k]);
	}
	free(ea); free(eb);
	g->oofs = calloc(n + 1, sizeof(int)); g->iofs = calloc(n + 1, sizeof(int));
	cnt = calloc(n + 1, sizeof(int));
	for (i = 0; i < n; i++) {
		int j, co = 0, ci = 0;
		for (j = 0; j < n; j++) { if (isarc(g, i, j)) co++; if (g->directed && isarc(g, j, i)) ci++; }
		g->oofs[i + 1] = g->oofs[i] + co; g->iofs[i + 1] = g->iofs[i] + ci;
	}
	g->oadj = malloc((g->oofs[n] + 1) * sizeof(int)); g->iadj = malloc((g->iofs[n] + 1) * sizeof(int));
	for (i = 0; i < n; i++) {
		int j, po = g->oofs[i], pi = g->iofs[i];
		for (j = 0; j < n; j++) { if (isarc(g, i, j)) g->oadj[po++] = j; if (g->directed && isarc(g, j, i)) g->iadj[pi++] = j; }
	}
	free(cnt);
	g->y = malloc(n);
	for (i = 0; i < n; i++) { ST_double v; SF_vdata(1, (ST_int)(++pos), &v); g->y[i] = (signed char)(v < -0.5 ? -1 : (v > 0.5 ? 1 : 0)); }
	g->cov = calloc(g->ncov + 1, sizeof(double *));
	for (k = 0; k < g->ncov; k++) {
		g->cov[k] = malloc(n * sizeof(double));
		for (i = 0; i < n; i++) { ST_double v; SF_vdata(1, (ST_int)(++pos), &v); g->cov[k][i] = v; }
	}
	for (e = 0; e < P; e++) if (g->ecov[e] > g->ncov) return 198;
	if (pos > L) return 198;
	return 0;
}

static int eligible(const Model *g, int **elig) {
	int i, ne = 0;
	*elig = malloc(g->n * sizeof(int));
	for (i = 0; i < g->n; i++) if (g->y[i] >= 0) (*elig)[ne++] = i;
	return ne;
}

STDLL stata_call(int argc, char *argv[]) {
	char *buf, *tok;
	int code, rc = 0;
	long L = 0;
	Model g;
	if (argc < 1) { SF_error("alaam_sim: missing argument string\n"); return 198; }
	buf = malloc(strlen(argv[0]) + 1);
	if (!buf) return 909;
	strcpy(buf, argv[0]);
	tok = strtok(buf, " \t");
	code = tok ? (int)atof(tok) : 0;
	if (code == 99) { free(buf); SF_vstore(2, 1, (ST_double)ALAAM_PROTOCOL); return 0; }
	tok = strtok(NULL, " \t");
	L = tok ? (long)atof(tok) : 0;
	rc = model_read(&g, L);
	if (rc) { free(buf); SF_error("alaam_sim: the model block is invalid (plugin and package versions do not match?)\n"); return rc; }
	{
		int n = g.n, P = g.P, i, e, ne, *elig = NULL;
		signed char *A = malloc(n);
		double *cs = calloc(P, sizeof(double)), *S1 = calloc(P, sizeof(double)), *S0 = calloc(P, sizeof(double));
		memcpy(A, g.y, n);
		ne = eligible(&g, &elig);
		if (code == 10) {
			/* missing outcomes count as 0 in every statistic */
			signed char *W = malloc(n); int j;
			for (j = 0; j < n; j++) W[j] = (A[j] == 1) ? 1 : 0;
			for (i = 0; i < n; i++) {
				W[i] = 0;
				for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + (long)i * P + e), change_stat(&g, W, i, e));
				W[i] = (A[i] == 1) ? 1 : 0;
			}
			free(W);
		}
		else if (code == 11) {
			signed char *W = calloc(n, 1);
			double *z = calloc(P, sizeof(double));
			for (i = 0; i < n; i++) if (A[i] == 1) {
				for (e = 0; e < P; e++) z[e] += change_stat(&g, W, i, e);
				W[i] = 1;
			}
			for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + e), z[e]);
			free(W); free(z);
		}
		else if (code == 31 || code == 32) {
			char *t1 = strtok(NULL, " \t"), *t2 = strtok(NULL, " \t"), *t3 = strtok(NULL, " \t"), *t4 = strtok(NULL, " \t");
			Rng rng;
			double *theta = calloc(P, sizeof(double)), *z = calloc(P, sizeof(double));
			signed char *W = malloc(n), *Z0 = calloc(n, 1);
			long nsim, burn, interval, s, row = 0;
			if (!t1 || !t2 || !t3 || !t4) { SF_error("alaam_sim: short argument string\n"); rc = 198; goto done2; }
			rng_seed(&rng, (uint64_t)atof(t1)); cur_rng = &rng;
			nsim = (long)atof(t2); burn = (long)atof(t3); interval = (long)atof(t4);
			for (e = 0; e < P; e++) { char *x = strtok(NULL, " \t"); if (!x) { rc = 198; goto done2; } theta[e] = atof(x); }
			if (ne < 1) { SF_error("alaam_sim: no node with an observed outcome\n"); rc = 198; goto done2; }
			for (i = 0; i < n; i++) W[i] = (A[i] == 1) ? 1 : 0;
			/* statistics of the starting outcomes */
			for (i = 0; i < n; i++) if (A[i] == 1) {
				for (e = 0; e < P; e++) z[e] += change_stat(&g, Z0, i, e);
				Z0[i] = 1;
			}
			if (burn > 0) { sampler(&g, W, elig, ne, theta, 1, (int)burn, S1, S0, cs); for (e = 0; e < P; e++) z[e] += S1[e] - S0[e]; }
			for (s = 0; s < nsim; s++) {
				if (interval > 0) { sampler(&g, W, elig, ne, theta, 1, (int)interval, S1, S0, cs); for (e = 0; e < P; e++) z[e] += S1[e] - S0[e]; }
				if (code == 31) { for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + row + e), z[e]); row += P; }
				else { for (i = 0; i < n; i++) SF_vstore(2, (ST_int)(1 + row + i), A[i] < 0 ? SV_missval : (ST_double)W[i]); row += n; }
			}
		done2:
			free(theta); free(z); free(W); free(Z0);
		}
		else if (code == 21) {
			char *t1 = strtok(NULL, " \t"), *t2 = strtok(NULL, " \t"), *t3 = strtok(NULL, " \t"), *t4 = strtok(NULL, " \t");
			char *t5 = strtok(NULL, " \t"), *t6 = strtok(NULL, " \t"), *t7 = strtok(NULL, " \t");
			Rng rng;
			double *theta = calloc(P, sizeof(double)), *dz = calloc(P, sizeof(double)), *mu = calloc(P, sizeof(double));
			double *stv = calloc(P, sizeof(double)), *Cm = calloc((size_t)P * P, sizeof(double)), *Di = calloc((size_t)P * P, sizeof(double)), *avg = calloc(P, sizeof(double));
			double gain;
			int mprop, n1, n2, nsub, n3, k, it, j;
			signed char *W = malloc(n);
			if (!t1 || !t2 || !t3 || !t4 || !t5 || !t6 || !t7) { SF_error("alaam_sim: short argument string\n"); rc = 198; goto done3; }
			rng_seed(&rng, (uint64_t)atof(t1)); cur_rng = &rng;
			mprop = (int)atof(t2); n1 = (int)atof(t3); n2 = (int)atof(t4); nsub = (int)atof(t5); n3 = (int)atof(t6); gain = atof(t7);
			for (e = 0; e < P; e++) { char *x = strtok(NULL, " \t"); if (!x) { rc = 198; goto done3; } theta[e] = atof(x); }
			if (ne < 1) { SF_error("alaam_sim: no node with an observed outcome\n"); rc = 198; goto done3; }
			for (i = 0; i < n; i++) W[i] = (A[i] == 1) ? 1 : 0;
			/* phase 1: derivative matrix = covariance of the statistics at theta (the exponential family identity dE z / dtheta = Cov z) */
			for (it = 0; it < n1; it++) {
				sampler(&g, W, elig, ne, theta, 1, mprop, S1, S0, cs);
				for (e = 0; e < P; e++) dz[e] += S1[e] - S0[e];
				for (e = 0; e < P; e++) { mu[e] += dz[e]; for (j = 0; j < P; j++) Cm[e * P + j] += dz[e] * dz[j]; }
			}
			for (e = 0; e < P; e++) for (j = 0; j < P; j++) Cm[e * P + j] = Cm[e * P + j] / n1 - (mu[e] / n1) * (mu[j] / n1);
			/* Levenberg damping: the covariance is poorly estimated, and can be nearly singular, when the chain sits in one state */
			for (e = 0; e < P; e++) for (j = 0; j < P; j++) if (e != j) Cm[e * P + j] *= 0.7;
			for (e = 0; e < P; e++) if (Cm[e * P + e] < 1e-3) Cm[e * P + e] = 1e-3;
			if (!invert(Cm, Di, P)) {
				/* singular: use the diagonal */
				for (e = 0; e < P; e++) for (j = 0; j < P; j++) Di[e * P + j] = 0.0;
				for (e = 0; e < P; e++) Di[e * P + e] = Cm[e * P + e] > 1e-12 ? 1.0 / Cm[e * P + e] : 1.0;
			}
			/* phase 2 */
			for (k = 0; k < nsub; k++) {
				int nk = n2 << k;
				double a = gain / (double)(1 << k);
				for (e = 0; e < P; e++) avg[e] = 0.0;
				for (it = 0; it < nk; it++) {
					sampler(&g, W, elig, ne, theta, 1, mprop, S1, S0, cs);
					for (e = 0; e < P; e++) dz[e] += S1[e] - S0[e];
					{
						double mx = 0.0, sc = 1.0, *st = stv;
						for (e = 0; e < P; e++) {
							st[e] = 0.0;
							for (j = 0; j < P; j++) st[e] += Di[e * P + j] * dz[j];
							st[e] *= a;
							if (fabs(st[e]) > mx) mx = fabs(st[e]);
						}
						if (mx > 0.25) sc = 0.25 / mx;
						for (e = 0; e < P; e++) { theta[e] -= sc * st[e]; avg[e] += theta[e]; }
					}
				}
				if (k == nsub - 1) for (e = 0; e < P; e++) theta[e] = avg[e] / nk;
			}
			/* phase 3: at the estimate, from the state reached (observed statistics are tracked through dz) */
			for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + e), theta[e]);
			for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + P + e), 0.0);
			for (it = 0; it < n3; it++) {
				sampler(&g, W, elig, ne, theta, 1, mprop, S1, S0, cs);
				for (e = 0; e < P; e++) { dz[e] += S1[e] - S0[e]; SF_vstore(2, (ST_int)(1 + (long)(2 + it) * P + e), dz[e]); }
			}
		done3:
			free(theta); free(dz); free(mu); free(stv); free(Cm); free(Di); free(avg); free(W);
		}
		else if (code == 12 || code == 20) {
			char *t1 = strtok(NULL, " \t"), *t2 = strtok(NULL, " \t"), *t3 = strtok(NULL, " \t");
			Rng rng;
			double *theta = calloc(P, sizeof(double)), *dzA = calloc(P, sizeof(double));
			int M, mprop, thin = 1, t;
			double lr = 0.01;
			signed char *W = malloc(n);
			if (!t1 || !t2 || !t3) { SF_error("alaam_sim: short argument string\n"); rc = 198; goto done; }
			rng_seed(&rng, (uint64_t)atof(t1)); cur_rng = &rng;
			M = (int)atof(t2); mprop = (int)atof(t3);
			for (i = 0; i < n; i++) W[i] = (A[i] == 1) ? 1 : 0;
			if (ne < 1) { SF_error("alaam_sim: no node with an observed outcome\n"); rc = 198; goto done; }
			if (code == 12) {
				double *D0 = calloc(P, sizeof(double));
				for (t = 0; t < M; t++) {
					sampler(&g, W, elig, ne, theta, 0, mprop, S1, S0, cs);
					for (e = 0; e < P; e++) {
						double dz = S0[e] - S1[e], sc = S0[e] + S1[e], step, da = 0.0;
						D0[e] += dz * dz;
						if (sc != 0.0) da = 0.1 / (sc * sc);
						step = (dz > 0 ? 1.0 : (dz < 0 ? -1.0 : 0.0)) * da * dz * dz;
						if (step > 0.1) step = 0.1;
						theta[e] += step;
					}
				}
				free(D0);
				for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + e), theta[e]);
			}
			else {
				char *tt = strtok(NULL, " \t"), *tl = strtok(NULL, " \t");
				long row = 0;
				if (!tt || !tl) { SF_error("alaam_sim: short argument string\n"); rc = 198; goto done; }
				thin = (int)atof(tt); lr = atof(tl);
				for (e = 0; e < P; e++) { char *x = strtok(NULL, " \t"); if (!x) { rc = 198; goto done; } theta[e] = atof(x); }
				for (t = 0; t < M; t++) {
					double acc = sampler(&g, W, elig, ne, theta, 1, mprop, S1, S0, cs);
					for (e = 0; e < P; e++) {
						double step;
						dzA[e] += S1[e] - S0[e];
						step = fabs(theta[e]) > 0.01 ? fabs(theta[e]) : 0.01;
						step *= lr * (dzA[e] > 0 ? 1.0 : (dzA[e] < 0 ? -1.0 : 0.0));
						theta[e] -= step;
					}
					if (t % thin == 0) {
						for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + row + e), theta[e]);
						for (e = 0; e < P; e++) SF_vstore(2, (ST_int)(1 + row + P + e), dzA[e]);
						SF_vstore(2, (ST_int)(1 + row + 2 * P), acc);
						row += 2 * P + 1;
					}
				}
			}
		done:
			free(theta); free(dzA); free(W);
		}
		else { SF_error("alaam_sim: unknown algorithm code (plugin and package versions do not match)\n"); rc = 198; }
		free(A); free(cs); free(S1); free(S0); free(elig);
	}
	model_free(&g);
	free(buf);
	return rc;
}
