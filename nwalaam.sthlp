{smcl}
{* *! version 0.2.0  18nov2026 author: Thomas Grund}{...}
{marker topic}
{title:Title}

{p2colset 9 17 24 2}{...}
{p2col :nwalaam  {hline 2}}Autologistic actor attribute model (ALAAM) for a binary node outcome on a fixed network{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:nwalaam} {it:outcome} [{it:netname}]{cmd:,} {opth effects(string)} [{it:options}]

{p 8 16 2}
{cmd:nwalaam}  [{cmd:,} {opt ci} {opt noci} {opt level(#)} {it:display_options}]  (replay after a fit)

{p 8 16 2}
{cmd:nwalaam version}

{pstd}
{it:outcome} is a numeric variable of the current dataset with the values 0 and 1 (missing values are allowed: these actors are
never changed by the sampler and count as 0 in every statistic), one row per actor. When the dataset carries the node labels
(variable {cmd:_nwnode}, as after {helpb nwload} with {cmd:xvars}) the rows are matched to the actors by label; otherwise the
dataset must have one row per actor, in actor order. {it:netname} is the network (default: the current network); it may be
undirected or directed, one-mode. Valued networks are used as binary (a tie where the value is not 0).

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt eff:ects(list)}}the effects of the model (required); see below{p_end}
{synopt :{opt mple}}maximum pseudo-likelihood estimate only (a logit of the outcome on the change statistics); no standard errors from sampling{p_end}
{synopt :{opt runs(#)}}independent EE runs; default {cmd:30}{p_end}
{synopt :{opt iter:ations(#)}}iterations of every EE run, each of {cmd:proposals()} toggle proposals; default {cmd:50000}{p_end}
{synopt :{opt thin(#)}}every {it:#}th iteration is stored; default {cmd:100}{p_end}
{synopt :{opt burn:in(#)}}iterations discarded at the start of every run; default {cmd:1000}{p_end}
{synopt :{opt prop:osals(#)}}proposals per iteration; default {cmd:1000}{p_end}
{synopt :{opt lrate(#)}}learning rate of the EE update; default {cmd:0.01}{p_end}
{synopt :{opt start(matname)}}1 x {it:P} starting values (default: Algorithm S, 100 iterations){p_end}
{synopt :{opt seed(#)}}random-number seed of the runs (default: drawn from Stata's generator, so {cmd:set seed} makes a fit reproducible){p_end}
{synopt :{opt ci}, {opt noci}, {opt level(#)}, {it:display_options}}display of the coefficient table; see below{p_end}
{synoptline}


{marker desc}{...}
{title:Description}

{pstd}
{cmd:nwalaam} fits the autologistic actor attribute model for a binary outcome A on the nodes of a fixed network:
P(A = a) is proportional to exp(sum_k theta_k z_k(a)), where the statistics z_k(a) count configurations of the outcome on the network
(for example the nodes with the outcome that have a neighbour with the outcome, {cmd:contagion}). The log odds of A_i = 1 given the
others is the sum of theta_k times the change statistic of node i (evaluated at A_i = 0, with the current outcomes of the others).
The model follows ALAAMEE (Stivala, Wang and Lomi) and MPNet.

{pstd}
The default estimator is the equilibrium-expectation (EE) algorithm of Borisenko, Byshkin and Lomi in the form used by ALAAMEE: every
iteration proposes {cmd:proposals()} toggles of the outcome of a random node, accepts with the Metropolis probability, and moves
every parameter by a relative step {cmd:lrate()} against the sign of the accumulated difference between the simulated and the observed
statistics. Starting values come from Algorithm S. {cmd:runs()} independent runs are made; the estimate is their mean. The covariance
matrix is the inverse of the batch-means covariance of the accumulated differences (the information part, averaged over the runs) plus
the Monte-Carlo variance of the mean of the runs. The column {cmd:e(tratio)} holds the convergence t-ratios (mean over standard
deviation of the accumulated differences): a ratio above 0.3 in absolute value is flagged, as in ALAAMEE.

{pstd}
All computing is done in a compiled plugin; there is no slower fallback (error 601: plugin missing; 3495: plugin too old).
The change statistics of the plugin were compared with the definitions of the statistics.


{marker effects}{...}
{title:Effects}

{pstd}Effects are listed separated by blanks; those with an argument are written with the argument in parentheses.

{pstd}{bf:Undirected networks}{p_end}
{p2colset 9 33 35 2}{...}
{p2col :{cmd:density}}number of nodes with the outcome (intercept){p_end}
{p2col :{cmd:activity}}sum of the degrees of the nodes with the outcome{p_end}
{p2col :{cmd:twostar}, {cmd:threestar}}sum over the nodes with the outcome of C(degree, 2), C(degree, 3){p_end}
{p2col :{cmd:gwactivity}[{cmd:(}{it:a}{cmd:)}]}geometrically weighted activity, sum of exp(-a degree); {it:a} > 0, default log 2{p_end}
{p2col :{cmd:contagion}}ties between two nodes with the outcome{p_end}
{p2col :{cmd:alter2star1}}sum over the nodes with the outcome of (degree of the neighbour - 1) (partner activity two-path){p_end}
{p2col :{cmd:alter2star2}}two-paths between two nodes with the outcome (indirect partner attribute){p_end}
{p2col :{cmd:partneractivity}}ties between two nodes with the outcome, weighted by degree sum - 2{p_end}
{p2col :{cmd:partnerpartner}}two times the number of two-paths with the outcome at all three nodes{p_end}
{p2col :{cmd:triangle1}, {cmd:triangle2}, {cmd:triangle3}}triangles: weighted by the number of nodes with the outcome, by the pairs of them, only those with all three{p_end}
{p2col :{cmd:gwcontagion}[{cmd:(}{it:a}{cmd:)}]}sum over nodes with the outcome of exp(-a n) with n the number of neighbours with the outcome{p_end}
{p2col :{cmd:logcontagion}}sum over nodes with the outcome of log(n + 1){p_end}
{p2col :{cmd:powercontagion}[{cmd:(}{it:b}{cmd:)}]}sum over nodes with the outcome of n^(1/b); default {it:b} = 2{p_end}
{p2col :{cmd:ego(}{it:var}{cmd:)}}sum of the numeric covariate over the nodes with the outcome (centre the covariate){p_end}
{p2col :{cmd:alter(}{it:var}{cmd:)}}sum over the nodes with the outcome of the covariate of their neighbours{p_end}
{p2col :{cmd:same(}{it:var}{cmd:)}, {cmd:diff(}{it:var}{cmd:)}}neighbours with the same / a different category{p_end}
{p2col :{cmd:samecontagion(}{it:var}{cmd:)}}ties between two nodes with the outcome and the same category{p_end}
{p2colreset}{...}

{pstd}{bf:Directed networks}: {cmd:density}, {cmd:sender}, {cmd:receiver}, {cmd:reciprocity}, {cmd:instar2}, {cmd:instar3}, {cmd:outstar2},
{cmd:outstar3}, {cmd:mixedstar}, {cmd:mixedsource}, {cmd:mixedsink}, {cmd:contagion}, {cmd:contagionrecip}, {cmd:transitive1},
{cmd:transitive3}, {cmd:transitived1}, {cmd:transitiveu1}, {cmd:cyclic1}, {cmd:cyclic3}, {cmd:alterinstar2}, {cmd:alteroutstar2},
{cmd:gwsender}[{cmd:(}{it:a}{cmd:)}], {cmd:gwreceiver}[{cmd:(}{it:a}{cmd:)}], {cmd:ego(}{it:var}{cmd:)}, {cmd:alter(}{it:var}{cmd:)}.
The definitions are those of ALAAMEE (Stivala, Wang and Lomi).

{pstd}Not available yet: two-mode networks, setting and geographic homophily, distance effects, binary-covariate two-stars.
Geometrically weighted statistics are the remedy for near-degeneracy of the activity and star statistics on larger networks
(Stivala 2023); read their sign with care.


{marker display}{...}
{title:Display}

{pstd}The coefficient names are shown in full. At a line size below 100 the confidence intervals are left out so that names of up to 32
characters fit; {opt ci} (or an explicit {opt level(#)}) shows them (long names are then cut), {opt noci} leaves them out. {opt ci} and {opt noci}
together are an error. Other display options of Stata's estimation commands ({opt nopvalues}, {opt cformat()}, ...) are passed on.


{marker results}{...}
{title:Stored results}

{pstd}{cmd:nwalaam} stores in {cmd:e()}: {cmd:e(b)}, {cmd:e(V)}; {cmd:e(obs_stats)} (observed statistics), {cmd:e(b_mple)} (pseudo-likelihood
estimate; missing when it does not exist, for example under perfect separation), {cmd:e(b_start)}, {cmd:e(tratio)}, {cmd:e(diag)};
scalars {cmd:e(N)}, {cmd:e(n_ones)}, {cmd:e(nodes)}, {cmd:e(runs)}, {cmd:e(runs_used)}, {cmd:e(iterations)}, {cmd:e(thin)}, {cmd:e(burnin)},
{cmd:e(proposals)}, {cmd:e(acceptance)}, {cmd:e(seed)}; macros {cmd:e(cmd)}, {cmd:e(version)}, {cmd:e(depvar)}, {cmd:e(network)}, {cmd:e(effects)},
{cmd:e(directed)}, {cmd:e(method)}, {cmd:e(title)}.


{marker ref}{...}
{title:References}

{p 4 8 2}Borisenko, A., Byshkin, M., Lomi, A. 2019. A simple algorithm for scalable Monte Carlo inference. arXiv:1901.00533.{p_end}
{p 4 8 2}Daraganova, G., Robins, G. 2013. Autologistic actor attribute models. In Lusher, Koskinen, Robins (eds.), Exponential Random Graph Models for Social Networks, ch. 9. Cambridge University Press.{p_end}
{p 4 8 2}Stivala, A., Wang, P., Lomi, A. 2024. ALAAMEE: Open-source software for fitting autologistic actor attribute models. PLOS Complex Systems 1(4): e0000021.{p_end}
{p 4 8 2}Stivala, A. 2023. Overcoming near-degeneracy in the autologistic actor attribute model. arXiv:2309.07338.{p_end}

{hline}
{pstd}{it:nwalaam version 0.2.0, 18 November 2026}{p_end}{* nwversion}
