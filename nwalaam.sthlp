{smcl}
{* *! version 1.4.0  21nov2026 author: Thomas Grund}{...}
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
{cmd:nwalaam simulate} [{it:outcome} [{it:netname}]]{cmd:,} {c -(}{opt gen:erate(stub)} | {opt fr:ame(name)}{c )-} [{opt eff:ects(list)} {opt theta(matname)} {opt nsim(#)} {opt burn:in(#)} {opt inter:val(#)} {opt seed(#)} {opt replace} {opt mode(#)} {opt zone(varname)} {opt sampl:er(zoo)}]

{p 8 16 2}
{cmd:estat gof} [{cmd:,} {opt nsim(#)} {opt extra(list)} {opt noex:tra} {opt burn:in(#)} {opt inter:val(#)} {opt seed(#)} {opt plot} {opt distplot} {opt all} {opt name(name)} {opt sav:ing(filename)} {opt exp:ort(filename [, options])} {it:graph_options}]

{p 8 16 2}
{cmd:estat mcmcdiag} [{cmd:,} {opt plot} {opt name(name)} {opt sav:ing(filename)} {opt exp:ort(filename [, options])} {it:graph_options} {opt thresh:old(#)}]

{p 8 16 2}
{cmd:nwalaam version}

{pstd}
{it:outcome} is a numeric variable of the current dataset with the values 0 and 1 (missing values are allowed: these actors are
never changed by the sampler and count as 0 in every statistic), one row per actor. When the dataset carries the node labels
(variable {cmd:_nwnode}, as after {helpb nwload} with {cmd:xvars}) the rows are matched to the actors by label; otherwise the
dataset must have one row per actor, in actor order. {it:netname} is the network (default: the current network); it may be
undirected or directed, one-mode, or a two-mode network (see {help nwalaam##twomode:Two-mode networks}). Valued networks are used as binary (a tie where the value is not 0).
Self-loops of a one-mode network are ignored by the statistics; the command prints a one-line note when the network has some.
Covariates named in an effect such as {cmd:ego(}{it:var}{cmd:)} must not have missing values for any actor (of the level of the outcome, for a two-mode network): otherwise
the command stops with error 416 and names the variable and the number of missing values. Impute or recode them first (for a centred covariate, the mean is 0).

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt eff:ects(list)}}the effects of the model (required); see below{p_end}
{synopt :{opt meth:od(sa|ee)}}estimator: {cmd:sa}, stochastic approximation (maximum likelihood; the default), or {cmd:ee}, the equilibrium-expectation update of ALAAMEE{p_end}
{synopt :{opt mple}}maximum pseudo-likelihood estimate only (a logit of the outcome on the change statistics); no standard errors from sampling{p_end}
{synopt :{opt sait:er(#)}}length of the stochastic approximation: iterations of each phase; default {cmd:500}{p_end}
{synopt :{opt prop:osals(#)}}proposals per iteration; default {cmd:1000}{p_end}
{synopt :{opt start(matname)}}1 x {it:P} starting values (default: Algorithm S, 100 iterations){p_end}
{synopt :{opt mode(#)}}two-mode networks: the level (1 or 2) on which the outcome is defined; default {cmd:1}{p_end}
{synopt :{opt zone(varname)}}snowball sample: the zone (0, 1, ..., outermost) of every actor; the conditional sampler holds the outermost zone fixed{p_end}
{synopt :{opt sampl:er(zoo)}}zero-or-one sampler for rare (or very common) outcomes; default: the toggle sampler{p_end}
{synopt :{opt runs(#)}}{cmd:method(ee)}: independent runs; default {cmd:30}{p_end}
{synopt :{opt thr:eads(#)}}{cmd:method(ee)}: number of threads that execute the independent runs; default: one per core of the computer, at most one per run; the result does not depend on it{p_end}
{synopt :{opt iter:ations(#)}}{cmd:method(ee)}: iterations of every run; default {cmd:50000}{p_end}
{synopt :{opt thin(#)}}{cmd:method(ee)}: every {it:#}th iteration is stored; default {cmd:100}{p_end}
{synopt :{opt burn:in(#)}}{cmd:method(ee)}: iterations discarded at the start of every run; default {cmd:1000}{p_end}
{synopt :{opt lrate(#)}}{cmd:method(ee)}: learning rate of the update; default {cmd:0.01}{p_end}
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
{bf:Estimators.} The default, {cmd:method(sa)}, is stochastic approximation in the manner of Snijders (the algorithm of MPNet): starting from the
Algorithm S values, a Robbins-Monro procedure moves the coefficients so that the expected statistics of the model equal the observed
ones, which is the maximum likelihood condition, using the covariance of the simulated statistics as the derivative matrix and
averaging the coefficients over the last subphase. The procedure is repeated from its result, up to 20 times, until the
convergence t-ratios (mean over standard deviation of the difference between simulated and observed statistics in a final
phase at the estimate) are all below 0.1 in absolute value. The covariance matrix of the estimate is the inverse of the
covariance of the simulated statistics at the estimate (the information); on small networks, where the exact maximum likelihood estimate
can be computed, the standard errors agree with the exact ones within a few percent. {cmd:saiter()} lengthens the phases.

{pstd}
{cmd:method(ee)} is the equilibrium-expectation (EE) algorithm of Borisenko, Byshkin and Lomi in the form used by ALAAMEE: every
iteration proposes {cmd:proposals()} toggles of the outcome of a random node, accepts with the Metropolis probability, and moves
every parameter by a relative step {cmd:lrate()} against the sign of the accumulated difference between the simulated and the observed
statistics. {cmd:runs()} independent runs are made; the estimate is their mean. The covariance matrix is the inverse of the
batch-means covariance of the accumulated differences (the information part, averaged over the runs) plus the Monte-Carlo
variance of the mean of the runs. {cmd:e(tratio)} holds the convergence t-ratios; a ratio above 0.3 in absolute value is flagged, as in ALAAMEE.
{cmd:threads()} sets how many runs are computed at the same time (default: one thread per core, never more threads than runs, so with the default
{cmd:runs(30)} a computer with 20 cores works on 20 runs at once and the fit takes about the time of two runs). Every run has its own random-number stream, derived from {cmd:seed()}
and the number of the run, and its own memory, so the estimates, the standard errors and the diagnostics are exactly the same for {cmd:threads(1)}, {cmd:threads(4)}
and the default. Stochastic approximation, {cmd:mple}, {cmd:simulate} and {cmd:estat gof} follow one chain each and use one thread.

{pstd}
{bf:Which estimator.} Stochastic approximation gives standard errors that agree with the exact ones. In checks on small networks, where the exact maximum likelihood estimate can be computed by enumeration, the equilibrium-expectation estimates were up to 0.5 standard errors away from the exact estimates and its standard errors were between 0.7 and 1.2 times the exact ones; on larger networks they can be too small (check them against a parametric bootstrap before relying on them). The two estimators also differ in the estimates: where the
exact maximum likelihood estimate can be computed by enumerating all outcome vectors, stochastic approximation reproduces it,
while the equilibrium-expectation update stops at slightly different values. On the karate club (34 actors; Density, Activity and Contagion
with the faction as outcome) the exact estimates are 0.535, -0.832 and 1.461; stochastic approximation gives 0.52, -0.83, 1.46 and the
equilibrium-expectation update gives -0.03, -0.75 and 1.53 (the same as ALAAMEE's own program), a Density difference of 0.6 standard
errors; on a random network of 16 actors the equilibrium-expectation estimates lie 0.15 to 0.24 standard errors from the exact ones. For
that reason stochastic approximation is the default; {cmd:method(ee)} reproduces ALAAMEE and is the choice when results are to be compared with it.
The pseudo-likelihood estimate ({cmd:mple}) treats the log odds of each actor as independent and its standard errors are too small; use it
for starting values and for a quick look.

{pstd}
{bf:Starting values and divergence.} The default start is the one of Snijders' Algorithm S (a short simulation, 100 iterations) whose step for every coefficient is limited to the range -0.1 to 0.1, so
that the start stays moderate. (The ALAAMEE program limits only the upward step; the symmetric limit is a deliberate difference of this package.) Should the stochastic approximation
still diverge, the command restarts from the same start limited to the range -10 to 10, and
if needed from zero; the start used is in {cmd:e(b_start)}. If it still diverges it stops with error 430; try another {opt seed()}, fewer effects or a smaller {opt lrate()}. {cmd:start(}{it:matname}{cmd:)} gives the start by hand.

{pstd}
All computing is done in a compiled plugin; there is no slower fallback (error 601: plugin missing; 3495: plugin too old).
The change statistics of the plugin were compared with ALAAMEE's own for random undirected, directed and two-mode networks with missing
values (differences below 1e-15).


{pstd}
{bf:Interpretation.} A cross-sectional ALAAM describes how an outcome is distributed over a given network: a positive contagion
coefficient says that actors with the outcome are tied to each other more often than the other effects explain. It does not
identify social influence: the same pattern arises when actors choose friends who are like them (selection) or share a
setting (confounding). The network is taken as given; use a longitudinal model such as {helpb nwsaom} to separate influence and selection.
The contagion coefficient alone cannot tell whether one exposed neighbour or several are needed: {cmd:thresholdcontagion}, {cmd:fractioncontagion} and
{cmd:sharecontagion} let a specification represent exposure thresholds (complex contagion), but a significant coefficient of such a term is still a
statement about how the outcome is distributed over the network, not about the process that produced it.
The exposure effects {cmd:thresholdexposure}, {cmd:fractionexposure} and {cmd:shareexposure} apply the same statistics to the nodes {it:without}
the outcome: they count nodes that are exposed but have not adopted. A negative coefficient means that there are fewer such exposed nodes than the
other effects imply, that is, exposed nodes tend to have the outcome. Together with the contagion effects they let a specification
represent exposure thresholds more fully; as with all effects, they describe the distribution of the outcome, not the process.


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
{p2col :{cmd:thresholdcontagion}[{cmd:(}{it:k}{cmd:)}]}number of nodes with the outcome that have at least {it:k} neighbours with the outcome; {it:k} a positive integer, default 2{p_end}
{p2col :{cmd:fractioncontagion}[{cmd:(}{it:q}{cmd:)}]}number of nodes with the outcome of which at least the share {it:q} of the neighbours have the outcome (nodes without neighbours do not count); 0 < {it:q} <= 1, default 0.5{p_end}
{p2col :{cmd:sharecontagion}}sum over nodes with the outcome of the share of their neighbours that have the outcome (0 for a node without neighbours){p_end}
{p2col :{cmd:thresholdexposure}[{cmd:(}{it:k}{cmd:)}]}number of nodes {it:without} the outcome that have at least {it:k} neighbours with the outcome; {it:k} a positive integer, default 2{p_end}
{p2col :{cmd:fractionexposure}[{cmd:(}{it:q}{cmd:)}]}number of nodes without the outcome of which at least the share {it:q} of the neighbours have the outcome (nodes without neighbours do not count); 0 < {it:q} <= 1, default 0.5{p_end}
{p2col :{cmd:shareexposure}}sum over nodes without the outcome of the share of their neighbours that have the outcome (0 for a node without neighbours){p_end}
{p2col :{cmd:ego(}{it:var}{cmd:)}}sum of the numeric covariate over the nodes with the outcome (centre the covariate){p_end}
{p2col :{cmd:alter(}{it:var}{cmd:)}}sum over the nodes with the outcome of the covariate of their neighbours{p_end}
{p2col :{cmd:same(}{it:var}{cmd:)}, {cmd:diff(}{it:var}{cmd:)}}neighbours with the same / a different category{p_end}
{p2col :{cmd:samecontagion(}{it:var}{cmd:)}}ties between two nodes with the outcome and the same category{p_end}
{p2col :{cmd:samepartneractivity(}{it:var}{cmd:)}, {cmd:diffpartneractivity(}{it:var}{cmd:)}}two-paths from a node with the outcome to a node of the same / a different category{p_end}
{p2col :{cmd:sameindirect(}{it:var}{cmd:)}, {cmd:diffindirect(}{it:var}{cmd:)}}two-paths between two nodes with the outcome of the same / different categories{p_end}
{p2col :{cmd:alterbinarytwostar1(}{it:var}{cmd:)}}sum over the nodes with the outcome of (degree - 1) of the neighbours with the binary attribute 1{p_end}
{p2col :{cmd:alterbinarytwostar2(}{it:var}{cmd:)}}two-paths between two nodes with the outcome whose middle node has the binary attribute 1{p_end}
{p2col :{cmd:settinghomophily(}{it:netname}{cmd:)}}ties in a second (setting) network between two nodes with the outcome; the main network is not required to have the tie{p_end}
{p2col :{cmd:geographichomophily(}{it:matname}{cmd:)}}sum over all other nodes with the outcome of their distance to the node (a distance matrix, see below){p_end}
{p2col :{cmd:contagiondist(}{it:matname}{cmd:)}}sum over the neighbours with the outcome of their distance to the node{p_end}
{p2colreset}{...}

{pstd}{bf:Directed networks}{p_end}
{p2colset 9 33 35 2}{...}
{p2col :{cmd:density}}number of nodes with the outcome{p_end}
{p2col :{cmd:sender}, {cmd:receiver}}out-degree, in-degree of the nodes with the outcome{p_end}
{p2col :{cmd:reciprocity}}mutual ties of the nodes with the outcome{p_end}
{p2col :{cmd:outstar2}, {cmd:outstar3}, {cmd:instar2}, {cmd:instar3}}out- and in-stars at the nodes with the outcome{p_end}
{p2col :{cmd:mixedstar}, {cmd:mixedsource}, {cmd:mixedsink}}mixed two-stars (in and out), and with the other end source or sink{p_end}
{p2col :{cmd:contagion}}ties to nodes with the outcome, in and out; a mutual tie counts twice{p_end}
{p2col :{cmd:contagionrecip}}mutual ties between two nodes with the outcome{p_end}
{p2col :{cmd:gwcontagion}[{cmd:(}{it:a}{cmd:)}], {cmd:logcontagion}, {cmd:powercontagion}[{cmd:(}{it:b}{cmd:)}]}as in undirected networks, with the number of ties (in and out) to nodes with the outcome{p_end}
{p2col :{cmd:thresholdcontagion}[{cmd:(}{it:k}{cmd:)}], {cmd:fractioncontagion}[{cmd:(}{it:q}{cmd:)}], {cmd:sharecontagion}, {cmd:thresholdexposure}[{cmd:(}{it:k}{cmd:)}], {cmd:fractionexposure}[{cmd:(}{it:q}{cmd:)}], {cmd:shareexposure}}as in undirected networks; {it:n} is the number of ties (in and out; a mutual tie counts twice) to nodes with the outcome and the share is {it:n} divided by the number of all ties (in and out) of the node{p_end}
{p2col :{cmd:transitive1}, {cmd:transitive3}}transitive triples (a->b, b->c, a->c): the number of triples in which the node with the outcome is the middle node b; all three nodes with the outcome{p_end}
{p2col :{cmd:transitived1}, {cmd:transitiveu1}}transitive triples in which the node with the outcome is the source a, the sink c{p_end}
{p2col :{cmd:cyclic1}, {cmd:cyclic3}}cyclic triples (a->b->c->a) with the node with the outcome in it; all three nodes with the outcome{p_end}
{p2col :{cmd:alterinstar2}, {cmd:alteroutstar2}}two-stars between a node with the outcome and another node with the outcome through a common target / source{p_end}
{p2col :{cmd:gwsender}[{cmd:(}{it:a}{cmd:)}], {cmd:gwreceiver}[{cmd:(}{it:a}{cmd:)}]}geometrically weighted out- and in-degree{p_end}
{p2col :{cmd:ego(}{it:var}{cmd:)}, {cmd:alter(}{it:var}{cmd:)}}covariate of the node itself, of its neighbours (in and out){p_end}
{p2col :{cmd:same(}{it:var}{cmd:)}, {cmd:diff(}{it:var}{cmd:)}, {cmd:samecontagion(}{it:var}{cmd:)}}ties (in and out) to nodes of the same / a different category; with the outcome and the same category{p_end}
{p2col :{cmd:sendermatch(}{it:var}{cmd:)}, {cmd:sendermismatch(}{it:var}{cmd:)}}out-ties to nodes of the same / a different category{p_end}
{p2col :{cmd:receivermatch(}{it:var}{cmd:)}, {cmd:receivermismatch(}{it:var}{cmd:)}}in-ties from nodes of the same / a different category{p_end}
{p2col :{cmd:reciprocitymatch(}{it:var}{cmd:)}, {cmd:reciprocitymismatch(}{it:var}{cmd:)}}mutual ties with nodes of the same / a different category{p_end}
{p2col :{cmd:settinghomophily(}{it:netname}{cmd:)}, {cmd:geographichomophily(}{it:matname}{cmd:)}}as in undirected networks{p_end}
{p2colreset}{...}
The definitions are those of ALAAMEE (Stivala, Wang and Lomi), including its change statistics of the directed gw, log and power contagion, which sum over the
out-neighbours with the outcome (counting their in-neighbours with the outcome) and the in-neighbours with the outcome (counting their out-neighbours with the
outcome). The threshold, fraction and share contagion effects are not part of ALAAMEE; they are exact counts of the nodes with the outcome by their exposure (the
number or the share of their neighbours with the outcome). They are not defined for two-mode networks, whose ties join nodes of different levels. Continuous covariates are used as they are: centre them.

{pstd}{bf:Setting network and distance matrix.} {it:netname} of {cmd:settinghomophily()} is another network of the same actors (one-mode; direction is ignored). {it:matname}
is a Stata matrix (or, for large networks, a Mata matrix) with one row and one column for each actor of the network, in the order of the network; a
distance matrix is not changed or centred. {cmd:geographichomophily} visits every actor for every proposal and is slow on large networks.

{marker twomode}{...}
{pstd}{bf:Two-mode networks.} The outcome is defined on the actors of one level, chosen with {opt mode(#)}; the other level is the network of
events, groups or affiliations the actors are tied to. The outcome and the covariates are variables of the current dataset: either one row for each actor of
the outcome level (in the order of the network), or, with node labels in the dataset, matched by label, or one row for every node of both levels. Covariates of
the other level (for {cmd:alter()} and the binary two-stars) need rows for all nodes. The effects are {cmd:density}, {cmd:activity}, {cmd:twostar}, {cmd:threestar},
{cmd:gwactivity}, {cmd:alter2star1}, {cmd:alter2star2}, {cmd:fourcycle1} (four-cycles through the actor and another actor of the same level, whatever its outcome),
{cmd:fourcycle2} (the other actor also has the outcome), {cmd:ego()}, {cmd:alter()}, {cmd:samepartneractivity()}, {cmd:diffpartneractivity()}, {cmd:sameindirect()},
{cmd:diffindirect()}, {cmd:alterbinarytwostar1()} and {cmd:alterbinarytwostar2()}, {cmd:settinghomophily()} and {cmd:geographichomophily()}; contagion effects do
not exist because there are no ties within a level.

{pstd}{bf:Samplers.} The default toggle sampler proposes to change the outcome of a random actor. With rare (or very common) outcomes most proposals are wasted;
{cmd:sampler(zoo)} proposes a move to 1 or a move to 0 with equal probability and then a random actor with the other value, correcting the acceptance probability
for the numbers of actors with each value. (The correction is derived exactly, including the cases with all or no actor having the outcome, and checked against the exact
maximum likelihood on a small network.) {cmd:zone(varname)} declares snowball-sampled data: the zones 0 (seeds), 1, ..., and the actors of the outermost zone keep their observed outcome while
the others change (the conditional sampler of Stivala et al. 2020); not with a two-mode network or {cmd:sampler(zoo)}.
Geometrically weighted statistics are the remedy for near-degeneracy of the activity and star statistics on larger networks
(Stivala 2023); read their sign with care.


{marker display}{...}
{title:Display}

{pstd}The coefficient names are shown in full. At a line size below 100 the confidence intervals are left out so that names of up to 32
characters fit; {opt ci} (or an explicit {opt level(#)}) shows them (long names are then cut), {opt noci} leaves them out. {opt ci} and {opt noci}
together are an error. Other display options of Stata's estimation commands ({opt nopvalues}, {opt cformat()}, ...) are passed on.


{marker simulate}{...}
{title:nwalaam simulate}

{pstd}Draws outcome vectors from the model by the toggle sampler (the same sampler as the estimation) and stores them: after a fit
without further arguments, with the fitted coefficients; otherwise with {it:outcome}, {it:netname}, {opt effects()} and the coefficients
{opt theta()} (a 1 x {it:P} matrix). {opt generate(stub)} writes the draws to the variables {it:stub}1, {it:stub}2, ... (one row per actor, like the
outcome; {opt replace} overwrites), {opt frame(name)} to a frame with the variable {cmd:actor} and one variable {cmd:sim}{it:j} per draw. The
chain starts at the observed outcome; {opt burnin()} proposals (default the larger of 10,000 and 100 times the number of
actors) precede the first draw and {opt interval()} proposals (default the larger of 1,000 and 10 times the number of actors) lie
between draws. Missing outcomes stay missing. The same {opt seed()} gives the same draws. Returned in {cmd:r()}: {cmd:r(nsim)}, {cmd:r(seed)},
{cmd:r(burnin)}, {cmd:r(interval)} and the observed statistics {cmd:r(observed)}.


{marker gof}{...}
{title:estat gof}

{pstd}Simulates {opt nsim()} outcome vectors (default 1000) from the fitted model and compares the statistics of the simulated vectors with
the observed ones: for each statistic the t-ratio (simulated mean - observed) / simulated standard deviation, and the Mahalanobis
distance of the observed statistics from the simulated ones (pseudo-inverse of the covariance matrix if it is singular). The statistics of the model
come first (their t-ratios should be near 0 after a converged fit), then statistics that are not in the model ({opt extra(list)}: any effects;
default, as in ALAAMEE, two-stars, three-stars, partner-activity two-paths, triangles, indirect-partner and partner-partner effects for
undirected networks and the sender, receiver, reciprocity, star, mixed star, contagion and triangle effects for directed ones;
{opt noextra} none). A t-ratio beyond 0.3 in absolute value is marked. {opt plot} draws the t-ratios; {opt distplot} draws, for every statistic of the model ({opt all}: also the other statistics), the histogram of the simulated values with the observed value as a red line
(not part of ALAAMEE). {opt name()} names the graph (default {cmd:nwalaam_gof} and {cmd:nwalaam_gofdist}), {opt saving()} saves it as a {cmd:.gph} file, {opt export()} exports it with {helpb graph export}
(the extension of the file name sets the format), and further options go to {helpb twoway} ({opt plot}) or {helpb graph combine} ({opt distplot}); the plots need nwcommands (core) 3.29.0 or later. The titles of {opt distplot} are black; {opt title()} replaces the default title. {cmd:r(gof)} holds the table
(observed, simulated mean, standard deviation, t-ratio), {cmd:r(mahalanobis)} the distance. The simulation can be misleading
when the distribution of the outcome has two modes (all or none of the actors with the outcome, a common near-degeneracy of
contagion models on small networks): the simulated mean statistics then lie between the modes and no t-ratio is small.


{marker mcmcdiag}{...}
{title:estat mcmcdiag}

{pstd}Convergence diagnostics of the fit. After {cmd:method(sa)}: the estimate, the convergence t-ratio of each coefficient in the final phase
(below 0.1 in absolute value is converged) and the autocorrelation of the simulated differences; after {cmd:method(ee)}: per
coefficient the mean and the standard deviation over the runs, the mean and maximum convergence t-ratio, the number of runs beyond the
threshold ({opt threshold()}, default 0.3) and the lag-1 autocorrelation of the thinned values of the coefficient and of the differences. {opt plot} draws the
trace of the first EE run (coefficients and differences against the iteration) or the differences of the final phase of the
stochastic approximation (default graph name {cmd:nwalaam_trace}; {opt saving()}, {opt export()} and graph options as for {cmd:estat gof}, passed to {helpb graph combine}). Not after {cmd:mple} (r(498)).


{marker results}{...}
{title:Stored results}

{pstd}{cmd:nwalaam} stores in {cmd:e()}: {cmd:e(b)}, {cmd:e(V)}; {cmd:e(obs_stats)} (observed statistics), {cmd:e(b_mple)} (pseudo-likelihood
estimate; missing when it does not exist, for example under perfect separation), {cmd:e(b_start)}, {cmd:e(tratio)}, {cmd:e(diag)};
scalars {cmd:e(N)}, {cmd:e(n_ones)}, {cmd:e(nodes)}, {cmd:e(runs)}, {cmd:e(runs_used)}, {cmd:e(iterations)}, {cmd:e(thin)}, {cmd:e(burnin)},
{cmd:e(proposals)}, {cmd:e(acceptance)}, {cmd:e(seed)} (EE), {cmd:e(repetitions)}, {cmd:e(sa_converged)}, {cmd:e(saiter)} (SA); matrices {cmd:e(run_diag)}, {cmd:e(trace1)}; macros {cmd:e(cmd)}, {cmd:e(version)}, {cmd:e(depvar)}, {cmd:e(network)}, {cmd:e(effects)},
{cmd:e(directed)}, {cmd:e(twomode)}, {cmd:e(modeopt)}, {cmd:e(zone)}, {cmd:e(sampler)}, {cmd:e(method)} ({cmd:sa}, {cmd:ee} or {cmd:mple}), {cmd:e(effects_spec)}, {cmd:e(title)}.


{marker ref}{...}
{title:Manual}

{pstd}
A user manual (PDF in the {cmd:manual} folder of the package source) is archived at Zenodo, {browse "https://doi.org/10.5281/zenodo.24397901":doi:10.5281/zenodo.24397901}, all versions. Please cite it by that address and report the package version printed by {cmd:nwalaam version}.
{p_end}


{title:References}

{p 4 8 2}Borisenko, A., Byshkin, M., Lomi, A. 2019. A simple algorithm for scalable Monte Carlo inference. arXiv:1901.00533.{p_end}
{p 4 8 2}Daraganova, G., Robins, G. 2013. Autologistic actor attribute models. In Lusher, Koskinen, Robins (eds.), Exponential Random Graph Models for Social Networks, ch. 9. Cambridge University Press.{p_end}
{p 4 8 2}Snijders, T. A. B. 2002. Markov chain Monte Carlo estimation of exponential random graph models. Journal of Social Structure 3(2).{p_end}
{p 4 8 2}Stivala, A., Wang, P., Lomi, A. 2024. ALAAMEE: Open-source software for fitting autologistic actor attribute models. PLOS Complex Systems 1(4): e0000021.{p_end}
{p 4 8 2}Stivala, A., Gallagher, H. C., Rolls, D. A., Wang, P., Robins, G. L. 2020. Using sampled network data with the autologistic actor attribute model. arXiv:2002.00849.{p_end}
{p 4 8 2}Stivala, A. 2023. Overcoming near-degeneracy in the autologistic actor attribute model. arXiv:2309.07338.{p_end}

{hline}
{pstd}{it:nwalaam version 1.4.0, 21 November 2026}{p_end}{* nwversion}
