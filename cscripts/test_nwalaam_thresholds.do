cscript
* Threshold, fraction and share contagion (thresholdcontagion(k), fractioncontagion(q), sharecontagion) and the exposure effects of the
* actors without the outcome (thresholdexposure(k), fractionexposure(q), shareexposure): the observed statistics against
* a direct count, the names and defaults, the refused arguments, fits on an undirected and a directed network, and the independence of the
* result of the number of threads. (The change statistics are checked against the definition in test_nwalaam_changestats.)
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
_nwsyntax karate

* ---- the observed statistics of the karate club, counted directly
mata:
A = *(`netobj'->get_matrix())
Y = st_data(., "faction")
n = rows(A)
A0 = editmissing(A, 0)
nj = A0 * Y
dj = rowsum(A0)
zthr2 = sum(Y :* (nj :>= 2))
zthr3 = sum(Y :* (nj :>= 3))
zfr = sum(Y :* (dj :> 0) :* (nj :>= (0.3 :* dj :- 1e-9)))
zsh = sum(Y :* (dj :> 0) :* nj :/ (dj + (dj :== 0)))
end
mata: AlaamChangeStats(A, 0, Y, (43, 43, 44, 45), (2, 3, 0.3, 0), J(1, 4, 0), J(0, 34, .), D = ., z = .)
mata: st_matrix("OBS", z)
mata: st_matrix("REF", (zthr2, zthr3, zfr, zsh))
matrix list OBS
matrix list REF
assert mreldif(OBS, REF) < 1e-12
assert OBS[1, 1] > 0 & OBS[1, 3] > 0 & OBS[1, 4] > 0
di "PASS: the observed statistics equal a direct count"
* the exposure effects count the actors WITHOUT the outcome
mata: zte2 = sum((1 :- Y) :* (nj :>= 2)); zte3 = sum((1 :- Y) :* (nj :>= 3))
mata: zfe = sum((1 :- Y) :* (dj :> 0) :* (nj :>= (0.3 :* dj :- 1e-9)))
mata: zse = sum((1 :- Y) :* (dj :> 0) :* nj :/ (dj + (dj :== 0)))
mata: AlaamChangeStats(A, 0, Y, (46, 46, 47, 48), (2, 3, 0.3, 0), J(1, 4, 0), J(0, 34, .), D = ., z = .)
mata: st_matrix("OBS", z)
mata: st_matrix("REF", (zte2, zte3, zfe, zse))
assert mreldif(OBS, REF) < 1e-12
assert OBS[1, 1] > 0 & OBS[1, 3] > 0
di "PASS: the observed exposure statistics equal a direct count (actors without the outcome)"

* ---- names and defaults
local fast "seed(11)"
mata: alaam_effects("density thresholdcontagion fractioncontagion sharecontagion", 0, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(invtokens(nm) == "density thresholdcontagion fractioncontagion sharecontagion")
mata: assert(et == (1, 43, 44, 45)); assert(ep == (0, 2, 0.5, 0))
mata: alaam_effects("density thresholdcontagion(3) fractioncontagion(0.3) sharecontagion", 0, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(invtokens(nm) == "density thresholdcontagion_3 fractioncontagion_0p3 sharecontagion"); assert(ep == (0, 3, 0.3, 0))
mata: alaam_effects("thresholdcontagion fractioncontagion(0.25) sharecontagion", 1, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(et == (137, 138, 139)); assert(invtokens(nm) == "thresholdcontagion fractioncontagion_0p25 sharecontagion")
mata: alaam_effects("density thresholdexposure fractionexposure shareexposure", 0, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(invtokens(nm) == "density thresholdexposure fractionexposure shareexposure"); assert(et == (1, 46, 47, 48)); assert(ep == (0, 2, 0.5, 0))
mata: alaam_effects("thresholdexposure(3) fractionexposure(0.3) shareexposure", 1, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(et == (140, 141, 142)); assert(invtokens(nm) == "thresholdexposure_3 fractionexposure_0p3 shareexposure")
di "PASS: names, codes and defaults (k = 2, q = 0.5), undirected and directed"
* the default parameters are the ones of the statistics
mata: AlaamChangeStats(A, 0, Y, (43, 44), (2, 0.5), J(1, 2, 0), J(0, 34, .), D = ., z = .)
mata: st_matrix("OBS", z)
mata: st_matrix("REF", (zthr2, sum(Y :* (dj :> 0) :* (nj :>= (0.5 :* dj :- 1e-9)))))
assert mreldif(OBS, REF) < 1e-12

* ---- refused arguments
foreach bad in "thresholdcontagion(0)" "thresholdcontagion(1.5)" "thresholdcontagion(-1)" "fractioncontagion(0)" "fractioncontagion(1.5)" "fractioncontagion(-0.2)" "sharecontagion(2)" "thresholdexposure(0)" "thresholdexposure(2.5)" "fractionexposure(0)" "fractionexposure(1.2)" "shareexposure(1)" {
	capture nwalaam faction karate, effects(density `bad') `fast'
	assert _rc == 198
}
// fraction exactly 1 is allowed
mata: alaam_effects("fractioncontagion(1)", 0, 0, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
mata: assert(ep == 1)
// not for two-mode networks
capture mata: alaam_effects("thresholdcontagion", 0, 1, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
assert _rc == 198
capture mata: alaam_effects("fractionexposure", 0, 1, et = J(1, 0, .), ep = J(1, 0, .), ec = J(1, 0, .), nm = J(1, 0, ""), cn = J(1, 0, ""), sn = J(1, 0, ""), dn = J(1, 0, ""))
assert _rc == 198
di "PASS: bad arguments and two-mode networks -> r(198)"

* ---- an undirected fit with all three, and EE with threads
nwalaam faction karate, effects(density contagion thresholdcontagion(2)) `fast'
assert "`: colnames e(b)'" == "density contagion thresholdcontagion_2"
capture nwalaam faction karate, effects(density contagion thresholdexposure(2) shareexposure) `fast'
assert inlist(_rc, 0, 430)
if _rc == 0 assert "`: colnames e(b)'" == "density contagion thresholdexposure_2 shareexposure"
local eff "density contagion thresholdcontagion(2) fractioncontagion(0.5)"
nwalaam faction karate, effects(`eff') method(ee) runs(6) iterations(3000) thin(50) burnin(500) seed(7) threads(1)
matrix B1 = e(b)
matrix V1 = e(V)
foreach th in 2 6 {
	nwalaam faction karate, effects(`eff') method(ee) runs(6) iterations(3000) thin(50) burnin(500) seed(7) threads(`th')
	assert mreldif(e(b), B1) == 0
	assert mreldif(e(V), V1) == 0
}
di "PASS: undirected fits; EE independent of threads"

* ---- gof and simulate with the new effects
nwalaam faction karate, effects(density contagion thresholdcontagion(2)) `fast'
estat gof, nsim(50) seed(1) noextra
di "PASS: estat gof"

* ---- a directed network
set seed 123
clear
set obs 140
gen from = ceil(runiform() * 18)
gen to = ceil(runiform() * 18)
drop if from == to
nwfromedge from to, directed name(dnet)
clear
set obs 18
gen byte yd = runiform() < 0.45
mata: AD = *(`netobj'->get_matrix())
_nwsyntax dnet
mata: AD = *(`netobj'->get_matrix())
mata: YD = st_data(., "yd")
mata: AD0 = editmissing(AD, 0)
mata: nn = AD0 * YD + AD0' * YD
mata: dd = rowsum(AD0) + colsum(AD0)'
mata: AlaamChangeStats(AD, 1, YD, (137, 138, 139), (2, 0.5, 0), J(1, 3, 0), J(0, 18, .), D = ., z = .)
mata: st_matrix("OBS", z)
mata: st_matrix("REF", (sum(YD :* (nn :>= 2)), sum(YD :* (dd :> 0) :* (nn :>= (0.5 :* dd :- 1e-9))), sum(YD :* (dd :> 0) :* nn :/ (dd + (dd :== 0)))))
assert mreldif(OBS, REF) < 1e-12
mata: AlaamChangeStats(AD, 1, YD, (140, 141, 142), (2, 0.5, 0), J(1, 3, 0), J(0, 18, .), D = ., z = .)
mata: st_matrix("OBS", z)
mata: st_matrix("REF", (sum((1 :- YD) :* (nn :>= 2)), sum((1 :- YD) :* (dd :> 0) :* (nn :>= (0.5 :* dd :- 1e-9))), sum((1 :- YD) :* (dd :> 0) :* nn :/ (dd + (dd :== 0)))))
assert mreldif(OBS, REF) < 1e-12
capture nwalaam yd dnet, effects(density thresholdcontagion(2) sharecontagion) saiter(50) seed(3)
assert inlist(_rc, 0, 430)
if _rc == 0 assert "`: colnames e(b)'" == "density thresholdcontagion_2 sharecontagion"
di "PASS: directed network (statistics counted from in and out ties, a mutual tie twice)"
di "test_nwalaam_thresholds: ALL OK"
