cscript
* There is no Mata fallback: without the plugin file nwalaam stops with r(601); with a file that is not a plugin with r(3495); with the
* plugin of this tree it fits. Works in a scratch copy of the ado files (without lib/plugins) and puts a file of the plugin's name
* first on the adopath. The plugin directories are removed from the adopath whatever slash style Stata stores them in.
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
local root `"`c(pwd)'"'
quietly do `"`nwcore'/unw_core.do"'
mata: mata mlib index
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
mata: st_local("pname", AlaamPluginFilename())
mata: st_local("plugsrc", AlaamPluginPath())
confirm file `"`plugsrc'"'
tempfile __scx
local sc `"`__scx'_d"'
capture mkdir `"`sc'"'
capture erase `"`sc'/`pname'"'
foreach f in nwalaam.ado _nwalaam_libcheck.ado nwalaam_version.ado {
	quietly copy `"`root'/`f'"' `"`sc'/`f'"', replace
}
local r2 = subinstr(`"`root'"', "/", "\", .)
local r3 = subinstr(`"`root'"', "\", "/", .)
foreach r in `"`root'"' `"`r2'"' `"`r3'"' {
	foreach sub in "/lib/plugins/macos" "/lib/plugins/unix" "/lib/plugins/windows" "\lib\plugins\macos" "\lib\plugins\unix" "\lib\plugins\windows" "" {
		capture adopath - `"`r'`sub'"'
	}
}
adopath ++ `"`sc'"'
capture program drop nwalaam
capture program drop _nwalaam_libcheck
capture noisily nwalaam faction karate, effects(density contagion) runs(2) iterations(1000)
assert _rc == 601
assert "`c(frame)'" == "default"
di as text "  PASS: no plugin file -> r(601)"
file open fh using `"`sc'/`pname'"', write text replace
file write fh "this is not a plugin" _n
file close fh
capture noisily nwalaam faction karate, effects(density contagion) runs(2) iterations(1000)
assert _rc == 3495
assert "`c(frame)'" == "default"
capture frame __nwalaam_native: describe
assert _rc != 0
di as text "  PASS: a file that cannot be loaded -> r(3495), no frame left behind"
erase `"`sc'/`pname'"'
local sc2 `"`sc'/new"'
capture mkdir `"`sc2'"'
quietly copy `"`plugsrc'"' `"`sc2'/`pname'"', replace
capture shell chmod +x "`sc2'/`pname'"
adopath ++ `"`sc2'"'
nwalaam faction karate, effects(density contagion) runs(2) iterations(4000) thin(50) burnin(500) seed(1)
assert e(runs) == 2
di as text "  PASS: with the plugin the model fits"
di "test_nwalaam_noplugin: ALL OK"
