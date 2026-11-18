* The plugin's platform follows the operating system, as the core's NwNativePlatform does: StataNow's console binary on a Mac reports
* c(os) = "Unix" with c(machine_type) starting "Mac" and must load the macOS plugin, never the Linux file.
local nwcore : env NWCORE
if `"`nwcore'"' == "" {
	di as error "set the environment variable NWCORE to the directory of the nwcommands (core) package"
	exit 198
}
mata: mata mlib index
mata: st_local("plat", AlaamPlatform()); st_local("fname", AlaamPluginFilename()); st_local("sub", AlaamPluginSubdir())
di "c(os) = `c(os)', c(machine_type) = `c(machine_type)' -> platform `plat', file `fname', directory `sub'"
if "`c(os)'" == "Windows" local want windows
else if "`c(os)'" == "MacOSX" | substr("`c(machine_type)'", 1, 3) == "Mac" local want macos
else local want unix
assert "`plat'" == "`want'"
assert "`fname'" == "alaam_sim_`want'.plugin"
assert "`sub'" == "`want'"
di "PASS: platform `plat' (c(os) `c(os)')"
