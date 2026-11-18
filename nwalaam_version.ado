*! version 0.1.0  18nov2026 nwalaam
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
capture program drop nwalaam_version
program nwalaam_version, rclass
	version 14
	_nwalaam_libcheck
	local ver "0.1.0"
	local date "18nov2026"
	mata: st_local("__p", AlaamPluginPath())
	mata: st_local("__ok", strofreal(AlaamNativeAvailable()))
	capture which nwversion
	di as text "nwalaam version `ver', `date'"
	if "`__ok'" == "1" di as text "  native plugin           : `__p' (ships with the package)"
	else di as text "  native plugin           : missing"
	di as text "  reference               : ALAAMEE (Stivala, Wang and Lomi) and MPNet"
	di as text "  license                 : GNU GPL-3"
	return local version "`ver'"
	return local date "`date'"
end
