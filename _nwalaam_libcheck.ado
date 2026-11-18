*! version 1.1.1  18nov2026 nwalaam
* Copyright (C) 2026 Thomas U. Grund. Part of nwalaam; GNU GPL-3 (see LICENSE).
capture program drop _nwalaam_libcheck
program _nwalaam_libcheck
	version 14
	capture mata: mata which AlaamFit()
	if _rc {
		capture mata: mata mlib index
		capture mata: mata which AlaamFit()
		if _rc {
			di as error "{p 0 0 2}nwalaam: the Mata library lnwalaam.mlib was not found. Reinstall the package with {bf:net install nwalaam, replace}.{p_end}"
			exit 3499
		}
	}
end
