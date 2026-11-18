# NOTICE

<!--version-->**Version 1.0.0** (released 18 November 2026)<!--/version-->

nwalaam is Copyright (C) 2026 Thomas U. Grund and is distributed under the GNU General Public License, version 3 (file `LICENSE`).

## What this package follows

`nwalaam` fits autologistic actor attribute models (ALAAMs; Robins, Pattison and Elliott 2001; Daraganova and Robins 2013) with the
equilibrium-expectation (EE) estimator of Borisenko, Byshkin and Lomi (2019), as implemented in the software ALAAMEE.
* **ALAAMEE** (Python and R; Stivala, Wang and Lomi). License: MIT. https://github.com/stivalaa/ALAAMEE
  Reference: Stivala, A., Wang, P., and Lomi, A. (2024). ALAAMEE: Open-source software for fitting autologistic actor attribute
  models. PLOS Complex Systems 1(4): e0000021.
* **MPNet / PNet** (Wang, Robins and Pattison) is the other reference implementation of the model.
The code of nwalaam is written independently from the published method; no code of ALAAMEE is included. The change statistics follow
the documented definitions of ALAAMEE. ALAAMEE's MIT licence requires its copyright notice and permission notice to be kept where its
code is copied; none is copied.

## Validation

Every change statistic is compared with its definition (a count over the whole outcome vector), the pseudo-likelihood estimate with
Stata's `logit`, and the EE estimates and standard errors with the exact maximum likelihood estimate and Fisher information of a small
network (all 2^16 outcome vectors enumerated); see `cscripts/test_nwalaam_*.do`. The karate club example of ALAAMEE (`dev/karate`,
from the ALAAMEE repository, MIT licence) is reproduced.

## AI assistance

The software was developed with Claude Code (Anthropic), an AI coding assistant, working under the author's direction. The author
specified the scope, reviewed the code, and is responsible for the package. Its correctness rests on validation, not on the assistant.

## How to cite

Grund, T. U. (2026). nwalaam: a Stata package (version 1.0.0).
Please also cite the core package (Grund, T. U., nwcommands) and the ALAAMEE paper above.
