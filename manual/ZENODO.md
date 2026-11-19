# Zenodo record for the nwalaam manual

Metadata ready to paste into a **new** Zenodo record (the first for this manual; it is not a new version of the nwsaom or nwergm record). Nothing has been uploaded.

**DOI.** None yet. On zenodo.org start a new upload, press *Reserve DOI*, write the reserved all-versions DOI into `\manualdoi` in `nwalaam_manual.tex` (it is empty now, and an empty DOI prints nothing), rebuild with `make -C manual pdf` and upload the rebuilt PDF to the reserved record. The citation text always uses the all-versions DOI (master RELEASING.md B.10).

**Labelling.** The title and the citation of the manual contain no version number. The Zenodo version label is the series "X.Y" (here "1.2"; internal to Zenodo, not part of the citation), and the description says that the manual describes nwalaam X.Y.x. The PDF title page has one line saying which series it describes (`\manualseries`); the exact package version of the examples is `\pkgversion`. The software version is cited separately, as printed by `nwalaam version`.

**Patch releases.** The manual of the series 1.2 holds for 1.2.x; a patch release changes no documented option or number. Title and citation stay unchanged; no new Zenodo version is needed for a patch. A new Zenodo version only for 1.2 or 2.0 with changed documented behaviour (then change `\manualseries`).

## Metadata

**Title.** nwalaam: Autologistic actor attribute models in Stata. Manual

**Version.** 1.2 (the series label, internal to Zenodo: the manual describes nwalaam 1.2.x; not part of the citation)

**Publication date.** The day of the upload.

**Resource type.** Publication, sub-type Software documentation (as the nwsaom and nwergm records).

**Creators.** Grund, Thomas U. (RWTH Aachen University), ORCID 0000-0001-6751-5429 (as on the other records; remove it in the form if the author does not want it).

**License.** Creative Commons Attribution 4.0 International (CC-BY-4.0), as the other manual records. The code stays under GPL-3.0-only.

**Access.** Open.

**Description** (plain text, one paragraph):

> This is the user manual and reference of nwalaam, series 1.2 (examples produced with version 1.2.0 and the nwcommands core 3.27.0), a Stata package for autologistic actor attribute models (ALAAMs): models for a binary outcome on the actors of one fixed network (who smokes, who adopted a practice, who took a side), with effects such as contagion, activity, triangles, covariates, a setting network and a distance matrix. nwalaam builds on the nwcommands core package and follows ALAAMEE (Stivala, Wang and Lomi) in its effects, the definitions of their change statistics, its samplers (toggle, zero-or-one, snowball-conditional) and the equilibrium-expectation estimator, and MPNet in its default estimator, stochastic approximation. The manual describes the installation, the syntax and every option, the effects for undirected, directed and two-mode networks, the three estimators and their accuracy (including a comparison with the exact maximum likelihood estimate of the karate club), the special data features, goodness of fit, convergence diagnostics and simulation from the model, worked examples (the karate club, smoking among Glasgow pupils, the Lazega law firm, the southern women), the validation against ALAAMEE and the exact distribution, and the limits of the package. nwalaam requires the nwcommands core package and is meant for researchers who work in Stata. It is not a replacement for ALAAMEE or MPNet, which remain the reference implementations, and it is not affiliated with or endorsed by their authors. The software is free software under the GNU General Public License, version 3. The package and the manual were prepared with the help of an AI coding assistant working under the author's direction, as stated in the manual. Please cite this manual (all-versions DOI) and report the exact version printed by nwalaam version.

**Keywords.** Stata; nwcommands; nwalaam; autologistic actor attribute model; ALAAM; social network analysis; social contagion; network autocorrelation; stochastic approximation; equilibrium expectation; ALAAMEE; MPNet

**Language.** English

**Related identifiers.** (relation, identifier, type)

| Relation | Identifier | Type |
|---|---|---|
| Is supplement to | https://github.com/thomas-u-grund/nwalaam/releases/tag/v1.2.0 | URL |
| Is supplement to | https://github.com/thomas-u-grund/nwalaam | URL |
| Requires | https://github.com/thomas-u-grund/nwcommands | URL |
| Is related to (manual of the sibling package nwsaom) | 10.5281/zenodo.23623807 | DOI |
| Is related to (manual of the sibling package nwergm) | 10.5281/zenodo.23661331 | DOI |

The first identifier points to the GitHub release v1.2.0: confirm that this release is published before the upload (if the package has a later 1.2.x release by then, use its tag). The two sibling links are optional.

**Other fields.** Funding: none. Communities: none.

## Steps for the owner

1. Read the PDF once (`manual/nwalaam_manual.pdf`, 42 pages); the examples were produced with Stata 17 on macOS with version 1.2.0 of the package.
2. On zenodo.org create the new record, reserve the DOI, put it into `\manualdoi`, rebuild (`make -C manual pdf`), upload the PDF and fill the fields above.
3. Commit the rebuilt PDF (`git add -f manual/nwalaam_manual.pdf`).
4. Afterwards: the DOI goes into the citation texts (README of the package, `CITATION.cff`, the website page *How to cite* of nwalaam and the manual link of its menu), as for nwergm; none of these is changed yet.
