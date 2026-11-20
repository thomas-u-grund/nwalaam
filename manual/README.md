# The nwalaam manual

`nwalaam_manual.pdf` is the user manual and reference of nwalaam (title "Manual", no version number in the citation). The LaTeX source is `nwalaam_manual.tex` with one file per chapter (`ch*.tex`, `app_translation.tex`); the listings of the examples are `ex/out/*.txt`.

## Bumping the version (one-line change)

The manual's title and citation carry no version number; the title page states the series it describes (`\manualseries`, valid for 1.3.x). The exact version and date of the release that produced the listings are two more macros in `nwalaam_manual.tex`, written by `tools/set_version.py` together with the tag in the installation sentence of `ch01_intro.tex`:

    \newcommand{\manualseries}{1.3}
    \newcommand{\pkgversion}{1.3.0}
    \newcommand{\pkgdate}{20 November 2026}

Change `\manualseries` only for a new minor or major release whose documented behaviour changed. For a patch release rebuild (`make -C manual pdf`) and rerun the examples if their output changed.

## The DOI

`\manualdoi` is empty until the Zenodo record exists; an empty DOI prints nothing in the citation. `ZENODO.md` holds the metadata of the record and the steps to do.

## Building

    make -C manual pdf        # pdflatex, 3 passes
    make -C manual check      # undefined references, overfull boxes
    NWCORE=/path/to/nwcommands make -C manual examples   # reruns the do-files of ex/ in Stata 17 (several minutes; Glasgow, Lazega and Davis come from the web)

The do-files `ex/ex01_karate.do` to `ex/ex04_special.do` write one text log per listing into `ex/out` with `log using`; `ex/mklistings.py` removes the log header and footer lines. The karate club data are in `ex/data` (edges and attributes from the ALAAMEE repository, MIT licence). The listings with random numbers use fixed seeds; the two timings in `ex/out/lazega_threads.txt` depend on the computer and its load.
