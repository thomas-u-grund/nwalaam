# Build of the manual.  make pdf     : the PDF is built (3 passes)
#                      make examples: reruns the do-files of ex/ in Stata 17 (needs NWCORE=/path/to/nwcommands; about 5 minutes; Glasgow, Lazega and Davis are loaded from the web)
#                      make check   : undefined references, overfull boxes
# The PDF is committed with git add -f (the LaTeX build files are ignored).
all: pdf

pdf:
	pdflatex -interaction=nonstopmode nwalaam_manual.tex
	pdflatex -interaction=nonstopmode nwalaam_manual.tex
	pdflatex -interaction=nonstopmode nwalaam_manual.tex

examples:
	./run_examples.sh

check:
	@grep -n "undefined\|Undefined\|Overfull\|Citation.*undefined" nwalaam_manual.log || echo "no undefined references, no overfull boxes"

clean:
	rm -f *.aux *.log *.toc *.out *.bbl *.blg *.lof *.lot
