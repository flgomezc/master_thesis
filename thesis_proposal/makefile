FIGDIR := figures/

all: document.tex ${FIGDIR}
	rm -f document.pdf
	pdflatex document.tex
	pdflatex document.tex
	rm -f document.log
	rm -f document.aux
	evince document.pdf

${FIGDIR}:
	mkdir figures
