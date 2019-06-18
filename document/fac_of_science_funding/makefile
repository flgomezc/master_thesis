
FIGDIR := figures/

all: document.pdf

document.pdf: document.tex ${FIGDIR}
	pdflatex document.tex
	pdflatex document.tex

view:
	evince document.pdf

${FIGDIR}:
	mkdir figures

clean:
	rm -f document.pdf
	rm -f document.log
	rm -f document.aux
