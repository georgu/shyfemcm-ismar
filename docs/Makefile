
MANUAL = shyfem_main042025

default:

manual: docs

docs:
	pdflatex $(MANUAL).tex
	bibtex $(MANUAL)
	pdflatex $(MANUAL).tex
	pdflatex $(MANUAL).tex



clean:
	-rm -f *.dvi
	-rm -f *.blg *.bbl
	-rm -f *.toc *aux
	-rm -f *.log

cleanall: clean
	-rm -f *.pdf

