
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
MANUAL = shyfem_manual

default:

manual: doc

#--------------------------------

doc:
	./femdoc.sh
	pdflatex $(MANUAL).tex
	bibtex $(MANUAL)
	pdflatex $(MANUAL).tex
	pdflatex $(MANUAL).tex


#--------------------------------

clean:
	-rm -f *.dvi
	-rm -f *.blg *.bbl
	-rm -f *.toc *aux
	-rm -f *.log
	-rm -f S_*.tex P_*.tex

cleanall: clean
	-rm -f *.pdf

