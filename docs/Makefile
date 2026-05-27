
MANUAL = shyfem_main042025

default:

manual: docs

docs:
	pdflatex $(MANUAL).tex
	bibtex $(MANUAL)
	pdflatex $(MANUAL).tex
	pdflatex $(MANUAL).tex



clean:
	-rm -f *.aux *.dvi *.log *.toc *.bbl *.blg
	-rm -f S_*.tex P_*.tex
	-rm -f S_*.txt P_*.txt
	-rm -f $(DOC)*.html
	-rm -f femdoc.zip
	-rm -f ggg hhh
	-rm -f *.tmp *.bak
	-rm -f *~

cleanps:
	-rm -f *.ps *.pdf

#--------------------------------

zip:	cleanall
	zip -r femdoc *

save:	zip

cleanall: clean
	-rm -f *.pdf

