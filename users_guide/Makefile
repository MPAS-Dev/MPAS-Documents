all:
	pdflatex mpas_$(CORE)_users_guide
	bibtex mpas_$(CORE)_users_guide
	pdflatex mpas_$(CORE)_users_guide
	pdflatex mpas_$(CORE)_users_guide
clean:
	rm -f *.aux *.dvi *.pdf *.log *.out *.toc *.bbl *.blg
	rm -f shared/*.aux
ifdef CORE
	rm -f $(CORE)/*.aux
endif

