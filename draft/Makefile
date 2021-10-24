# Makefile for compiling RMarkdown documents
# Document name (no file extension)
DOC := dissertation

# Strings to search for in log file for rerunning pdflatex until all
# cross references and other issues are resolved
RERUN := "(undefined references|Rerun to get (cross-references|the bars|point totals) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed|has been referenced but does not)"
# Same for bibliography
RERUNBIB := "No file.*\.bbl|Citation.*undefined"

# Top level
all: doc

# Doc requires the pdf file
doc: $(DOC).pdf
	open $<

# Main level
# pdf file requires compiling the tex file
# using pdflatex and bibtex
$(DOC).pdf: $(DOC).tex
	@echo "======= Running PDFLaTeX ======="
	@pdflatex $<
	# if the document has a glossary
	@makeglossaries $(DOC)
	# if the document has a nomenclature
	@makeindex $(DOC).nlo -s nomencl.ist -o $(DOC).nls
	# Rerun bibtex until there are no issues
	@while(egrep -q $(RERUNBIB) $(DOC).log && bibtex $(DOC) && pdflatex $<) do \
		echo "======= Rerunning BibTeX =======" ; done
	# Rerun pdflatex until there are no issues
	@while(egrep -q $(RERUN) $(DOC).log && pdflatex $<) do \
		echo "======= Rerunning PDFLaTeX =======" ; done

# Lowest level
# tex document requires knitting by calling
# rmarkdown::render in R
$(DOC).tex: $(DOC).Rmd
	@R -e 'rmarkdown::render("$(DOC).Rmd")'

# pdflatex command
pdflatex:
	pdflatex --shell-escape -interaction=nonstopmode -file-line-error -synctex=1 %O %S

# purge auxillary files
purge:
	@rm -f *.{acr,alg,glg,gls,ilg,nls,acn,glo,ist,lof,lot,nlo,aux,dvi,log,bbl,blg,brf,fls,toc,thm,out,fdb_latexmk}

# clean directory (excluding .Rmd file)
clean: purge
	@rm -f $(DOC).pdf $(DOC).tex

.PHONY: all purge clean