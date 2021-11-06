# Makefile for compiling RMarkdown documents
# Document name (no file extension)
DOC = dissertation
DIFFDOC = dissertation_v2
DIFF = diff

# Strings to search for in log file for rerunning pdflatex until all
# cross references and other issues are resolved
RERUN := "(undefined references|Rerun to get (cross-references|the bars|point totals) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed|has been referenced but does not)"

# Same for bibliography
RERUNBIB := "No file.*\.bbl|Citation.*undefined"

# Top level
all: doc diffdoc diff

# Doc requires the pdf file
doc: $(DOC).pdf
	@open $<

diffdoc: $(DIFFDOC).pdf
	@open $<

diff: $(DIFF).pdf
	@open $<

# Main level
# pdf file requires compiling the tex file
# using pdflatex and bibtex
$(DOC).pdf: $(DOC).tex
	@echo "\n======= Running PDFLaTeX =======\n"
	@pdflatex $<
	# if the document has a glossary
	@makeglossaries $(DOC)
	# if the document has a nomenclature
	@makeindex $(DOC).nlo -s nomencl.ist -o $(DOC).nls
	# Rerun bibtex until there are no issues
	@while(egrep -q $(RERUNBIB) $(DOC).log && bibtex $(DOC) && pdflatex $<) do \
		echo "\n======= Rerunning BibTeX =======\n" ; done
	# Rerun pdflatex until there are no issues
	@while(egrep -q $(RERUN) $(DOC).log && pdflatex $<) do \
		echo "\n======= Rerunning PDFLaTeX =======\n" ; done

$(DIFFDOC).pdf: $(DIFFDOC).tex
	@echo "\n======= Running PDFLaTeX =======\n"
	@pdflatex $<
	# if the document has a glossary
	@makeglossaries $(DIFFDOC)
	# if the document has a nomenclature
	@makeindex $(DIFFDOC).nlo -s nomencl.ist -o $(DIFFDOC).nls
	# Rerun bibtex until there are no issues
	@while(egrep -q $(RERUNBIB) $(DIFFDOC).log && bibtex $(DIFFDOC) && pdflatex $<) do \
		echo "\n======= Rerunning BibTeX =======\n" ; done
	# Rerun pdflatex until there are no issues
	@while(egrep -q $(RERUN) $(DIFFDOC).log && pdflatex $<) do \
		echo "\n======= Rerunning PDFLaTeX =======\n" ; done

$(DIFF).pdf: $(DIFF).tex
	@pdflatex $<

# Lowest level
# tex document requires knitting by calling
# rmarkdown::render in R
$(DOC).tex: $(DOC).Rmd
	@R -e 'rmarkdown::render("$(DOC).Rmd")'

$(DIFFDOC).tex: $(DIFFDOC).Rmd
	@R -e 'rmarkdown::render("$(DIFFDOC).Rmd")'

$(DIFF).tex: $(DOC).tex $(DIFFDOC).tex
	@latexdiff -Vt CFONT $(DOC).tex $(DIFFDOC).tex > $(DIFF).tex
	@sed -i '' 's:\\usepackage{mathptmx}:\\usepackage{lmodern}:g' $(DIFF).tex

# pdflatex command
pdflatex:
	pdflatex --shell-escape -interaction=nonstopmode -file-line-error -synctex=1 %O %S

# purge auxillary files
purge:
	@rm -f *.{acr,alg,gz,glg,gls,ilg,nls,acn,glo,ist,lof,lot,nlo,aux,dvi,log,bbl,blg,brf,fls,toc,thm,out,fdb_latexmk}

# clean directory (excluding .Rmd file)
clean: purge
	@rm -f $(DOC).pdf $(DOC).tex $(DIFFDOC).pdf $(DIFFDOC).tex $(DIFF).tex $(DIFF).pdf

.PHONY: all diff purge clean