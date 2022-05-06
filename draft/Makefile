# Makefile for compiling RMarkdown documents
# Document name (no file extension)
DOC = dissertation_v1-2-3
REVISION = dissertation_v1-2-4
DIFF = dissertation_diff

# pdflatex command
LTXFLAGS = "--shell-escape -interaction=nonstopmode -file-line-error -synctex=1 %O %S"
# Strings to search for in log file for rerunning pdflatex until all
# cross references and other issues are resolved
RERUN := "(Rerun to get (cross-references|the bars|point totals) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed|has been referenced but does not)"
RERUNDIFF := "(Rerun to get (cross-references|the bars|point totals) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed|has been referenced but does not)"
# Same for bibliography
RERUNBIB := "No file.*\.bbl|Citation.*undefined"

# Top level
all: doc

# Doc requires the pdf file
doc: $(DOC).pdf
	@open $<

revision: $(REVISION).pdf
	@open $<

diff: $(DIFF).pdf
	@open $<

# Main level
# pdf file requires compiling the tex file
# using pdflatex and bibtex
$(DOC).pdf: $(DOC).tex
	@printf "\nPDFLaTeX for %s ..." "$(DOC)"
	@pdflatex $(LTXFLAGS) $< > /dev/null 2>&1
	@echo "\nMaking Glossaries ..."
	@makeglossaries $(DOC) > /dev/null 2>&1
	@makeindex $(DOC).nlo -s nomencl.ist -o $(DOC).nls > /dev/null 2>&1
	@while(egrep -q $(RERUNBIB) $(DOC).log > /dev/null 2>&1 && bibtex $(DOC) > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "\nRerunning BibTeX ..." ; done
	@while(egrep $(RERUN) $(DOC).log > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "Rerunning PDFLaTeX ..." ; done

$(REVISION).pdf: $(REVISION).tex
	@printf "\nPDFLaTeX for %s ..." "$(REVISION)"
	@pdflatex $(LTXFLAGS) $< > /dev/null 2>&1
	@echo "\nMaking Glossaries ..."
	@makeglossaries $(REVISION) > /dev/null 2>&1
	@makeindex $(REVISION).nlo -s nomencl.ist -o $(REVISION).nls > /dev/null 2>&1
	@while(egrep -q $(RERUNBIB) $(REVISION).log > /dev/null 2>&1 && bibtex $(REVISION) > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "\nRerunning BibTeX ..." ; done
	@while(egrep $(RERUN) $(REVISION).log > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "Rerunning PDFLaTeX ..." ; done

$(DIFF).pdf: $(DIFF).tex
	@printf "\nPDFLaTeX for %s ..." "$(DIFF)"
	@pdflatex $(LTXFLAGS) $< > /dev/null 2>&1
	@echo "\nMaking Glossaries ..."
	@makeglossaries $(DIFF) > /dev/null 2>&1
	@makeindex $(DIFF).nlo -s nomencl.ist -o $(DIFF).nls > /dev/null 2>&1
	@while(egrep -q $(RERUNBIB) $(DIFF).log > /dev/null 2>&1 && bibtex $(DIFF) > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "\nRerunning BibTeX ..." ; done
	@while(egrep $(RERUNDIFF) $(DIFF).log > /dev/null 2>&1 && pdflatex $(LTXFLAGS) $< > /dev/null 2>&1) do \
		echo "Rerunning PDFLaTeX ..." ; done

# Lowest level
# tex document requires knitting by calling
# rmarkdown::render in R
$(DOC).tex: $(DOC).Rmd
	@printf "\nRmarkdown Render for %s ..." "$(DOC)"
	@R -q -e 'rmarkdown::render("$(DOC).Rmd", quiet = T)' > /dev/null 2>&1
	@rm -f $(DOC).pdf

$(REVISION).tex: $(REVISION).Rmd
	@printf "\nRmarkdown Render for %s ..." "$(REVISION)"
	@R -q -e 'rmarkdown::render("$(REVISION).Rmd", quiet = T)' > /dev/null 2>&1
	@rm -f $(REVISION).pdf

$(DIFF).tex: $(DOC).tex $(REVISION).tex
	@echo "\nLaTeX Diff for %s ..." "$(DIFF)"
	@latexdiff -V -t CFONT $(DOC).tex $(REVISION).tex > $(DIFF).tex
	@sed -i '' 's:\\usepackage{mathptmx}:\\usepackage{lmodern}:g' $(DIFF).tex

# purge auxillary files
purge:
	@rm -f *.{acr,alg,gz,glg,gls,ilg,nls,acn,glo,ist,lof,lot,nlo,aux,dvi,log,bbl,blg,brf,fls,toc,thm,out,fdb_latexmk}

# clean directory (excluding .Rmd file)
clean: purge
	@rm -f $(DOC).pdf $(DOC).tex $(REVISION).pdf $(REVISION).tex $(DIFF).tex $(DIFF).pdf

.PHONY: all diff doc revision purge clean