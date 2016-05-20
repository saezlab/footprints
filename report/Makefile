RMDs = $(wildcard *.Rmd)

include ../Makefile.inc

all: $(RMDs:%.Rmd=%.pdf)

Supplementary.pdf: Supplementary.Rmd
	Rscript -e "bookdown::render_book('$<', 'bookdown::pdf_book')"
	mv _book/_main.pdf $@
	rm -rf _book _bookdown_files

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
