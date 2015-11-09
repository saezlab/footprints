include ../../Makefile.inc

RMDS = $(wildcard *.Rmd)
PDFs = $(RMDS:%.Rmd=%.pdf)

all: $(PDFS)

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
