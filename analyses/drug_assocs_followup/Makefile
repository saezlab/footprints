include ../../Makefile.inc

RMDS = $(wildcard *.Rmd)
PDFS = $(RMDS:%.Rmd=%.pdf)

all: $(PDFS)

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
