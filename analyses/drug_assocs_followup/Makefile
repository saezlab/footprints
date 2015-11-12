include ../../Makefile.inc

RMDS = $(wildcard *.Rmd)
PDFS = $(RMDS:%.Rmd=%.pdf)

all: $(PDFS)

%.pdf: %.Rmd util.r
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
