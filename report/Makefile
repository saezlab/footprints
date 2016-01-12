RMDs = $(wildcard *.Rmd)

include ../Makefile.inc

all: $(RMDs:%.Rmd=%.pdf)

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
