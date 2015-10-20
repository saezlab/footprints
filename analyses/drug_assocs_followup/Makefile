include ../../Makefile.inc
	
all: SKCM_assocs.pdf COREAD_assocs.pdf

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', 'pdf_document', output_file='$@')"
