include ../../Makefile.inc

all: geneset_overlap.pdf bars.pdf model_significance.pdf

geneset_overlap.pdf: geneset_overlap.r
	Rscript $< $@

numbers.txt: numbers.r
	Rscript $< $@

bars.pdf: plot_bars.r numbers.txt
	Rscript $< $@

model_significance.pdf: model_significance.r
	Rscript $< $@

treemap.pdf: treemap.r
	Rscript $< $@
