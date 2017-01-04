include ../../Makefile.inc

all: geneset_overlap.pdf dataset_size_bars.pdf gene_significance.pdf treemap.pdf

geneset_overlap.pdf: geneset_overlap.r
	Rscript $< $@

dataset_size.txt: dataset_size.r
	Rscript $< $@

dataset_size_bars.pdf: dataset_size_bars.r dataset_size.txt
	Rscript $< $@

gene_significance.pdf: model_significance.r
	Rscript $< $@

treemap.pdf: treemap.r
	Rscript $< $@
