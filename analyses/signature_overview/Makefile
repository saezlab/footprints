include ../../Makefile.inc

all: geneset_overlap.pdf geneset_overlap_venn.pdf dataset_size_bars.pdf gene_significance.pdf treemap.pdf

%.pdf: %.r
	Rscript $< $@

dataset_size_bars.pdf: dataset_size_bars.r dataset_size.txt
	Rscript $< $@
