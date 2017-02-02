include ../../Makefile.inc

all: heatmap.pdf

heatmap.pdf: heatmap.r assocs.RData
	Rscript $^ $@

assocs.RData: assocs.r ../../scores/tcga/pathways_mapped/speed_matrix.RData
	Rscript $^ $@
