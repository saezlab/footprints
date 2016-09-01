include ../../Makefile.inc

all: assocs.RData

assocs.RData: assocs.r ../../scores/tcga/pathways_mapped/speed_matrix.RData
	Rscript $^ $@

correlation.RData: correlation.r ../../scores/tcga/pathways_mapped/speed_matrix.RData
	Rscript $^ $@
