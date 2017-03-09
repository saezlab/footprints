include ../../Makefile.inc

all: bootstrap_matrix.pdf bootstrap_matrix_robust.pdf gene_overlap.pdf

.SECONDARY:

resample_%.RData: model_resample.r ../../data/zscores.RData ../../model/model_%.r
	Rscript $^ $@

gdsc_%.RData: score_gdsc.r resample_%.RData
	Rscript $^ $@

bootstrap_%.pdf: bootstrap_plot.r gdsc_%.RData
	Rscript $^ $@

gene_overlap.pdf: gene_overlap.r resample_matrix.RData
	Rscript $^ $@
