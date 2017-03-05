include ../../Makefile.inc

all: stability_matrix.pdf stability_matrix_robust.pdf gene_overlap.pdf

.SECONDARY:

resample_%.RData: model_resample.r ../../data/zscores.RData ../../model/model_%.r
	Rscript $^ $@

gdsc_%.RData: score_gdsc.r resample_%.RData
	Rscript $^ $@

stability_%.pdf: stability_plot.r gdsc_%.RData
	Rscript $^ $@

gene_overlap.pdf: gene_overlap.r resample_matrix.RData
	Rscript $^ $@
