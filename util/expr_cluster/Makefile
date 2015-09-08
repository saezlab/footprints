include ../../Makefile.inc

TISSUES = BLCA BRCA CESC COREAD GBM HNSC KIRC LAML LGG LIHC LUAD LUSC OV PAAD SKCM THCA

all: corrected_expr.h5 speed_cluster.RData expr_cluster.RData

expr_cluster.RData: cluster_expr_assemble.r $(TISSUES:%=.cluster_%.RData)
	Rscript $< $@

speed_cluster.RData: cluster_speed.r
	Rscript $< ../../model/model_linear.RData $@

.cluster_%.RData: cluster_expr.r
	$(call bsub_log,10240) Rscript $< $* $@

corrected_expr.h5: expr.r
	$(call bsub_log,10240) Rscript expr.r $@
