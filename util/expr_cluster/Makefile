include ../../Makefile.inc

TISSUES = $(shell ../../config/get.sh merge/tissues)

all: corrected_expr.RData speed_cluster.RData expr_cluster.RData

expr_cluster.RData: cluster_expr_assemble.r $(TISSUES:%=.cluster_%.RData)
	Rscript $< $@

speed_cluster.RData: cluster_speed.r
	Rscript $< ../../model/model_linear.RData $@

.cluster_%.RData: cluster_expr.r
	$(call bsub_log,10240) Rscript $< $* $@

corrected_expr.RData: expr.r
	Rscript expr.r $@
