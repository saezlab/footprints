TISSUES = BLCA BRCA CESC COREAD GBM HNSC KIRC LAML LGG LIHC LUAD LUSC OV PAAD SKCM THCA

all: corrected_expr.h5 speed_cluster.RData expr_cluster.RData

expr_cluster.RData: cluster_expr_assemble.r $(TISSUES:%=.cluster_%.RData)
	Rscript $< $@

speed_cluster.RData: cluster_speed.r
	Rscript $< ../model/model_linear.RData $@

.cluster_%.RData: cluster_expr.r
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
		Rscript $< $* $@

corrected_expr.h5: expr.r
	Rscript expr.r $@
