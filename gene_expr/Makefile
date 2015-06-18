TISSUES = BLCA BRCA CESC COREAD GBM HNSC KIRC LAML LGG LIHC LUAD LUSC OV PAAD SKCM THCA

all: corrected_expr.h5 $(TISSUES:%=cluster_%.RData)

cluster_%.RData: cluster.r
	bsub -K -M 10240 -R "rusage[mem=10240]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
		Rscript $< $* $@

corrected_expr.h5: expr.r
	Rscript expr.r $@
