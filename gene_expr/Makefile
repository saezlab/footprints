TISSUES = BLCA BRCA CESC COREAD GBM HNSC KIRC LAML LGG LIHC LUAD LUSC OV PAAD SKCM THCA

all: corrected_expr.h5 $(TISSUES:%=cluster_%.RData)

cluster_%.RData: cluster.r
	bsub -K -M 102400 -R "rusage[mem=102400]" -R "select[panfs_nobackup_research]" -oo $(@:%.RData=%.log) \
		Rscript $< $* $@

corrected_expr.h5: expr.r
	Rscript expr.r $@
