TISSUES = HNSC ESCA LAML DLBC LCML MM ALL BRCA COREAD LIHC STAD KIRC LUAD LUSC SCLC GBM LGG NB PAAD SKCM THCA BLCA CESC OV

all: $(addsuffix .RData,$(TISSUES))

%.RData: cluster.r
	bsub -K -M 20480 -R "rusage[mem=20480]" -oo $(@:%.RData=%.log) \
		-R 'select[panfs_nobackup_research]' \
		Rscript $< $* $@
