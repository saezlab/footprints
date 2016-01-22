include ../../../Makefile.inc

# TODO: pull this in over config/yaml
TISSUES = ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD

all: $(TISSUES:%=path_%.out)

proteins.txt:
	grep ^protein < SuperPathway.txt | grep -v abstract | cut -f2 > $@

paradigm_%_expr.txt: proteins.txt
	Rscript preprocess_expr.r $* $@

path_%.out: paradigm_%_expr.txt
	$(call bsub_log,10240) paradigm -m 10 -c config.txt -p SuperPathway.txt -b $*
