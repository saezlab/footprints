include ../../../Makefile.inc

TISSUES = $(shell ../../../config/get.sh tcga/tissues)

all: $(TISSUES:%=path_%.out)

expr: $(TISSUES:%=paradigm_%_expr.txt)

proteins.txt:
	grep ^protein < SuperPathway.txt | grep -v abstract | cut -f2 > $@

paradigm_%_expr.txt: preprocess_expr.r proteins.txt
	Rscript preprocess_expr.r $* $@

path_%.out: paradigm_%_expr.txt
	$(call bsub,10240) -oo $@ paradigm -m 10 -c config.txt -p SuperPathway.txt -b paradigm_$*
