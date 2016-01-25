include ../../../Makefile.inc

TISSUES = $(shell ../../../config/get.sh tcga/tissues)
SAMPLES = $(TISSUES:%=%/.dirstamp)

all: $(SAMPLES)  #$(TISSUES:%=path_%.out)

expr: $(TISSUES:%=paradigm_%_expr.txt)

proteins.txt:
	grep ^protein < SuperPathway.txt | grep -v abstract | cut -f2 > $@

# this could be improved by:
# http://www.gnu.org/savannah-checkouts/gnu/make/manual/html_node/Prerequisite-Types.html
$(SAMPLES): %/.dirstamp: preprocess_expr.r proteins.txt
	mkdir -p $*
	Rscript preprocess_expr.r $*
	touch $*/.dirstamp

path_%.out: paradigm_%_expr.txt
	$(call bsub,10240) -oo $@ paradigm -m 10 -c config.txt -p SuperPathway.txt -b paradigm_$*
