include ../../Makefile.inc
	
TCGA_SCORES = $(wildcard ../../scores/tcga/*.RData)
SURV = $(TCGA_SCORES:../../scores/tcga/%.RData=%.pdf)
CONT = $(SURV:%=cont_%)
DISCR = $(SURV:%=discr_%)

all: prereqs $(CONT) $(DISCR)

prereqs:;
#	make -C ../../scores/tcga

discr_%.pdf: survival_discrete.r ../../scores/tcga/%.RData
	Rscript $^ $@

cont_%.pdf: survival_cont.r ../../scores/tcga/%.RData
	Rscript $^ $@
