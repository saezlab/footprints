include ../../Makefile.inc
	
TCGA_SCORES = $(wildcard ../../scores/tcga/*.RData)
SURV = $(TCGA_SCORES:../../scores/tcga/%.RData=%.pdf)
CONT = $(SURV:%=cont_%)
DISCR = $(SURV:%=discr_%)

all: prereqs $(CONT) $(DISCR)

prereqs:;
#	make -C ../../scores/tcga

discr_%.pdf: ../../scores/tcga/%.RData survival_discrete.r
	Rscript survival_discrete.r $< $@

cont_%.pdf: ../../scores/tcga/%.RData survival_cont.r
	Rscript survival_cont.r $< $@
