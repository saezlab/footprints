TCGA_SCORES = $(wildcard ../../scores/tcga/*.RData)
SURV = $(TCGA_SCORES:../../scores/tcga/%.RData=%.pdf)
CONT = $(SURV:%=cont_%)
DISCR = $(SURV:%=discr_%)

include ../../Makefile.inc
$(call ext_dep,../../scores/tcga,$(notdir $(CONT) $(DISCR)))


all: $(CONT) $(DISCR)

discr_%.pdf: survival_discrete.r ../../scores/tcga/%.RData
	Rscript $^ $@

cont_%.pdf: survival_cont.r ../../scores/tcga/%.RData
	Rscript $^ $@
