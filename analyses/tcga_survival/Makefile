TCGA_SCORES = $(shell find ../../scores/tcga -name "*.RData")
SURV = $(TCGA_SCORES:../../scores/tcga/pathways_%.RData=%.pdf)
CONT = $(SURV:%=cont/%)
DISCR = $(SURV:%=discr/%)

include ../../Makefile.inc
#$(call ext_dep,../../scores/tcga,$(notdir $(CONT) $(DISCR)))


all: $(CONT) $(DISCR)

discr/%.pdf: survival_discrete.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@

cont/%.pdf: survival_cont.r ../../scores/tcga/pathways_%.RData
	@mkdir -p $(dir $@)
	Rscript $^ $@
