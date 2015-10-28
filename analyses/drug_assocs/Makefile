# this only works if scores target was run before
GSCORES = $(wildcard ../../scores/gdsc/*.RData)
GPLOTS = $(GSCORES:../../scores/gdsc/%.RData=gdsc/%.pdf)
CSCORES = $(wildcard ../../scores/ccle/*.RData)
CPLOTS = $(CSCORES:../../scores/ccle/%.RData=ccle/%.pdf)

include ../../Makefile.inc
$(call ext_dep,../../scores/gdsc,$(notdir $(GSCORES)))
$(call ext_dep,../../scores/ccle,$(notdir $(CSCORES)))


all: $(GPLOTS) #$(CPLOTS)

# use each model script to create a model file
gdsc/%.pdf: drug_assocs.r ../../scores/gdsc/%.RData
	Rscript $< $(word 2,$^) $@

ccle/%.pdf: drug_assocs.r ../../scores/ccle/%.RData
	Rscript $< $(word 2,$^) $@
