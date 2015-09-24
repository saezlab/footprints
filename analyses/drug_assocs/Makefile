include ../../Makefile.inc

# this only works if scores target was run before
GSCORES = $(wildcard ../../scores/gdsc/*.RData)
GPLOTS = $(GSCORES:../../scores/gdsc/%.RData=gdsc/%.pdf)
CSCORES = $(wildcard ../../scores/ccle/*.RData)
CPLOTS = $(CSCORES:../../scores/ccle/%.RData=ccle/%.pdf)

.PHONY: prereqs

all: prereqs $(GPLOTS) #$(CPLOTS)

prereqs: ;
#	make -C ../../scores/gdsc $(notdir $(GSCORES))
##	make -C ../../scores/ccle $(notdir $(CSCORES))

# use each model script to create a model file
gdsc/%.pdf: drug_assocs.r ../../scores/gdsc/%.RData
	Rscript $< $(word 2,$^) $@

ccle/%.pdf: drug_assocs.r ../../scores/ccle/%.RData
	Rscript $< $(word 2,$^) $@
