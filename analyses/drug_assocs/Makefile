# this only works if scores target was run before
GDSC_SCORES = $(wildcard ../../scores/gdsc/*.RData)
GDSC_ASSOCS = $(GDSC_SCORES:../../scores/gdsc/%.RData=%.RData)
GDSC_PLOTS = $(GDSC_SCORES:../../scores/gdsc/%.RData=%.pdf) mutation.pdf
#CSCORES = $(wildcard ../../scores/ccle/*.RData)
#CPLOTS = $(CSCORES:../../scores/ccle/%.RData=ccle_%.pdf)

include ../../Makefile.inc
$(call ext_dep,../../scores/gdsc,$(notdir $(GSCORES)))
$(call ext_dep,../../scores/ccle,$(notdir $(CSCORES)))


all: $(GDSC_PLOTS) #$(CPLOTS)

# use each model script to create a model file
$(GDSC_PLOTS): %.pdf: plot.r %.RData
	Rscript $^ $@

$(GDSC_ASSOCS): %.RData: gdsc_assocs.r ../../scores/gdsc/%.RData
	Rscript $^ $@

mutation.RData: mutation.r
	Rscript $^ $@

#ccle_%.pdf: drug_assocs.r ../../scores/ccle/%.RData
#	Rscript $^ $@
