# this only works if scores target was run before
GDSC_SCORES = $(shell find ../../scores/gdsc -name "*.RData")
GDSC_ASSOCS = $(GDSC_SCORES:../../scores/gdsc/pathways_%.RData=assocs_%.RData)
GDSC_PLOTS = $(GDSC_ASSOCS:assocs_%.RData=plots_%.pdf)

include ../../Makefile.inc
$(call ext_dep,../../scores/gdsc,$(notdir $(GSCORES)))
$(call ext_dep,../../scores/ccle,$(notdir $(CSCORES)))


all: drug_ranges.pdf $(GDSC_PLOTS) assocs_mapped/mutation.RData drug_ranges.pdf

# use each model script to create a model file
$(GDSC_PLOTS): plots_%.pdf: plot.r assocs_%.RData
	Rscript $^ $@

$(GDSC_ASSOCS): assocs_%.RData: gdsc_assocs.r ../../scores/gdsc/pathways_%.RData
	$(call bsub_log,10240) Rscript $^ $@

mutation.RData: mutation.r
	Rscript $^ $@

assocs_mapped/mutation.RData: mutation.RData
	ln -s $^ $@

drug_ranges.pdf: drug_ranges.r
	Rscript $^ $@
