include ../../Makefile.inc

SCORES = $(wildcard ../../scores/tcga/*.RData)
STEMS = $(filter-out speed_fullmat.RData,$(SCORES:../../scores/tcga/%.RData=%.RData))
PLOTS = $(STEMS:%.RData=box+cor/%.pdf) $(STEMS:%.RData=heatmap/%.pdf)

$(call ext_dep,../../scores/tcga,$(STEMS))
$(call ext_dep,../../scores/merge,$(STEMS))

all: $(PLOTS) tact_merge.pdf tissue_act.pdf

box+cor/%.pdf: tissue_act.r ../../scores/tcga/pathways_mapped/%.RData
	Rscript $^ $@

heatmap/%.pdf: heatmap.r ../../scores/tcga/pathways_mapped/%.RData
	Rscript $^ $@

tact_merge.pdf: tact_merge.r ../../scores/merge/speed_matrix.RData
	Rscript $^ $@

tissue_act.pdf: tissue_act.r ../../scores/tcga/pathways_mapped/speed_matrix.RData
	Rscript $^ $@
