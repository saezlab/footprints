include ../../Makefile.inc

SCORES = $(wildcard ../../scores/tcga/*.RData)
STEMS = $(filter-out speed_fullmat.RData,$(SCORES:../../scores/tcga/%.RData=%.RData))
PLOTS = $(STEMS:%.RData=box+cor/%.pdf) $(STEMS:%.RData=heatmap/%.pdf)

$(call ext_dep,../../scores/tcga,$(STEMS))
$(call ext_dep,../../scores/merge,$(STEMS))

all: $(PLOTS) tact_merge.pdf

box+cor/%.pdf: tissue_act.r ../../scores/tcga/%.RData
	Rscript $^ $@

heatmap/%.pdf: pancan_heatmap.r ../../scores/tcga/%.RData
	Rscript $^ $@

tact_merge.pdf: tact_merge.r ../../scores/merge/speed_linear.RData
	Rscript $^ $@
