SPEED = $(wildcard ../../scores/speed/*.RData)
LINCS = $(wildcard ../../scores/lincs/*.RData)
HEATMAPS = $(SPEED:../../scores/speed/%.RData=heatmap_%.pdf)

plot: $(HEATMAPS) roc.pdf roc_lincs.pdf

heatmap_%.pdf: ../../scores/speed/%.RData heatmap.r
	Rscript heatmap.r $< $@

roc.pdf: roc.r $(SCORES) ../../scores/speed/speed_matrix.RData
	Rscript $< $@

roc_lincs.pdf: roc_lincs.r $(LINCS) ../../scores/lincs/speed_matrix.RData
	Rscript $< $@
