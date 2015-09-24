SCORES = $(wildcard ../../scores/speed/*.RData)
HEATMAPS = $(SCORES:../../scores/speed/%.RData=heatmap_%.pdf)

plot: $(HEATMAPS)

heatmap_%.pdf: ../../scores/speed/%.RData heatmap.r
	Rscript heatmap.r $< $@
