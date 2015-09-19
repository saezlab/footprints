SCORES = $(wildcard ../../scores/speed_nocv/*.RData)
HEATMAPS = $(SCORES:../../scores/speed_nocv/%.RData=heatmap_%.pdf)

plot: $(HEATMAPS)

heatmap_%.pdf: ../../scores/speed_nocv/%.RData heatmap.r
	Rscript heatmap.r $< $@
