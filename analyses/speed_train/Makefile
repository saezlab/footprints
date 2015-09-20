SCORES = $(wildcard ../../scores/speed_train/*.RData)
HEATMAPS = $(SCORES:../../scores/speed_train/%.RData=heatmap_%.pdf)

plot: $(HEATMAPS)

heatmap_%.pdf: ../../scores/speed_train/%.RData heatmap.r
	Rscript heatmap.r $< $@
