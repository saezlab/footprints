SCORES = $(wildcard ../../scores/speed_test/*.RData)
HEATMAPS = $(SCORES:../../scores/speed_test/%.RData=heatmap_%.pdf)

plot: $(HEATMAPS)

heatmap_%.pdf: ../../scores/speed_test/%.RData heatmap.r
	Rscript heatmap.r $< $@
