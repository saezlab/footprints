SCORES = $(wildcard ../../scores/speed/*.RData)
PLOTS = $(SCORES:../../scores/speed/%.RData=%.pdf)

plot: $(PLOTS)

%.pdf: ../../scores/speed/%.RData heatmap.r
	Rscript heatmap.r $< $@
