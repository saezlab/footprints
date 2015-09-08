SCORES = $(wildcard ../../scores/speed/*.RData)
PLOTS = $(SCORES:../../scores/speed/%.RData=%.pdf) zscores.pdf

plot: $(PLOTS)

zscores.pdf: zscores.r
	Rscript $<

%.pdf: ../../scores/speed/%.RData heatmap.r
	Rscript heatmap.r $< $@
