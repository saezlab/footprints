SUBS = fdr_de.RData top100_de.RData top100_z.RData speed_matrix.RData

all: tsne_df.RData

fdr_de.RData top100_de.RData: %.RData: %.r ../../data/expr.RData
	Rscript $^ $@

top100_z.RData: %.RData: %.r ../../data/zscores.RData
	Rscript $^ $@

speed_matrix.RData: %.RData: speed.r ../../scores/speed/%.RData
	Rscript $^ $@

tsne_df.RData: tsne_df.r $(SUBS)
	Rscript $< $@ $(SUBS)
