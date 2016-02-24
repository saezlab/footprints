all: tsne_fullmat.RData tsne_top100.RData

tsne_fullmat.RData: tsne_speed.r ../../scores/speed/speed_fullmat.RData
	Rscript $< $@

tsne_top100.RData: tsne_speed.r ../../scores/speed/speed_matrix.RData
	Rscript $< $@
