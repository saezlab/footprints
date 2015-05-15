all: speed_linear.RData speed_matrix.RData

speed_linear.RData: speed.r ../../model/model_linear.RData
	Rscript $< $(word 2,$^) $@

speed_matrix.RData: speed.r ../../model/model_matrix.RData
	Rscript $< $(word 2,$^) $@
