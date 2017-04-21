all: expr.RData plot.pdf hek_expr.pdf hek_phospho.pdf

expr.RData: expr.r validation.yaml
	Rscript expr.r

plot.pdf: do.r expr.RData
	Rscript do.r

hek_expr.pdf: hek_expr.r speed2_raw_counts_hek293.tsv ../../model/model_matrix.RData
	Rscript $<

hek_phospho.pdf: hek_phospho.r signalling_data_log2fc_to_bsa_ctrl.txt
	Rscript $<
