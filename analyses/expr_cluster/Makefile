all: speed_cluster.pdf expr_cluster.pdf

%_cluster.pdf: ../../gene_expr/%_cluster.RData drug_assocs.r
	Rscript drug_assocs.r $< $@
