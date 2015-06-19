all: drug_speed.pdf drug_expr.pdf surv_speed.pdf surv_expr.pdf

drug_%.pdf: ../../gene_expr/%_cluster.RData drug_assocs.r
	Rscript drug_assocs.r $< $@

surv_%.pdf: ../../gene_expr/%_cluster.RData surv_assocs.r
	Rscript surv_assocs.r $< $@
