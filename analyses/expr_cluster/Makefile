all: drug_speed.pdf drug_expr.pdf surv_speed.pdf surv_expr.pdf tsne_both.pdf tsne_gdsc.pdf

drug_%.pdf: ../../expr_cluster/%_cluster.RData drug_assocs.r
	Rscript drug_assocs.r $< $@

surv_%.pdf: ../../expr_cluster/%_cluster.RData surv_assocs.r
	Rscript surv_assocs.r $< $@

tsne_%.pdf: tsne_%.r
	Rscript $<
