all: fdr_de.RData top100_de.RData top100_z.RData

fdr_de.RData top100_de.RData: %.RData: %.r ../../data/expr.RData
	Rscript $^ $@

top100_z.RData: %.RData: %.r ../../data/zscores.RData
	Rscript $^ $@
