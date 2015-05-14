all: gatza.RData go.RData

gatza.RData: gatza2014.r ng.3073-S2.csv
	Rscript $< $(word 2,$^) $@

go.RData: go.r
	Rscript $< $@
