all: gatza.RData go.RData biocarta.RData

gatza.RData: gatza2014.r ng.3073-S2.csv
	Rscript $< $(word 2,$^) $@

go.RData: go.r
	Rscript $< $@

biocarta.RData: biocarta.r Enrichr/src/main/resources/BioCarta.gmt
	Rscript $< $(word 2,$^) $@
