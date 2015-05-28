all: gatza.RData go.RData biocarta.RData reactome.RData

gatza.RData: gatza2014.r ng.3073-S2.csv
	Rscript $< $(word 2,$^) $@

go.RData: go.r
	Rscript $< $@

biocarta.RData: biocarta.r Enrichr/src/main/resources/BioCarta.gmt
	Rscript $< $(word 2,$^) $@

reactome.RData: reactome.r ReactomePathways.gmt
	Rscript $< $(word 2,$^) $@

ReactomePathways.gmt:
	wget -N http://www.reactome.org/download/current/ReactomePathways.gmt.zip
	unzip ReactomePathways.gmt.zip
	rm ReactomePathways.gmt.zip
