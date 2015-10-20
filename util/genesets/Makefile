all: gatza.RData go.RData biocarta.RData reactome.RData chea.RData

gatza.RData: gatza2014.r ng.3073-S2.csv
	Rscript $^ $@

go.RData: go.r
	Rscript $^ $@

biocarta.RData: biocarta.r Enrichr/src/main/resources/BioCarta.gmt
	Rscript $^ $@

reactome.RData: reactome.r ReactomePathways.gmt
	Rscript $^ $@

chea.RData: chea.r Enrichr/src/main/resources/ChEA.gmt
	Rscript $^ $@

ReactomePathways.gmt:
	wget -N http://www.reactome.org/download/current/ReactomePathways.gmt.zip
	unzip ReactomePathways.gmt.zip
	rm ReactomePathways.gmt.zip
