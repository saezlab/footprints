LATEXMK = latexmk -f -xelatex #-c # no pdf created with -c
RNWs = $(sort $(wildcard *.Rnw))

FIGS = $(filter-out S_%, $(RNWs))
SUPP = $(filter S_%, $(RNWs))

#include ../Makefile.inc

#.INTERMEDIATE:


all: Figures.pdf Supplementary.pdf PROGENy_model.xlsx

clean: 
	rm -f {Figures,Supplementary}.{tex,aux,bbl,bcf,blg,fdb_latexmk,fls,log,out,run.xml,lof,lot,toc,nlo}
	rm -f references.bib

Figures.pdf: Figures.tex config.tex $(FIGS)
	$(LATEXMK) Figures

Supplementary.pdf: Supplementary.tex config.tex $(SUPP)
	$(LATEXMK) Supplementary

PROGENy_model.xlsx: excel_model.r ../model/model_matrix.RData
	Rscript $^

%.pdf: %.tex config.tex
	$(LATEXMK) $*

%.tex: %.Rnw
	R -e 'library(knitr);knit("$<")'
