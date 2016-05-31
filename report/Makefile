LATEXMK = latexmk -xelatex #-c # no pdf created with -c
RNWs = $(sort $(wildcard *.Rnw))

FIGS = $(filter-out S_%, $(RNWs))
SUPP = $(filter S_%, $(RNWs))

#include ../Makefile.inc

#.INTERMEDIATE:


all: Figures.pdf Supplementary.pdf


Figures.pdf: Figures.tex config.tex $(FIGS)
	$(LATEXMK) Figures

Supplementary.pdf: Supplementary.tex config.tex $(SUPP)
	$(LATEXMK) Supplementary


%.tex: %.Rnw
	R -e 'library(knitr);knit("$<")'
