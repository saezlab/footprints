LATEXMK = latexmk -xelatex -c
RNWs = $(sort $(wildcard *.Rnw))

FIGS = $(patsubst %.Rnw, %.tex, $(filter-out S_%, $(RNWs)))
SUPP = $(patsubst %.Rnw, %.tex, $(filter S_%, $(RNWs)))

include ../Makefile.inc

.INTERMEDIATE:


all: Figures.pdf Supplementary.pdf


Figures.pdf: Figures.tex $(FIGS)
	$(LATEXMK) Figures

Supplementary.pdf: Supplementary.tex $(SUPP)
	$(LATEXMK) Supplementary


%.tex: %.Rnw
	R -e 'library(knitr);knit("$<")'
