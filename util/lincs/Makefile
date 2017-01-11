TXT = $(wildcard *.txt)

all: index.RData

index.RData: index.r $(TXT)
	Rscript index.r $@ $^
