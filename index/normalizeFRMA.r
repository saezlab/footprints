#!/usr/bin/env Rscript

###
### command-line arguments
###

library(optparse)
optionList = list(
    make_option("--indexfile", help="Index file for .CEL data"),
    make_option("--outfile", help="RObject file to save the result in"),
    make_option("--platform", help="The GPLxxx string for which platform FRMA should be performed"),
    make_option("--normalize", help="Normalisation: none, quantile"),
    make_option("--summarize", help="Summarisation: median_polish, average, median, weighted_average, robust_weighted_average, random_effect)")
)
opt = parse_args(OptionParser(option_list=optionList))

###
### load files
###

# load index
fc = read.table(opt$indexfile, sep="\t", as.is=T)
index = fc[-1,-1]
rownames(index) = fc[-1,1]
colnames(index) = as.vector(t(fc[1,-1])) # whyy??
rm(fc)

# choose appropriate array subsets
gsms = index[index$GPL==opt$platform & !is.na(index$FILE),]
gsms = gsms[!duplicated(gsms$FILE),]

# read into affy batch object
if (opt$platform == 'GPL6244') {
    library(oligo)
    dataObj = read.celfiles(gsms$FILE)
} else {
    library(affy)
    dataObj = ReadAffy(filenames=gsms$FILE, verbose=TRUE)
}

###
### frma-normalize 
###

library(frma)
object = frma(dataObj, 
              summarize=opt$summarize, 
              background="rma", 
              normalize=opt$normalize, 
              verbose=T)

# make sure order is the same as we put in
stopifnot(colnames(pData(dataObj)) == gsms$FILE)

# save exprs only
object = exprs(object)
colnames(object) = gsms$GSM
save(object, file=opt$outfile)

