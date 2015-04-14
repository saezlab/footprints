#!/usr/bin/env Rscript
library(RSQLite)
library(GEOquery)
library(modules)
library(optparse)
gn = import('general')
ar = import('array')
ma = import('microarray')

# takes a comma-separated string and explodes it to a vector
explC = function(gstr) unlist(strsplit(na.omit(gstr), split=","))

# takes a character vector of GSM identifiers and returns expression data
getSoftExpression = function(index) {
    pair = unique(index[c('GSE','GPL')])
    lexpr = lapply(1:nrow(pair), function(i) {
                   print(i)
        gsms = index[index$GSE == pair[i,'GSE'] & index$GPL == pair[i,'GPL'], 'GSM']
        ma$GSMsToExpression(gsms, cachedir=file.path(module_file(), "soft_cache"))
    } %catch% NA) #TODO: illumina arrays+single array failures

    ar$stack(lexpr[!is.na(lexpr)], along=2)
}

if (is.null(module_name())) {
    optionList = list(
        make_option("--db", help="SPEED.db to use", default="../SPEED-API/other_scripts/SPEED1.db"),
        make_option("--outfile", help=".RData to save expr object to")
    )
    opt = parse_args(OptionParser(option_list=optionList))

    # create db connection
    db = dbConnect(SQLite(), dbname=opt$db)
    exps = dbGetQuery(db, "SELECT * FROM Experiments")
    index = gn$read.table("../SPEED-Data/index.txt")
    index = index[index$SUBSET %in% exps$e_id,]

    # normalize and HGNC-map arrays
    expr = getSoftExpression(index)
    expr[is.nan(expr)] = NA

    save(expr, file=opt$outfile)
}

