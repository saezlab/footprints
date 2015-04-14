# tissue-specific EN
#
# ...
library(modules)

feat2score = function(targets, feats, tissues, randomize=F, bootstrap=F, ...) { # ...: bjw fix
    require(modules)
    require(glmnet)
    b = import('base', attach_operators=F)
    ar = import('array')
    an = import('anova')

    f2s_tissue_drug = function(target, fmat) {
        # only look at samples where response variable is available
        if (sum(!is.na(target)) < 10)
            return(NULL)
        
        fmat = fmat[!is.na(target),]
        fmat = fmat[,apply(fmat, 2, function(x) length(unique(x))!=1)]
        target = target[!is.na(target)]

        if (randomize)
            target = setNames(sample(target, size=length(target), replace=T), names(target))

        mod = cv.glmnet(fmat, c(target), alpha=1, nfolds=5)

        # look at models with lambda.min or more regularized
        weight = as.matrix(coef(mod, s=mod$lambda))
        nfeat = unname(colSums(weight[colnames(fmat),] != 0))
        valid = mod$lambda >= mod$lambda.min

        # best models per number of features
        best.per.nfeat = ar$map(mod$cvm, along=1, function(x) x==min(x), subsets=nfeat)
        best.per.nfeat[!valid | is.na(best.per.nfeat) | nfeat==0] = F

        # compare predicted vs observed
        if (any(best.per.nfeat)) {
            pred = predict(mod, s=mod$lambda[best.per.nfeat], newx=fmat)
            rsq = an$calcAssocs(pred, target, sumsq=T)

            # remove intercept and all-0 weights
            weight = weight[-1,best.per.nfeat,drop=F]
            weight = weight[rowSums(abs(weight))!=0,,drop=F]
            colnames(weight) = nfeat[best.per.nfeat]

            if (randomize) # limit size of background models
                weight = NA

            list(nfeat = nfeat[best.per.nfeat],
                 weights = weight,
                 rsquared = setNames(rsq$main, nfeat[best.per.nfeat]))
        } else NULL
    }

    # assemble array, split feats by tissues, build model for each
    f2s_tissue = function(targets, feats, randomize) {
        targets = ar$split(targets, along=2)
        b$omit$null(lapply(targets, function(t) f2s_tissue_drug(t, feats)))
    }

    targets = ar$split(targets, along=1, subsets=tissues)
    feats = ar$split(feats, along=1, subsets=tissues)
    mapply(f2s_tissue, targets, feats, MoreArgs=list(randomize=randomize))
#TODO: maybe stack here already, to put more load on the nodes as well?
}

rsq2array = function(result) {
    b = import('base')
    ar = import('array')

    ac = x -> as.character(x)
    fr = function(x) { # fill right: rsq(2 feat) >= rsq(1 feat)
        for (i in 2:length(x))
            if (is.na(x[i]))
                x[i] = x[i-1]
        x
    }

    repl2rsquared = rsq -> ar$stack(lapply(rsq,
        repl -> ar$stack(lapply(repl,
        tissue -> do.call(rbind, lapply(tissue,
        drug -> fr(setNames(drug$rsquared[ac(1:4)], 1:4))))))))
    repl2rsquared(result)
}

if (is.null(module_name())) {
    library(optparse)
    optionList = list(
        make_option("--randomize", help="whether or not to randomize drug response", default=F),
        make_option("--bootstrap", help="whether or not to bootstrap", default=F),
        make_option("--num", help="how often to build models", default=1),
        make_option("--sets", help=paste("'-' separated list of input feature names",
                    "including: pathways, mut, drivers"), default="pathways"),
        make_option("--outfile", help="file to save results to", default="out.RData")
    )
    opt = parse_args(OptionParser(option_list=optionList))
    opt$sets = strsplit(opt$sets, "-")[[1]]

    ar = import('array')
    io = import('io')
    sg = import('sanger_robject')
    bj = import('hpc/BatchJobsWrapper')

    feats = list(tissues = sg$getTissues(minN=10),
                 drug = sg$getDrugResponseForCellLines('IC50s')) # or AUC

    if ("pathways" %in% opt$sets)
        feats$pathways = io$load('pathways.RData')
    if ("mut" %in% opt$sets)
        feats$mut = sg$getMutatedGenes()
    if ("adrivers" %in% opt$sets)
        feats$drivers = sg$getMutatedGenes(intogen_all=T)
#    if ("tdrivers" %in% opt$sets)
#        feats$drivers = sg$getMutatedGenes(intogen_tissue=T) # add to feats after split

    feats = ar$intersect_list(feats, along=1)
    targets = feats$drug
    tissues = feats$tissues
    feats$drug = feats$tissues = NULL

    feats = do.call(cbind, feats)

    result = bj$Q(feat2score, 1:opt$num, more.args=list(targets=targets, feats=feats,
                  tissues=tissues, randomize=opt$randomize, bootstrap=opt$bootstrap))

    if (opt$randomize) # cut down file size
        result = rsq2array(result)

    save(result, file=opt$outfile)
}
