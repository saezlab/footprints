KEGG = list(
    MAPK = "04010", # MAPK signaling pathway
    EGFR = "04012", # ErbB signaling pathway
#    PPAR = "03320", # PPAR signaling pathway
    PI3K = "04150", # mTOR signaling pathway
    Wnt = "04310", # Wnt signaling pathway
#    Notch = "04330", # Notch signaling pathway
#    Insulin = "04910", # Insulin signaling pathway
    NFkB = "04064", # NF-kappa B signaling pathway
    VEGF = "04370", # VEGF signaling pathway
    `JAK-STAT` = "04630", # Jak-STAT signaling pathway
    TGFB = "04350", # TGF-beta signaling pathway
    Trail = "04210" # Apoptosis
)

tissue2scores = function(tissue, lookup, INDEX, EXPR, KEGG) {
    b = import('base')
    io = import('io')
    ar = import('array')
    spia = import_package('SPIA')

    # prepare data
    index = io$load(INDEX)
    expr = io$load(EXPR)
# tissues = pathway perturbation subsets
    tissues = index$pathway
# tumors = perturbed
    tumors = expr[,index$sign == "+" & index$pathway == tissue]
# normals = control
    normals = expr[,index$sign == "0" & index$pathway == tissue]
    data = cbind(tumors, normals)

    # map IDs to entrez gene IDs, SPIA needs this
    rownames(data) = lookup$entrezgene[match(rownames(data), lookup$hgnc_symbol)]
    data = limma::avereps(data[!is.na(rownames(data)),])

    sample2score = function(sample) {
        # create data for each sample
        message(sample)
        data2 = cbind(data[,sample,drop=FALSE], data[,colnames(normals)])
        types = c("tumor", rep("normal", ncol(normals)))

        # compute differential expression between tumor and normal
        design = model.matrix(~ 0 + as.factor(types))
        colnames(design) = sub("as\\.factor\\(types\\)", "", colnames(design))
        contrast = limma::makeContrasts("tumor-normal", levels=design)
        fit = limma::lmFit(data2, design) %>%
            limma::contrasts.fit(contrast) %>%
            limma::eBayes()
        result = dplyr::data_frame(gene = rownames(fit),
            p.value = fit$p.value[,1],
            fold_change = fit$coefficients[,1]) %>%
            dplyr::mutate(adj.p = p.adjust(p.value, method="fdr"))

        # calculate SPIA scores (relevant fields: Name [KEGG name], ID [KEGG ID], tA [score])
        re = spia$spia(
            de = setNames(result$fold_change, result$gene), #[result$adj.p < 0.1],
            all = result$gene,
            organism = "hsa",
            pathids = unlist(KEGG, use.names=FALSE),
            nB = 2000, # bootstraps
            plots = FALSE,
            verbose = FALSE
        ) %catch% data.frame()

        KEGGi = setNames(names(KEGG), KEGG)
        setNames(re$tA, KEGGi[re$ID])
    }
    
    sapply(colnames(tumors), sample2score, simplify=FALSE, USE.NAMES=TRUE) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale)
}

if (is.null(module_name())) {
    b = import('base')
    io = import('io')
    ar = import('array')
    hpc = import('hpc')

#    INFILE = commandArgs(TRUE)[1] %or% "../../genesets/reactome.RData"
    INDEX = commandArgs(TRUE)[1] %or% "../../data/lincs_perturbation_qc/index.RData"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/lincs_perturbation_qc/expr.RData"
    OUTFILE = commandArgs(TRUE)[3] %or% "spia.RData"

    # load pathway gene sets
#    genesets = io$load(INFILE)

#    # get all tissues which have a normal
#    tissues = io$h5load(EXPR, "/tissue")
#    tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))
    tissues = unique(io$load(INDEX)$pathway)

    # HGNC -> entrez gene lookup
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
        filter="hgnc_symbol", values=rownames(io$load(EXPR)), mart=.)

    # run spia in jobs
    result = hpc$Q(tissue2scores, tissue=tissues,
        more.args=list(INDEX=INDEX, EXPR=EXPR, KEGG=KEGG, lookup=lookup), memory=1024)

    result = ar$stack(result, along=1)

    # save results
    save(result, file=OUTFILE)
}
