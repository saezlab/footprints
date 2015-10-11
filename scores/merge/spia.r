b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
EXPR = commandArgs(TRUE)[2] %or% "../../util/expr_cluster/corrected_expr.h5"
OUTFILE = commandArgs(TRUE)[3] %or% "spia.RData"

tissue2scores = function(tissue, EXPR, spia, lookup) {
    io = import('io')

    # convert hgnc to entrez
    e = function(data) {
        rownames(data) = lookup$entrezgene[match(rownames(data), lookup$hgnc_symbol)]
        limma::avereps(data[!is.na(rownames(data)),])
    }

    tissues = io$h5load(EXPR, "/tissue")
    tumors = e(t(io$h5load(EXPR, "/expr", index=which(tissues == tissue))))
    normals = e(t(io$h5load(EXPR, "/expr", index=which(tissues == paste0(tissue, "_N")))))

    spia$spia(tumors, normals, per_sample=TRUE, pathids=spia$speed2kegg, verbose=TRUE)
}

# load pathway gene sets and tissues
EXPR = tools:::file_path_as_absolute(EXPR)
genesets = io$load(INFILE)
tissues = io$h5load(EXPR, "/tissue")
tissues = sub("_N", "", unique(tissues[grepl("_N", tissues)]))

# HGNC -> entrez gene lookup
lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
    filter="hgnc_symbol", values=io$h5names(EXPR, "/expr")[[2]], mart=.)

# run spia in jobs and save
result = hpc$Q(tissue2scores, tissue=tissues,
    const=list(EXPR=EXPR, spia=spia, lookup=lookup), memory=8192, n_jobs=length(tissues))

result = ar$stack(result, along=1)
colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
