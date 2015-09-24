b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "spia.RData"

tissue2scores = function(tissue, spia) {
    io = import('io')
    tcga = import('data/tcga')

    # convert hgnc to entrez
    expr = tcga$rna_seq(tissue)

    # HGNC -> entrez gene lookup
    lookup = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
        biomaRt::getBM(attributes=c("hgnc_symbol", "entrezgene"),
        filter="hgnc_symbol", values=rownames(expr), mart=.)

    rownames(expr) = lookup$entrezgene[match(rownames(expr), lookup$hgnc_symbol)]
    expr = limma::avereps(expr[!is.na(rownames(expr)),])

    is_normal = grepl("[Nn]ormal", tcga$barcode2index(colnames(expr))$Sample.Definition)
    tumors = expr[,!is_normal]
    normals = expr[,is_normal]

    spia$spia(tumors, normals, per_sample=TRUE, pathids=spia$speed2kegg, verbose=TRUE)
}

# load pathway gene sets
genesets = io$load(INFILE)
tissues = c("BLCA", "BRCA", "CESC", "COREAD", "ESCA", "HNSC",
            "KIRC", "LIHC", "LUAD", "LUSC", "PAAD")

# run spia in jobs and save
result = hpc$Q(tissue2scores, tissue=tissues,
    const=list(spia=spia), memory=8192, n_jobs=length(tissues))

result = ar$stack(result, along=1)
colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
