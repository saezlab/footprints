tissue2scores = function(tissue, genesets, INDEX) {
    io = import('io')
    lincs = import('data/lincs')
    pathifier = import_package('pathifier')

    # prepare data
    index = io$load(INDEX)
    cid_control = unique(index[index$pathway == "control",]$distil_id)
    cid_perturbed = unique(index[index$pathway == tissue,]$distil_id)
    tumors = lincs$get_z(cid=cid_perturbed, rid=lincs$projected, map.genes="hgnc_symbol")
    normals = lincs$get_z(cid=cid_control, rid=lincs$projected, map.genes="hgnc_symbol")
    data = cbind(tumors, normals)

    result = pathifier$quantify_pathways_deregulation(
        data = data,
        allgenes = rownames(tumors),
        syms = genesets,
        pathwaynames = names(genesets),
        normals = c(rep(FALSE, ncol(tumors)), rep(TRUE, ncol(normals))),
        attempts = 100,
        min_exp = -Inf,
        min_std = 0.4
    )

    re = do.call(cbind, lapply(result$scores, c))
    rownames(re) = colnames(data)
    re
}

b = import('base', attach_operators=FALSE)
import('base/operators')
io = import('io')
ar = import('array')
hpc = import('hpc')

INFILE = commandArgs(TRUE)[1] %or% "../../util/genesets/reactome.RData"
INDEX = "../../util/lincs_perturbation_qc/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "pathifier.RData"

# load pathway gene sets
genesets = io$load(INFILE)
pathways = intersect(names(genesets), unique(io$load(INDEX)$pathway))

# run pathifier in jobs
result = hpc$Q(tissue2scores, tissue=pathways,
    const=list(INDEX=INDEX, genesets=genesets), memory=8192, n_jobs=50)

result = ar$stack(result, along=3) #TODO: check why this is required
result = ar$map(result, along=3, function(x) mean(x,na.rm=TRUE))

# save results
save(result, file=OUTFILE)
