dplyr = import_package('dplyr')
b = import('base')
io = import('io')
ar = import('array')
spia = import('../../util/spia')
hpc = import('hpc')

INDEX = commandArgs(TRUE)[1] %or% "../../util/lincs_perturbation_qc/index.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "spia.RData"

pathway2scores = function(pathway, index, spia) {
    lincs = import('data/lincs')

    cid_control = unique(index[index$pathway == "control",]$distil_id)
    cid_perturbed = unique(index[index$pathway == pathway,]$distil_id)
    perturbed = lincs$get_z(cid=cid_perturbed, rid=lincs$projected, map.genes="entrezgene")
    control = lincs$get_z(cid=cid_control, rid=lincs$projected, map.genes="entrezgene")

    spia$spia(perturbed, control, per_sample=TRUE, pathids=spia$speed2kegg, verbose=TRUE)
}

# get index and pathways
index = io$load(INDEX)
pathway = intersect(names(spia$speed2kegg), unique(index$pathway))

# run spia in jobs and save
result = hpc$Q(pathway2scores, pathway=pathway,
    const=list(spia=spia, index=index), memory=8192, n_jobs=50)

result = ar$stack(result, along=1)
colnames(result) = spia$kegg2speed[colnames(result)]

save(result, file=OUTFILE)
