library(modules)
library(pheatmap)
b = import('base')
io = import('io')
gs = import('genesets')
hpc = import('hpc')

# plot heatmap of enrichment scores
drawHeatmap = function(full_mat, method) {
    mat = full_mat[,dimnames(full_mat)[[2]]==method,]
    rownames(mat) = sapply(rownames(mat), function(x) substr(x, 1, 90))

    lx = -log(mat)
    lx[lx==Inf] = 10
    rslx = rowSums(lx)
    lx = lx[rslx>=b$maxN(rslx, 200),]

    pheatmap(lx, fontsize_row=4, fontsize_col=6, main=method,
        breaks=c((1:90)/20, log(1:11)*(max(lx)+0.1-5)/log(11)+5))
}

if (is.null(module_name())) {
    # run enrichment
    mat = hpc$Q(piano::runGSA,
                geneLevelStats = io$load('models_4-8h.RData')[["50"]]$x,
                more.args = list(gsc = gs$piano('GO_projected', 'TF_iorio', 'gatza'),
                                 gsSizeLim = c(5,500)),
                memory=1024) %>% gs$piano_result2matrix()

    pdf("factor_enrichment.pdf", paper="a4", width=20, height=27)
    for (method in c('pAdjNonDirectional','pAdjDistinctDirUp','pAdjDistinctDirDn'))
        try(drawHeatmap(mat, method))
    dev.off()
}
