library(biomaRt)
library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
df = import('data_frame')
hpc = import('hpc')

# get index, expr data for test set
zfit = io$load("../../model/model_fullmat.RData")

# load gene lists for pathways
mapGO = biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") %>%
    biomaRt::getBM(attributes=c("hgnc_symbol", "go_id", "name_1006"),
                   filter="hgnc_symbol", values=rownames(zfit), mart=.)

# convert df that we get from biomart into list of vectors
genesets = mapGO %>%
    filter(go_id != "" & hgnc_symbol != "") %>%
    mutate(go_id=paste(go_id,name_1006)) %>%
    dplyr::select(hgnc_symbol, go_id) %>%
    group_by(go_id) %>%
    filter(n() >= 5 & n() <= 200) %>%
    ungroup() %>%
    unstack()

# run GSEA on fitted z-scores
do_gsea = function(zcol, set, zfit, genesets) {
    gsea = import('../../util/gsea')
    gsea$wGSEA(zfit[,zcol], genesets[[set]], significance=TRUE)
}
idx = df$create_index(zcol=colnames(zfit), set=names(genesets), expand_grid=TRUE,
                      args=list(zfit=zfit, genesets=genesets))

scores = df$call(idx, do_gsea, hpc_args=list(memory = 1024, n_jobs = 200)) %>%
    group_by(zcol) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup()

#xx = scores %>% filter(abs(NES) >=2 & adj.p < 1e-4) %>% transmute(zcol=zcol, NES=NES, adj.p=adj.p, set=set) %>% arrange(zcol, adj.p) %>% write.table("gocats.txt", sep="\t", quote=F)
