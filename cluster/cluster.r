library(dplyr)
b = import('base')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')
icgc = import('data/icgc')

TISSUE = commandArgs(TRUE)[1] %or% "BRCA"
OUTFILE = commandArgs(TRUE)[2] %or% "BRCA.RData"

# load TCGA, GDSC data for given tissue
avail = icgc$available(clinical=TRUE, rna_seq=TRUE, map_to="icgc_specimen_id")
clinical = icgc$clinical(avail) %>%
    select(icgc_specimen_id, project_code, tissue) %>% #TODO: separate cancer/normal
    filter(tissue == TISSUE)
tcga_expr = icgc$rna_seq(clinical$icgc_specimen_id, voom=TRUE)
tcga_expr = tcga_expr[,clinical$icgc_specimen_id] #TODO: duplicates, ordering @icgc module

tissues = gdsc$tissues(TISSUE)
gdsc_expr = gdsc$basal_expression()
gdsc_expr = gdsc_expr[,intersect(names(tissues), colnames(gdsc_expr))]

# perform batch correction
genes = intersect(rownames(tcga_expr), rownames(gdsc_expr))
expr = na.omit(cbind(tcga_expr[genes,], gdsc_expr[genes,]))
batches = c(clinical$project_code, rep("cell_line", ncol(gdsc_expr)))
corrected = st$batch$combat(expr, batches)

# nmf-cluster them together, get optimal number of clusters
# couple of hours w/o hpc
clust = st$nmf(corrected, k=2:10, max_iter=5000, rep=10) %>%
    arrange(k, cluster)
coph = select(clust, k, rho) %>% unique()
n = coph$k[which(diff(coph$rho) > 0) + 1]

# create ar$mask for tcga and gdsc separately
n2mask = function(n) {
    df = clust %>% filter(k == n) %>% select(sample, cluster)
    setNames(df$cluster, df$sample) %>% ar$mask() + 0
}
masks = lapply(n, n2mask)

# save object with matrix for TCGA, GDSC
tumors = lapply(masks, x -> x[grepl("^SP", rownames(x)),]) %>% ar$stack(along=2)
tumors = tumors[,colSums(tumors) != 0]
clines = lapply(masks, x -> x[!grepl("^SP", rownames(x)),]) %>% ar$stack(along=2)
clines = clines[,colSums(clines) != 0]
save(tumors, clines, file=OUTFILE)
