library(dplyr)
b = import('base')
io = import('io')
df = import('data_frame')

fnames = c("top100_z.RData", "top100_de.RData", "fdr_de.RData", "speed_matrix.RData")
file2name = setNames(c("Top 100 z-score", "Top 100 DE", "10% FDR", "speed_matrix"), fnames)

INFILES = unname(commandArgs(TRUE)[1:4] %or% fnames)
OUTFILE = commandArgs(TRUE)[4] %or% "tsne_df.RData"

file2tsne = function(fname) {
    mat = t(io$load(fname))
    mat = mat[!duplicated(mat),]
    pts = as.data.frame(Rtsne::Rtsne(mat)$Y)
    colnames(pts) = c('x', 'y')
    pts$pathway = sub("\\..*$", "", rownames(mat))
    pts
}

tsne = b$lnapply(fnames, file2tsne) %>%
    df$add_name_col("method", bind=TRUE) %>%
    mutate(method = file2name[method])

save(tsne, file=OUTFILE)
