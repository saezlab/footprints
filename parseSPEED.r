# get a list: entrez ids/zvalues from all output.out
library(parallel)
library(modules)
gn = import('general')
ar = import('array')
bm = import('biomart')

# define z-values text files
basedir = '../DATA/RAW/Transcriptomic/SPEED/INDEX'
dirs = list.files(path=basedir, pattern="^GSE", recursive=T, include.dirs=T)
zvals = file.path(basedir, dirs, "zvalues.dat")
descs = file.path(basedir, dirs, "Description.txt")
sets = file.path(basedir, dirs, "set.dat")

# put together meta-data data frame
desc = t(ar$stack(lapply(descs, function(f) gn$read.vector(f, sep="\t"))))
set = do.call(rbind, lapply(sets, gn$read.vector))
bto = gn$read.vector(paste(basedir, "index_bto.txt", sep="/"), sep="\t")
meta = data.frame(pathway = dirs %|% "grep -oP ^[^/]+", 
                  cells = desc[,'cell line'],
                  BTO = names(bto)[match(desc[,'cell line'], bto)],
                  treatment = desc[,'modification'],
                  effect = desc[,'effect'],
                  time = desc[,'time'],
                  comment = desc[,'comment'],
                  GSE = set[,1], 
                  norm = set[,2], 
                  GPL = set[,3], 
                  row.names = dirs %|% "grep -oP GSE[^/]+")

# load combine them to matrix (TODO?: handle multiple occurrences of gene IDs)
read2zvals = function(fname) {
    tab = read.table(fname, sep="\t", strip.white=T, stringsAsFactors=F)
    re = tab$zvalue
    names(re) = tab$GeneID
    re[names(re) %|% "grep -P ^[0-9]+$"]
}
scorelist = mclapply(zvals, read2zvals)
names(scorelist) = rownames(meta)

# map HGNC symbols to expression values for all experiments (TODO?:mapping for each patform enough)
allentrez = unique(unlist(mclapply(scorelist, names)))
hgnc = bm$getHGNCByEntrez(allentrez)
scorelist = mclapply(scorelist, function(f) bm$entrez2hgnc(f, hgnc))

# get mean pathway activation vectors
scores = array_stack(scorelist)
save(scores=scores, meta=meta, file="zval.ro")
write.table(meta, file="zval_meta.txt", sep="\t", quote=F)

