b = import('base')
io = import('io')
ar = import('array')

# read individual z-score files
OUTFILE = commandArgs(TRUE)[1]
INFILES = commandArgs(TRUE)[-1] %or% c("zscores/H2O2/E-GEOD-47739.RData",
                                       "zscores/MAPK/E-GEOD-14934.RData")

# put together data matrix + index df
contents = io$load(INFILES)

zscores = lapply(contents, function(x) x$zscores) %>% ar$stack(along=2)
colnames(zscores) = 1:ncol(zscores)

index = lapply(contents, function(x) {
    as.data.frame(x$index, stringsAsFactors=FALSE) %>%
        lapply(unlist) %>%
        as.data.frame(stringsAsFactors=FALSE)
}) %>% do.call(rbind, .)
rownames(index) = 1:nrow(index)

save(zscores, index, file=OUTFILE)
