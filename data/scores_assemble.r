b = import('base')
io = import('io')
ar = import('array')

# read individual z-score files
TYPE = commandArgs(TRUE)[1] %or% "z"
OUTFILE = commandArgs(TRUE)[2]
INFILES = commandArgs(TRUE)[c(-1,-2)] %or% c("scores/H2O2/E-GEOD-47739.RData",
                                             "scores/MAPK/E-GEOD-14934.RData")

# put together data matrix + index df
contents = io$load(INFILES)

scores = lapply(contents, function(x) x[[paste0(TYPE,"scores")]]) %>% ar$stack(along=2)
colnames(scores) = 1:ncol(scores)

index = lapply(contents, function(x) {
    as.data.frame(x$index, stringsAsFactors=FALSE) %>%
        lapply(unlist) %>%
        as.data.frame(stringsAsFactors=FALSE)
}) %>% do.call(rbind, .)
rownames(index) = 1:nrow(index)

save(scores, index, file=OUTFILE)
