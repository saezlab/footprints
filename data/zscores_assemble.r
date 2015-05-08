b = import('base')
io = import('io')
ar = import('array')

# read individual z-score files
OUTFILE = commandArgs(TRUE)[1]
INFILES = commandArgs(TRUE)[-1] %or% c("zscores/H2O2/E-GEOD-47739.RData",
                                       "zscores/MAPK/E-GEOD-14934.RData")

# put together data matrix + index df
contents = lapply(INFILES, io$load)

zscores = lapply(contents, function(x) x$zscores) %>% ar$stack(along=2)
index = lapply(contents, function(x) x$index) %>% do.call(rbind, .)

save(zscores, index, file=OUTFILE)
