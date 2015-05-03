io = import('io')

# read individual z-score files
OUTFILE = commandArgs(TRUE)[1]
INFILES = commandArgs(TRUE)[-1]

# put together data matrix + index df
contents = io$load()

zscores = lapply(contents, function(x) x$zscores) %>% cbind()
index = lapply(contents, function(x) x$index) %>% rbind()

save(zscores, index, file=OUTFILE)
