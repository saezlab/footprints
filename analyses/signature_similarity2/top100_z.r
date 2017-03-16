b = import('base')
io = import('io')
ar = import('array')

ZDATA = commandArgs(TRUE)[1] %or% "../../data/zscores.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "top100_z.RData"

zdata = io$load(ZDATA)

index = zdata$index
mat = t(t(zdata$zscores) * index$sign)
mat[apply(mat, 2, function(p) !b$min_mask(abs(p), 100))] = 0

save(mat, file=OUTFILE)
