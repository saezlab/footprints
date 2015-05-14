ar = import('array')
icgc = import('icgc')

# load expression, RPPA for all cancers where both available
avail = icgc$dataAvailable(RNASeq=T, RPPA=T, map_to="specimen")
expr = icgc$getRNASeq(specimen=avail, voom=TRUE)
prot = icgc$getRPPA(specimen=avail)
ar$intersect(expr, prot, along=2)

# calculate gatza, speed scores for expression


# plot linear fits w/ covar = exp.series for both, see what difference
