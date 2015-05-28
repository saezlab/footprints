io = import('io')
st = import('stats')
plt = import('plot')
icgc = import('data/icgc')

INFILE = "../../scores/tcga/speed_linear.RData"
OUTFILE = "speed_linear.pdf"

# load expression, RPPA for all cancers where both available
avail = icgc$available(clinical=TRUE, rna_seq=TRUE, rppa=TRUE, map_to="specimen")
clinical = icgc$clinical(specimen=avail)
prot = t(icgc$rppa(specimen=avail))

scores = io$load(INFILE)[avail,]
prot = prot[avail,] # get rid of duplicates

# plot linear fits w/ covar = exp.series for both, see what difference
#TODO: linear fits for pathway<-->respective phospho
mapk = scores[,"MAPK"]
mek = prot[,"MEK1_pS217_S221"] #/ prot[,"MEK1"]
plt$linear_fit(mapk ~ mek, subsets=)
