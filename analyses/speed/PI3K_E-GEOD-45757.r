library(dplyr)
b = import('base')
io = import('io')
ar = import('array')
st = import('stats')

sensitive = c(
    "BxPC-3" = F,
    "Capan-2" = T,
    "CFPAC-1" = F,
    "COLO357" = NA,
    "HPAF_II" = F,
    "Hs766T" = T,
    "L3.3" = NA,
    "L3.6pl" = F,
    "L3.6sl" = NA,
    "MIAPaCa-2" = T,
    "Mpanc-96" = NA,
    "Panc-10.05" = T,
    "PANC-1" = T,
    "Panc-2.03" = T,
    "Panc-2.13" = F,
    "Panc-3.27" = F,
    "Panc-5.04" = T,
    "Panc-6.03" = F,
    "Panc-8.13" = F,
    "PL45" = NA,
    "SU86.86" = F,
    "SW1990" = F
)

INFILE = commandArgs(TRUE)[1] %or% "../../scores/speed/gsea_reactome.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "perturb2expr.pdf"

# load data
data = io$load(INFILE)
index = data$index %>%
    filter(accession == "E-GEOD-45757")
scores = data$scores[index$id,]

cells = data.frame(
    name = b$grep("^([A-Za-z0-9\\-\\_\\.]+)", index$cells),
    PI3K = scores[,'PI3K']
)
cells$sensitive = sensitive[match(cells$name, names(sensitive))]
#sign = ifelse(index$effect == "activating", 1, -1)
#pathway = sign * t(ar$mask(index$pathway)) + 0

# compute associations
result = st$lm(PI3K ~ sensitive, data=cells)
