b = import('base')
io = import('io')

speed1 = list(
    "JAK-STAT" = "JAK-STAT",
    "MAPK" = "MAPK_only",
    "EGFR" = "MAPK_PI3K",
    "PI3K" = "PI3K_only",
    "TGFb" = "TGFB",
    "TNFa" = "TNFa",
    "VEGF" = "VEGF"
)

INFILE = commandArgs(TRUE)[1] %or% "../speed1_web.txt"
OUTFILE = commandArgs(TRUE)[2] %or% "speed1.RData"

lists = io$read_table(INFILE, header=TRUE) %>%
    unstack(Gene_Name ~ Pathway)
 
# if mapping to SPEED pathways
pathways = sapply(names(speed1), function(path) {
    unique(unlist(lists[speed1[[path]]]))
})

save(pathways, file=OUTFILE)
