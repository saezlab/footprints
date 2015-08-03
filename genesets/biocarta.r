b = import('base')
ar = import('array')
er = import('./enrichr')

BioCarta = list("H2O2" = "ARENRF2_PATHWAY",
          "IL-1" = "IL1R_PATHWAY",
#          "IL1" = "IL1R_PATHWAY",
          "JAK-STAT" = "STAT3_PATHWAY",
#          "MAPK_only" = "MAPK_PATHWAY",
          "MAPK" = "MAPK_PATHWAY",
          "MAPK_PI3K" = "AKT_PATHWAY",
#          "both" = "AKT_PATHWAY",
          "TLR" = "TOLL_PATHWAY",
#          "PI3K_only" = "PTEN_PATHWAY",
          "PI3K" = "PTEN_PATHWAY",
#          "TGFB" = "TGFB_PATHWAY",
          "TGFb" = "TGFB_PATHWAY",
          "TNFa" = c("TNFR1_PATHWAY", "TNFR2_PATHWAY"),
          "VEGF" = "VEGF_PATHWAY",
          "Wnt" = "WNT_PATHWAY",
          "Insulin" = "IGF1_PATHWAY",
          "Hypoxia" = "P53HYPOXIA_PATHWAY",
          "RAR" = "RARRXR_PATHWAY",
          "p53" = "P53_PATHWAY",
          "Estrogen" = "CARM_ER_PATHWAY",
          "Trail" = "DEATH_PATHWAY",
          "notch" = "NOTCH_PATHWAY",
          "NFkB" = "NFKB_PATHWAY",
          "PPAR" = c("PPARA_PATHWAY", "PPARG_PATHWAY"))

INFILE = commandArgs(TRUE)[1] %or% "./Enrichr/src/main/resources/BioCarta.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "biocarta.RData"

lists = er$parse_gmt(INFILE)
 
# if mapping to SPEED pathways
pathways = sapply(names(BioCarta), function(path) {
    unique(unlist(lists[BioCarta[[path]]]))
})

# if using the original pathways
#pathways = [unique(unlist(reactome))]

save(pathways, file=OUTFILE)
