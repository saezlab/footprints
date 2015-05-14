b = import('base')
ar = import('array')
er = import('./enrichr')

reactome = list(#"H2O2" = "",
#          "IL-1" = "",
#          "JAK-STAT" = "",
#          "MAPK" = "",
          "MAPK_PI3K" = "Signaling by EGFR",
#          "TLR" = "",
#          "PI3K" = "",
          "TGFb" = "Signaling by TGF beta",
#          "TNFa" = "",
          "VEGF" = "Signaling by VEGF",
          "Wnt" = "Signaling by Wnt",
          "Insulin" = "Signaling by Insulin receptor"
#          "Hypoxia" = "",
#          "RAR" = "",
#          "p53" = "",
#          "Estrogen" = "",
#          "Trail" = "",
#          "notch" = "",
#          "NFkB" = "",
#          "PPAR" = ""
)

INFILE = commandArgs(TRUE)[1] %or% "./Enrichr/src/main/resources/Reactome.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "reactome.RData"

lists = er$parse_gmt(INFILE)[unique(unlist(reactome))]
save(lists, file=OUTFILE)
