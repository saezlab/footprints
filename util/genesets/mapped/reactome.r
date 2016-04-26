b = import('base')
ar = import('array')
er = import('../enrichr')

reactome = list(#"H2O2" = "",
#    "IL-1" = "Interleukin-1 signaling",
    "JAK-STAT" = c("Signaling by Interleukins",
                   "Interferon Signaling",
                   "Signalling to STAT3"),
    "MAPK" = "Signalling to ERKs",
    "EGFR" = c("Signaling by EGFR",
               "Signaling by EGFR in Cancer"),
    "PI3K" = c("PI3K Cascade",
               "Constitutive Signaling by Aberrant PI3K in Cancer",
               "PI3K/AKT Signaling in Cancer",
               "PI3K/AKT activation"),
    "TGFb" = "Signaling by TGF-beta Receptor Complex",
    "TNFa" = "TNF signaling",
    "VEGF" = "Signaling by VEGF",
#    "Wnt" = c("Signaling by Wnt",
#              "Signaling by WNT in cancer"),
#    "Insulin" = "Signaling by Insulin receptor",
#   "RAR" = "Signaling by Retinoic Acid",
    "p53" = "Transcriptional Regulation by TP53",
#    "Estrogen" = "",
    "Trail" = "TRAIL  signaling",
#    "notch" = "Signaling by NOTCH",
    "NFkB" = c("TAK1 activates NFkB by phosphorylation and activation of IKKs complex",
               "RIP-mediated NFkB activation via ZBP1"),
#   "PPAR" = c("PPARA activates gene expression",
#              "Activation of PPARGC1A (PGC-1alpha) by phosphorylation"),
#   "SHH" = "Signaling by Hedgehog"
    "Hypoxia" = "Cellular response to hypoxia"
)

INFILE = commandArgs(TRUE)[1] %or% "../ReactomePathways.gmt"
OUTFILE = commandArgs(TRUE)[2] %or% "reactome.RData"

lists = er$parse_gmt(INFILE)
 
# if mapping to SPEED pathways
pathways = sapply(names(reactome), function(path) {
    unique(unlist(lists[reactome[[path]]]))
})

save(pathways, file=OUTFILE)
