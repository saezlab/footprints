library(modules)
gn = import('general')
plt = import('plots')
er = import('GSEA-Lists/enrichr')
bm = import('biomart')
sg = import('sanger_robject')
gs = import('GSEA')
ll = import('list')
an = import('anova')
ar = import('array')

# get gene lists for SPEED, GO, and pathway

GO = list("H2O2" = "GO:0042542", # response to hydrogen peroxide
          "IL-1" = "GO:0070498", # interleukin-1-mediated signalling pathway
          "IL1" = "GO:0070498",
          "JAK-STAT" = "GO:0007259", # JAK-STAT signal transduction
          "MAPK_only" = "GO:0000165", # MAPK cascade
          "MAPK" = "GO:0000165",
          "MAPK_PI3K" = "GO:0038127", # ERBB signaling pathway
          "both" = "GO:0038127",
          "TLR" = "GO:0002224", # toll-like receptor signaling pathway
          "PI3K_only" = "GO:0014065", # phosphatidylinositol 3-kinase signaling
          "PI3K" = "GO:0014065",
          "TGFB" = "GO:0007179", # transforming growth factor beta receptor signaling pathway
          "TGFb" = "GO:0007179",
          "TNFa" = "GO:0033209", # tumor necrosis factor-mediated signaling pathway
          "VEGF" = "GO:0038084", # vascular endothelial growth factor signaling pathway
          "Wnt" = "GO:0016055", # Wnt signaling pathway
          "Insulin" = "GO:0008286", # insulin receptor signaling pathway
          "Hypoxia" = "GO:0071456", # cellular response to hypoxia
          "RAR" = "GO:0048384", # retinoic acid receptor signaling pathway
          "p53" = "GO:0072331", # signal transduction by p53 class mediator
          "Estrogen" = "GO:0030520", # intracellular estrogen receptor signaling pathway
          "Trail" = "GO:0036462", # TRAIL-activated apoptotic signaling pathway
          "notch" = "GO:0007219", # Notch signaling pathway
          "NFkB" = "GO:0038061", # NIK/NF-kappaB signaling
          "PPAR" = "GO:0035357") # peroxisome proliferator activated receptor signaling pathway

BioCarta = list("H2O2" = "ARENRF2_PATHWAY",
          "IL-1" = "IL1R_PATHWAY",
          "IL1" = "IL1R_PATHWAY",
          "JAK-STAT" = "STAT3_PATHWAY",
          "MAPK_only" = "MAPK_PATHWAY",
          "MAPK" = "MAPK_PATHWAY",
          "MAPK_PI3K" = "AKT_PATHWAY",
          "both" = "AKT_PATHWAY",
          "TLR" = "TOLL_PATHWAY",
          "PI3K_only" = "PTEN_PATHWAY",
          "PI3K" = "PTEN_PATHWAY",
          "TGFB" = "TGFB_PATHWAY",
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

# load speed pathways and genes
#speed = gn$data('GSEA-Lists/speed')$speed1
speed = gn$data('TAC-y2/speed_lists_20140626')$opt50

go = bm$getHGNCByGO(unlist(GO[names(speed)]))

biocIndex = er$parseEnrichrFiles("BioCarta")$BioCarta
bc = sapply(BioCarta[names(speed)], function(x) biocIndex[x])
bc = lapply(bc, function(x) unique(unlist(x, use.names=F)))

genes = ll$transpose(list(SPEED=speed, GO=go, BioCarta=bc))

# plot the assocs here as well
expr = sg$getBASAL_EXPRESSION()
scores = lapply(ll$flatten(genes), function(l) gs$runGSEA(expr, l))
mat = ar$stack(scores)

Ys = sg$getDrugResponseForCellLines('IC50s') # or AUC
tissues = sg$getTissues(minN=5)
ar$intersect(mat, tissues, Ys, along=1)
Yf = sg$filterDrugResponse(Ys, tissues, top=0.1, abs=0, delta=2)
assocs.pan = an$calcAssocs(Ys, mat, covariate=tissues, p.adjust="fdr")
assocs.tissue = an$calcAssocs(Yf, mat, subsets=tissues, p.adjust="fdr", stack=T)

p.pan = plt$drawVolcano(assocs.pan, top=40, log='y', base.size=0.2) +
    ggtitle("Pan-cancer pathway associations, 5% FDR") + 
    xlab("Regression slope") + 
    ylab("FDR-adjusted p-value (log)")

p.tis = plt$drawVolcano(assocs.tissue, top=30, log='y', base.size=2, p=0.2) +
    ggtitle("Tissue-specific pathway associations, 20% FDR") + 
    xlab("Regression slope") + 
    ylab("FDR-adjusted p-value (log)")

# draw venn diagrams for each triple
pdf("VennPlots.pdf")
lapply(genes, function(x) {
    grid.newpage()
    plt$drawVenn(x)
})

#tsp = ar$map(mat[,grepl("\\.SPEED", colnames(mat))], along=1, mean, subsets=tissues)
#pheatmap(t(tsp), scale='none')
#tgo = ar$map(mat[,grepl("\\.GO", colnames(mat))], along=1, mean, subsets=tissues)
#pheatmap(t(tgo), scale='none')
#tbc = ar$map(mat[,grepl("\\.BioCarta", colnames(mat))], along=1, mean, subsets=tissues)
#pheatmap(t(tbc), scale='none')

print(p.pan)
print(p.tis)
dev.off()

#TODO: plot sp/go/bioc w+w/o names, covariate and subset

