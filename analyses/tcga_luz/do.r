io = import('io')
ar = import('array')
st = import('stats')
tcga = import('data/tcga')
plt = import('plot')

muts = c("KRAS", "BRAF", "HRAS", "EGFR", "FGFR3", "TP53", "PIK3CA", "PTEN", "CTNNB1", "VHL", "NOTCH1", "AKT1", "MAP2K1", "MAP3K1", "MTOR", "PIK3R1", "SMAD4")
# CASP8, SMARC4 not mutated where we have scores

tab = io$read_table("mutations_annotated_pathwayactivities_v3_mikeformat.txt", header=TRUE) %>%
    transmute(hgnc_symbol = GENE_NAME,
              tcga_barcode = substr(Tumor_Sample_Barcode, 1, 16),
              code = tcga$barcode2index(Tumor_Sample_Barcode)$Study.Abbreviation)

scores = io$load('../../scores/tcga/speed_matrix.RData')
tissues = tcga$barcode2index(rownames(scores))$Study.Abbreviation
mut_matrix = matrix(FALSE, nrow=nrow(scores), ncol=length(muts),
                    dimnames=list(rownames(scores), muts))
for (i in 1:nrow(tab))
    if (tab$hgnc_symbol[i] %in% colnames(mut_matrix) && tab$tcga_barcode[i] %in% rownames(scores))
        mut_matrix[tab$tcga_barcode[i],tab$hgnc_symbol[i]] = TRUE

# pan-cancer
pancan = st$lm(scores ~ tissues + mut_matrix) %>%
    filter(term == "mut_matrixTRUE") %>%
    select(-term, -size) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    filter(adj.p < 0.05)

pdf("do.pdf", width=10, height=10)

pancan %>%
    plt$matrix(estimate ~ scores + mut_matrix) +
    xlab("Gene mutated") +
    ylab("Pathway-response genes")

dev.off()

# tissue specific
do_plot = function(scores, mut_matrix, tissue) {
    scores = scores[tissues == tissue,]
    mut_matrix = mut_matrix[tissues == tissue,]
    mut_matrix = mut_matrix[,colSums(mut_matrix) != 0] # have at least one mutation

    tissue = st$lm(scores ~ mut_matrix) %>%
        filter(term == "mut_matrixTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        filter(adj.p < 0.1)

    tissue %>%
        plt$matrix(estimate ~ scores + mut_matrix) +
        xlab("Gene mutated") +
        ylab("Pathway-response genes")
}
# there are none for BRCA, READ, OV, CESC; only 1 for COAD????
ar$map(mut_matrix, along=1, function(x) sum(x), subsets=tissues)
