b = import('base')
io = import('io')
ar = import('array')
st = import('stats')
tcga = import('data/tcga')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% '../../scores/tcga/speed_matrix.RData'
OUTFILE = commandArgs(TRUE)[2] %or% 'pancan.pdf'
MUTFILE = "mutations_annotated_pathwayactivities_v3_mikeformat.txt"

tab = io$read_table(MUTFILE, header=TRUE) %>%
    transmute(hgnc_symbol = GENE_NAME,
              tcga_barcode = substr(Tumor_Sample_Barcode, 1, 16),
              study = tcga$barcode2study(Tumor_Sample_Barcode))

muts = unique(tab$hgnc_symbol)
scores = io$load(INFILE)
tissues = tcga$barcode2index(rownames(scores))$Study.Abbreviation
mut_matrix = matrix(FALSE, nrow=nrow(scores), ncol=length(muts),
                    dimnames=list(rownames(scores), muts))
for (i in 1:nrow(tab))
    if (tab$hgnc_symbol[i] %in% colnames(mut_matrix) &&
            tab$tcga_barcode[i] %in% rownames(scores))
        mut_matrix[tab$tcga_barcode[i],tab$hgnc_symbol[i]] = TRUE

# pan-cancer
pancan = st$lm(scores ~ tissues + mut_matrix) %>%
    filter(term == "mut_matrixTRUE") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"))

pdf(OUTFILE, width=26, height=20)

pancan %>%
    mutate(label = paste(mut_matrix, scores, sep=":")) %>%
    plt$color$p_effect(pvalue="adj.p") %>%
    plt$volcano(base.size=0.05, p=0.05)

pancan %>%
    plt$matrix(estimate ~ scores + mut_matrix) +
    xlab("Gene mutated") +
    ylab("Pathway-response genes")

dev.off()
