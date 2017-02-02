b = import('base')
io = import('io')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "assocs.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "heatmap.pdf"

assocs = io$load(INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)

for (subs in unique(assocs$subset)) {
    message(subs)
    r2 = assocs %>%
        filter(subset == subs) %>%
        plt$cluster(estimate ~ scores + rppa) %>%
        group_by(rppa) %>%
        mutate(p.sum = sum(-log(adj.p))) %>%
        ungroup() %>%
        filter(adj.p < 0.05)

    keep = unique(r2[c('rppa','p.sum')]) %>%
        top_n(40, p.sum)

    r3 = r2 %>%
        filter(rppa %in% keep$rppa) %>%
        mutate(label = ifelse(adj.p < 0.01, ".", "")) %>%
        mutate(label = ifelse(adj.p < 0.001, "*", label)) %>%
        plt$matrix(estimate ~ scores + rppa, color="estimate",
                   label="label", reverse_colors=TRUE) +
            ggtitle(sprintf("%s (FDR .<0.01, *<1e-3)", subs))

    print(r3) %catch% NULL
}
