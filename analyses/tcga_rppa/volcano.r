b = import('base')
io = import('io')
plt = import('plot')

INFILE = commandArgs(TRUE)[1] %or% "assocs.RData"
OUTFILE = commandArgs(TRUE)[2] %or% "volcano.pdf"

assocs = io$load(INFILE)

pdf(OUTFILE, paper="a4r", width=26, height=20)

for (subs in unique(assocs$subset)) {
    message(subs)

    if (subs == "pan") {
        base_size = 0.05
    } else {
        base_size = 1
    }

    p = assocs %>%
        filter(subset == subs) %>%
        mutate(label = paste(scores, rppa, sep=":")) %>%
        plt$color$p_effect() %>%
        plt$volcano() +
            ggtitle(subs)

    print(p) %catch% NULL
}
