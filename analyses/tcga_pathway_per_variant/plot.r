library(dplyr)
plt = import('plot')

INFILE = commandArgs(TRUE)[2] %or% "assocs/speed_matrix.RData"
OUTFILE = commandArgs(TRUE)[3] %or% "plots/speed_matrix.pdf"
load(INFILE) # variants, cnas

pdf(OUTFILE)#, width=26, height=20)

variants %>%
    filter(!is.na(size)) %>%
    plt$volcano(base.size=10, p=0.1) + ggtitle("pan-cancer variants")

variants %>%
    filter(is.na(size)) %>%
    mutate(size = 50) %>%
    plt$volcano(p=0.1) + ggtitle("tissue-specific variants (interaction terms)")

cnas %>%
    filter(!is.na(size)) %>%
    plt$volcano(base.size=0.01, p=0.1) + ggtitle("pan-cancer CNAs")

cnas %>%
    filter(is.na(size)) %>%
    mutate(size = 50) %>%
    plt$volcano(p=0.1) + ggtitle("tissue-specific CNAs (interaction terms)")

dev.off()
