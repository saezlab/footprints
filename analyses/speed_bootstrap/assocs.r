library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')
plt = import('plot')

data = io$load("../../scores/speed/speed_matrix.RData")
index = data$index
scores = data$scores

sign = ifelse(index$effect == "activating", 1, -1)
pathway = sign * t(ar$mask(index$pathway)) + 0
#pathway["EGFR",] = pathway["EGFR",] + pathway["MAPK",] + pathway["PI3K",]
#pathway["TNFa",] = pathway["TNFa",] + pathway["NFkB",]

result = st$lm(scores ~ pathway) %>%
    mutate(mlogp = -log(p.value)) %>%
    mutate(label = ifelse(p.value < 1e-5, "*", "")) %>%
    mutate(label = ifelse(p.value < 1e-10, "***", label))

pdf("assocs.pdf", width=8, height=6)

plt$matrix(result, mlogp ~ scores + pathway, palette="RdBu", limits=c(-280,280)) +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("- log(p.value)") +
    theme(axis.text.x = element_text(color="black", size=14),
          axis.text.y = element_text(color="black", size=14))

result$label = NULL

plt$matrix(result, estimate ~ scores + pathway, palette="RdBu", limits=c(-4,4)) +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("Effect size") +
    theme(axis.text.x = element_text(color="black", size=14),
          axis.text.y = element_text(color="black", size=14))

#TODO: figure out what 'statistic' is exactly and how to use it
plt$matrix(result, statistic ~ scores + pathway, palette="RdBu") +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("statistic")

dev.off()
