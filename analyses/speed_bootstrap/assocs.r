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
    mutate(mlogp = -log(p.value),
           abseff = abs(estimate),
           absstat = abs(statistic))

pdf("assocs.pdf", paper="a4r", width=20, height=20)

plt$matrix(result, mlogp ~ scores + pathway, palette="Blues") +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("- log(p.value)")

plt$matrix(result, abseff ~ scores + pathway, palette="Blues") +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("effect size")

#TODO: figure out what 'statistic' is exactly and how to use it
plt$matrix(result, absstat ~ scores + pathway, palette="Blues") +
    xlab("Pathway perturbed") +
    ylab("Assigned score") +
    ggtitle("statistic")

dev.off()
