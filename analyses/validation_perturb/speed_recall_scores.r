library(dplyr)
io = import('io')
ar = import('array')
st = import('stats')

# load all in ../scores/speed/*.RData
objs = io$load_regex("(.*)\\.RData", "../scores/speed")
index = objs$speed_linear$index # any full index would do, split this later so no duplicate data
scores = lapply(objs, function(o) {
    # scale the pathway scores for each sample to mean=0, sd=1
    re = ar$map(o$scores, along=2, scale)
    rownames(re) = o$index$id
    re
}) %>% ar$intersect_list(along=1)

# subset only to those where all scores are available
# for each pathway, do positive control for (+/-) perturbation

# MAPK activation: include MAPK, EGFR activation
perturb2assocs = function(path="MAPK", perturb=c("EGFR", "MAPK"), p_effect="activating") {
    sub_idx = index %>%
        filter(pathway %in% perturb & effect == p_effect) %>%
        filter(id %in% rownames(scores[[1]]))
    message(path, " : ", nrow(sub_idx))
    score = lapply(scores, function(x) x[,path] %catch% NULL) %>%
        ar$stack(along=2)
    positive = rownames(score) %in% sub_idx$id + 0
    st$lm(score ~ positive) %>%
        arrange(p.value) %>%
        mutate(pathway = path) %>%
        transmute(pathway, sign=p_effect, method=score, estimate, p.value, size=sum(positive))
}

checks = list(
    MAPKp = list(path="MAPK", perturb=c("EGFR", "MAPK"), p_effect="activating"),
    MAPKm = list(path="MAPK", perturb="MAPK", p_effect="inhibiting"),
    EGFRp = list(path="EGFR", perturb="EGFR", p_effect="activating"),
    EGFRm = list(path="EGFR", perturb="EGFR", p_effect="inhibiting"), # ONLY 6 CONTROL
    PI3Kp = list(path="PI3K", perturb="PI3K", p_effect="activating"), # ONLY 2
    PI3Km = list(path="PI3K", perturb="PI3K", p_effect="inhibiting"),
    TNFap = list(path="TNFa", perturb="TNFa", p_effect="activating"),
#    TNFam = list(path="TNFa", perturb="TNFa", p_effect="inhibiting"),
    NFkBp = list(path="NFkB", perturb=c("TNFa", "NFkB"), p_effect="activating"),
#    NFkBm = list(path="NFkB", perturb="TNFa", p_effect="inhibiting"),
    p53p = list(path="p53", perturb="p53", p_effect="activating"), # ONLY 5
#    p53m = list(path="p53", perturb="p53", p_effect="inhibiting")
    JAKSTATp = list(path="JAK-STAT", perturb="JAK-STAT", p_effect="activating"),
#    JAKSTATm = list(path="JAK-STAT", perturb="JAK-STAT", p_effect="inhibiting"),
#    TGFbp = list(path="TGFb", perturb="TGFb", p_effect="activating"),
#    TGFbm = list(path="TGFb", perturb="TGFb", p_effect="inhibiting"),
    VEGFp = list(path="VEGF", perturb="VEGF", p_effect="activating"),
    VEGFm = list(path="VEGF", perturb="VEGF", p_effect="inhibiting"), # ONLY 1
#    Wntp = list(path="Wnt", perturb="Wnt", p_effect="activating"),
#    Wntm = list(path="Wnt", perturb="Wnt", p_effect="inhibiting"),
    Hypoxiap = list(path="Hypoxia", perturb="Hypoxia", p_effect="activating"),
#    Hypoxiam = list(path="Hypoxia", perturb="Hypoxia", p_effect="inhibiting"),
    Trailp = list(path="Trail", perturb="Trail", p_effect="activating") # ONLY 3
#    Trailm = list(path="Trail", perturb="Trail", p_effect="inhibiting")
)
lapply(checks, function(c) do.call(perturb2assocs, c))
