library(dplyr)
library(magrittr)
b = import('base')
import('./speed', attach=TRUE)

if (is.null(module_name())) {
    MODEL = commandArgs(TRUE)[1] %or% "../../model/model_matrix.r"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
    ZSCORES = commandArgs(TRUE)[3] %or% "../../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[4] %or% "speed_matrix_50.RData"
    FRAC = commandArgs(TRUE)[5] %or% 50

    # load zscores, model building function, and expression for each experiment
    zdata = io$load(ZSCORES)
    index = zdata$index
    zdata2model = import_(sub("\\.r$", "", MODEL))$zscore2model
    expr = io$load(EXPR)

    keep = index %>%
        group_by(pathway) %>%
        sample_frac(as.numeric(FRAC)/100) %$%
        id

    zdata$zscores = zdata$zscores[,keep]
    zdata$index = zdata$index[match(keep, zdata$index$id),]
    stopifnot(zdata$index$id == colnames(zdata$zscores))

    scores = clustermq::Q(expr2scores, id=index$id, job_size=1, memory=10240,
              const = list(expr = expr,
                           zdata = zdata,
                           zdata2model = zdata2model)) %>% #,
#                           hpc_args = list(n_jobs=20, memory=2048))) %>%
        setNames(zdata$index$id) %>%
        ar$stack(along=1) %>%
        ar$map(along=1, scale) #%>% # scale each pathway across all experiments
#        ar$map(along=2, scale) # scale each experiment across all pathways

    stopifnot(zdata$index$id == rownames(scores))

    save(scores, index, file=OUTFILE)
}
