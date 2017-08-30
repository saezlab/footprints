downsampled_auc = function(n_sigs, expr, zdata, zdata2model) {
    library(dplyr)
    library(magrittr)
    index = zdata$index

    keep = zdata$index %>%
        group_by(pathway) %>%
        mutate(ni = sample(1:n(), size=n())) %>%
        filter(ni <= n_sigs) %$%
        id

    # build model
    zdata$zscores = zdata$zscores[,keep]
    zdata$index = zdata$index[match(keep, zdata$index$id),]
    stopifnot(zdata$index$id == colnames(zdata$zscores))
    model = zdata2model(zdata, hpc_args=list(n_jobs=0))$model

    # score experiments
    score_one = function(vecs, expr, index) {
        narray::intersect(vecs, expr, along=1)
        mat = t(expr) %*% vecs
        ctl = mat[index$control,,drop=FALSE]
        ptb = mat[index$perturbed,,drop=FALSE]
        colMeans(ptb) - colMeans(ctl)
    }
    scores = mapply(score_one, expr=expr$expr, index=expr$records,
                    MoreArgs = list(vecs=model)) %>%
        t() %>%
        narray::map(along=1, scale)

    # calculate ROC AUC
    st = import('stats')
    roc = import('./roc_util')
    auc = list(scores=scores, index=index[match(rownames(scores), index$id),]) %>%
        roc$scores2df() %>%
        na.omit() %>%
        mutate(method = "speed_matrix",
               inferred = signature,
               matched = perturbed == inferred) %>%
        group_by(inferred, signature) %>%
        do(st$roc(., "score", "matched")) %>%
        ungroup() %>%
        roc$roc2auc() %>%
        transmute(pathway=inferred, auc=PROGENy)
}

if (is.null(module_name())) {
    library(dplyr)
    library(magrittr)
    b = import('base')
    io = import('io')

    MODEL = commandArgs(TRUE)[1] %or% "../../model/model_matrix.r"
    EXPR = commandArgs(TRUE)[2] %or% "../../data/expr.RData"
    ZSCORES = commandArgs(TRUE)[3] %or% "../../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[4] %or% "speed_downsample_auc.RData"

    # load zscores, model building function, and expression for each experiment
    zdata = io$load(ZSCORES)
    zdata2model = import_(sub("\\.r$", "", MODEL))$zscore2model
    expr = io$load(EXPR)

    max_sigs = zdata$index %>%
        group_by(pathway) %>%
        summarize(n=n())

    n_sigs = rep(3:max(max_sigs$n), 3)
    aucs = clustermq::Q(downsampled_auc, n_sigs=n_sigs, job_size=1,
                        const = list(expr=expr, zdata=zdata, zdata2model=zdata2model))

    aucs = data_frame(n_sigs=n_sigs, auc=aucs) %>%
        tidyr::unnest() %>%
        left_join(max_sigs, by="pathway") %>%
        filter(n_sigs <= n) %>%
        select(-n)

    save(aucs, file=OUTFILE)
}
