# -> just take combinations of the different scores
# -> train a model, report some cross-val error
# -> see which {,2-way} combination is the best {max 4 features}
# --
# -> or also: full model, how much variability can scores/mut/etc explain

# Q: how much variance for each tissue-drug can be explained by tissue+<signature> features?

# drug response to 260 drugs (columns) <-> target
# pathway scores 19/50 pathways (columns) <-> features
# 19 different tissues (subset along rows)
import('base/operators')
io = import('io')
ar = import('array')
st = import('stats')
gdsc = import('data/gdsc')

args = commandArgs(TRUE)
OUTFILE = args[1] %or% 'model_vs_mut_NULL_1000rep_alldrugs.RData'

dset = list(
    speed = "../../scores/gdsc/speed_linear.RData",
    pathifier = "../../scores/merge/pathifier.RData",
    spia = "../../scores/merge/spia.RData",
    reactome = "../../scores/gdsc/gsea_reactome.RData",
    go = "../../scores/gdsc/gsea_go.RData"
) %>% io$load()

dset$tissues = gdsc$tissues(minN=10)
dset$mut = gdsc$mutated_genes(intogen=TRUE) + 0
dset$drugs = gdsc$drug_response('IC50')

dset = ar$intersect_list(dset, along=1)

tissues = dset$tissues
dset$tissues = NULL
drugs = dset$drugs
dset$drugs = NULL

#mut = list(mut=dset$mut)
#dset$mut = NULL

Qsquared = ...

result = st$ml(drugs ~ dset, subsets = tissues, xval=10, aggr=list(mlr::mse, mlr::mae, mlr::rmse),
               train_args = list("regr.glmnet"), atomic_class=NULL,
               hpc_args = list(n_jobs=200, memory=512)) # 200 jobs, 3 hrs per job

#@NEW mut + dset:
# 100 reps, 200 jobs: 80-85 min + 8 hours assembly

#@FIXED:
# for 1 rep, 20 jobs: 3% cpu @master / full @worker, 3 min runtime
# for 10 reps, 100 jobs: 10% cpu @master / full @worker, 14:22-14:32 (+ df processing time: 15:00)
# for 100 reps, 300 jobs: 80% cpu @master / ~30% @worker, 15:17-

save(result, file=OUTFILE)
