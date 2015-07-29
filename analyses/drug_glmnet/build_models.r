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
OUTFILE = args[1] %or% 'speed_linear.RData'
DATASETS = args[-1] %or% '../../scores/gdsc/speed_linear.RData'

tissues = gdsc$tissues(minN=10)
drugs = gdsc$drug_response()
dset = ar$stack(lapply(DATASETS, io$load), along=2)
ar$intersect(tissues, drugs, dset, along=1)

df = import('data_frame')
result = st$ml(drugs ~ dset, subsets = tissues, atomic = "dset", aggr=TRUE,
               train_args = list("regr.glmnet"),
               hpc_args = NULL)
#               hpc_args = list(chunk.size=200, memory=512))

save(result, file=OUTFILE)
