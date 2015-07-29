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
OUTFILE = args[1] %or% 'model_mse.RData'
#DATASETS = args[-1] %or% '../../scores/gdsc/speed_linear.RData'

speed = io$load("../../scores/gdsc/speed_linear.RData")
pathifier = io$load("../../scores/merge/pathifier.RData")
spia = io$load("../../scores/merge/spia.RData")

tissues = gdsc$tissues(minN=10)
drugs = gdsc$drug_response()

ar$intersect(tissues, drugs, speed, pathifier, spia, along=1)
dset = list(speed=speed, pathifier=pathifier, spia=spia)

drugs = drugs[,c(1,2)]
result = st$ml(drugs ~ dset, subsets = tissues, xval=10, aggr=TRUE,
               train_args = list("regr.glmnet"),
#               hpc_args = NULL)
               hpc_args = list(chunk.size=200, memory=512))

save(result, file=OUTFILE)
