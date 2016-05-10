##
# E-GEOD processing: single channel
##
library(annotate)
library(oligo)
library(limma)
library(genefilter)
library(ggplot2)

#Processing all arrays only in single channel
#setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all")
#agilent_files = (list.files())
#rawdata = read.maimages(agilent_files,source = "agilent",green.only = T)
#save(rawdata, file= "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all/rawexprdata.ro")

load("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all/rawexprdata.ro")

#Data processing
Bcorr = backgroundCorrect(rawdata,method="normexp")
exprdata = normalizeBetweenArrays(Bcorr,method = "quantile")
boxplot(exprdata$E)

#Data filtering
dim(exprdata)

#Filtering probes that are expressed - Limma suggestion
neg95 = apply(exprdata$E[exprdata$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
cutoff = matrix(1.1*neg95,nrow(exprdata),ncol(exprdata),byrow=TRUE)
isexpr = rowSums(exprdata$E > cutoff) >= 4
exprdata_filt = exprdata[exprdata$genes$ControlType==0 & isexpr,]
exprdata_gene = avereps(exprdata_filt,ID=exprdata_filt$genes[,"GeneName"])
expr_mat = exprdata_gene$E

#Modifying headers
samples = c()
for(sample in colnames(expr_mat)){
  samples = c(samples,strsplit(sample,"_")[[1]][1])
}
colnames(expr_mat)<-samples

#SPEED scores
load("/Users/admin/Documents/SPEED/root_files/model_matrix.RData")
useful_genes = intersect(rownames(zfit),rownames(expr_mat))
SPEED_scores_all = t(expr_mat[useful_genes,]) %*% zfit[useful_genes,]
save(SPEED_scores_all, file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all/SPEED_scores_all.ro")

#Plotting 
load("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all_targets.ro")
targets = all_samples
targets = targets[rownames(SPEED_scores_all),]
SPEED_scores_annotated = cbind(SPEED_scores_all,targets,CLASS = paste(targets$REPLICATE,targets$CELL,sep="."))

#Adjusting type of data
for(i in 1:11){
  SPEED_scores_annotated[,i] = as.numeric(SPEED_scores_annotated[,i])
}
SPEED_scores_annotated[,15] = as.numeric(SPEED_scores_annotated[,15])
SPEED_scores_annotated[,12] = as.factor(SPEED_scores_annotated[,12])
SPEED_scores_annotated[,13] = as.factor(SPEED_scores_annotated[,13])
SPEED_scores_annotated[,16] = as.factor(SPEED_scores_annotated[,16])

pdf(file = "/Users/admin/Documents/SPEED/E-GEOD_29917/Results/E-GEOD_29917_SingChan.pdf",width = 12,height = 10)
for(pathway in colnames(SPEED_scores_annotated)[1:11]){
  SPEED_plot <- ggplot(SPEED_scores_annotated, aes(x = TIME, y = SPEED_scores_annotated[,pathway], color = CELL, group = CLASS)) + 
    labs(y = pathway) + geom_point() + geom_line()
  print(SPEED_plot)
}
dev.off()

