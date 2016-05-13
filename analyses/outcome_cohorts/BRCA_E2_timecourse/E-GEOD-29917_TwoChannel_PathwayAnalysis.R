##
# E-GEOD processing: Analysis from scratch
##

library(annotate)
library(oligo)
library(limma)
library(genefilter)
library(ggplot2)

# MCF-7:WS8 CELL LINES
setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_WS8")
agilent_files <-(list.files())
rawdata<-read.maimages(agilent_files,source = "agilent")
save(rawdata,file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_WS8/MCF-7_WS8_rawdata.ro")
load("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_WS8/MCF-7_WS8_rawdata.ro")
#Modifying sample labels
samples = c()
for(sample in colnames(rawdata)){
  samples = c(samples,strsplit(sample,"_")[[1]][1])
}
colnames(rawdata)<-samples

#Plotting Red-Green raw intensities
boxplot(log2(rawdata$R),main="Raw intensities RED")
boxplot(log2(rawdata$G),main="Raw intensities GREEN")

#Background correction
Bcorr <- backgroundCorrect(rawdata,method="normexp")
MA <- normalizeWithinArrays(Bcorr, method="none")
boxplot(MA$M,main="Raw intensities")

#Within Array Normalisation
MA <- normalizeWithinArrays(Bcorr, method="loess")
#Bad quality arrays
#1,4,7,9,14?
colnames(MA$M)[c(1,4,7,9)]
##
#[1]"GSM740721_Ag_251485015684_S01_GE2_107_Sep09_1_1" 
#[4]"GSM740724_Ag_251485015234_S01_GE2_107_Sep09_1_1" 
#[7]"GSM740727_Ag_251485015727_S01_GE2_107_Sep09_1_1"
#[9]"GSM740729_Ag_251485015232_S01_GE2_107_Sep09_1_3"
MA = MA[,c(-1,-4,-7,-9)]
boxplot(MA$M,main="Raw intensities")

#Filter expression object
# Filtering genes that are expressed
MA_data = MA[MA$genes$ControlType==0,]
MA_data = avereps(MA_data,ID=MA_data$genes[,"GeneName"])
expr_mat = MA_data[["M"]]

#load SPEED_matrix
load("/Users/admin/Documents/SPEED/root_files/model_matrix.RData")
useful_genes = intersect(rownames(zfit),rownames(expr_mat))
SPEED_scores_WS8 = t(expr_mat[useful_genes,]) %*% zfit[useful_genes,]
save(SPEED_scores_WS8,file="/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_WS8/SPEED_scores_WS8.ro")

#####################
# MCF-7:5C CELL LINES
#####################
setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_5C")
agilent_files <-(list.files())
rawdata<-read.maimages(agilent_files,source = "agilent",green.only = F)
save(rawdata,file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_5C/MCF-7_5C_rawdata.ro")
load("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_5C/MCF-7_5C_rawdata.ro")

samples = c()
for(sample in colnames(rawdata)){
  samples = c(samples,strsplit(sample,"_")[[1]][1])
}
colnames(rawdata)<-samples

#Plotting Red-Green raw intensities
boxplot(log2(rawdata$R),main="Raw intensities RED")
boxplot(log2(rawdata$G),main="Raw intensities GREEN")

#Background correction
Bcorr <- backgroundCorrect(rawdata,method="normexp")
MA <- normalizeWithinArrays(Bcorr, method="none")
boxplot(MA$M,main="Raw intensities")
#18,20,24 & 35
MA <- normalizeWithinArrays(Bcorr, method="loess")
boxplot(MA$M[,c(-18, -20,-24,-35)],main="Raw intensities")
MA = MA[,c(-18, -20,-24,-35)]

#Filter expression object
# Filtering genes that are expressed
MA_data = MA[MA$genes$ControlType==0,]
MA_data = avereps(MA_data,ID=MA_data$genes[,"GeneName"])
expr_mat = MA_data[["M"]]

#load SPEED_matrix

load("/Users/admin/Documents/SPEED/root_files/model_matrix.RData")
useful_genes = intersect(rownames(zfit),rownames(expr_mat))
SPEED_scores_5C = t(expr_mat[useful_genes,]) %*% zfit[useful_genes,]
save(SPEED_scores_5C, file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_5C/SPEED_scores_5C.ro")

#####################
# MCF-7:2A CELL LINES
#####################

setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_2A")
agilent_files <-(list.files())
rawdata<-read.maimages(agilent_files,source = "agilent")
save(rawdata,file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_2A/MCF-7_2A_rawdata.ro")
load("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_2A/MCF-7_2A_rawdata.ro")

samples = c()
for(sample in colnames(rawdata)){
  samples = c(samples,strsplit(sample,"_")[[1]][1])
}
colnames(rawdata)<-samples

#Plotting Red-Green raw intensities
boxplot(log2(rawdata$R),main="Raw intensities RED")
boxplot(log2(rawdata$G),main="Raw intensities GREEN")

#Background correction
Bcorr <- backgroundCorrect(rawdata,method="normexp")
MA <- normalizeWithinArrays(Bcorr, method="none")
boxplot(MA$M,main="Raw intensities")
#4,10,11
MA <- normalizeWithinArrays(Bcorr, method="loess")
boxplot(MA$M[,c(-4,-10,-11)],main="Raw intensities")

MA = MA[,c(-4,-10,-11)]

# Filter expression object
# Filtering genes that are expressed
MA_data = MA[MA$genes$ControlType==0,]
MA_data = avereps(MA_data,ID=MA_data$genes[,"GeneName"])
expr_mat = MA_data[["M"]]

#load SPEED_matrix

load("/Users/admin/Documents/SPEED/root_files/model_matrix.RData")
useful_genes = intersect(rownames(zfit),rownames(expr_mat))
SPEED_scores_2A = t(expr_mat[useful_genes,]) %*% zfit[useful_genes,]
save(SPEED_scores_2A, file="/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/MCF-7_2A/SPEED_scores_2A.ro")

#################
# Plotting scores
#################
library(ggplot2)
setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/")
load("MCF-7_2A/SPEED_scores_2A.ro")
load("MCF-7_5C/SPEED_scores_5C.ro")
load("MCF-7_WS8/SPEED_scores_WS8.ro")

SPEED_scores_all = rbind(SPEED_scores_2A,SPEED_scores_5C,SPEED_scores_WS8)
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

#Plotting
pdf(file = "/Users/admin/Documents/SPEED/E-GEOD_29917/Results/E-GEOD_29917_RicProc.pdf",width = 12,height = 10)
for(pathway in colnames(SPEED_scores_annotated)[1:11]){
  SPEED_plot <- ggplot(SPEED_scores_annotated, aes(x = TIME, y = SPEED_scores_annotated[,pathway], color = CELL, group = CLASS)) + 
    labs(y = pathway) + geom_point() + geom_line()
  print(SPEED_plot)
}
dev.off()