#
# E-GEOD-29917 processed data analysis
#
library(biomaRt)
library(limma)
library(ggplot2)

# Generate raw expression matrix
setwd("/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/E-GEOD-29917.processed")
file_list = list.files()
raw_expmat = c()
for(file in file_list){
  print(file)
  sample_name = strsplit(file,"_")[[1]][1]
  processed_table = read.table(file, sep="\t",header = T,stringsAsFactors = F)
  Rnames = processed_table[,1]
  processed_table = processed_table[,2,drop=F]
  processed_table[,1] = as.numeric(processed_table[,1])
  processed_table = as.matrix(processed_table)
  rownames(processed_table) = Rnames
  #processed_table = processed_table[,2,drop=F]
  colnames(processed_table) = sample_name
  #processed_table[,1] = as.numeric(processed_table[,1])
  raw_expmat = cbind(raw_expmat,processed_table)
}

# repeated probes
rp = unique(row.names(raw_expmat)[duplicated(row.names(raw_expmat))])

# Generate annotation db
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
IDs = unique(rownames(raw_expmat))
IDs_hg = getBM(attributes=c('efg_agilent_wholegenome_4x44k_v1', 'hgnc_symbol'), 
               filters = 'efg_agilent_wholegenome_4x44k_v1', values = IDs, mart = ensembl)
IDs_hg = IDs_hg[IDs_hg[,2] != "",]
#save(IDs_hg,file = "/Users/admin/Documents/SPEED/root_files/agilent_hg.ro")
#load("/Users/admin/Documents/SPEED/root_files/agilent_hg.ro")

# 1st expr_mat filter
raw_expmat_filt = raw_expmat[rownames(raw_expmat) %in% IDs_hg[,1],]

#Transform expr matrix to an hgnc-exprmat
#Super inefficient 
expr_mat = c()
Rnames = c()
for(i in 1:nrow(raw_expmat_filt)){
  print(i)
  ID = rownames(raw_expmat_filt)[i]
  if(ID %in% IDs_hg[,1]){
    ID_df = IDs_hg[IDs_hg[,1] == ID,,drop=F]
    symbols = ID_df[,2]
    for(gene in symbols){
      expr_mat = rbind(expr_mat,raw_expmat_filt[i,])
      Rnames = c(Rnames,gene)
    }
  }
}

rownames(expr_mat) = Rnames
save(expr_mat,file="/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/Rawexpr_mat_hgnc.ro")

#Average similar probes
expr_mat = avereps(expr_mat,ID=rownames(expr_mat))

#SPEED scores
load("/Users/admin/Documents/SPEED/root_files/model_matrix.RData")
useful_genes = intersect(rownames(zfit),rownames(expr_mat))
SPEED_scores_all = t(expr_mat[useful_genes,]) %*% zfit[useful_genes,]
save(SPEED_scores_all, file = "/Users/admin/Documents/SPEED/E-GEOD_29917/root_files/all/SPEED_scores_all_processed.ro")

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

pdf(file = "/Users/admin/Documents/SPEED/E-GEOD_29917/Results/E-GEOD_29917_Processed.pdf",width = 12,height = 10)
for(pathway in colnames(SPEED_scores_annotated)[1:11]){
  SPEED_plot <- ggplot(SPEED_scores_annotated, aes(x = TIME, y = SPEED_scores_annotated[,pathway], color = CELL, group = CLASS)) + 
    labs(y = pathway) + geom_point() + geom_line()
  print(SPEED_plot)
}
dev.off()
