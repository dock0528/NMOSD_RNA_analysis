#install.packages('pheatmap')
#install.packages('magrittr')
#install.packages("gplots")
#install.packages("RColorBrewer")
#install.packages("devtools")
#install.packages("pheatmap")
library(devtools)
library(usethis)
library(gplots)
library(pheatmap)
library(readr)
library(magrittr)
library(dplyr)
library(tibble)

#################################【讀取 TPM matrix】#################################
#----{My}
my_protein.codong.tpm_df<-read.csv("./RNA_DATA/My_TPM_gene_samples_matrix(intersection gene)_v2.csv")
class(my_protein.codong.tpm_df) #"data.frame"

#----{External}
external_protein.codong.tpm_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene)_v2.csv")
class(external_protein.codong.tpm_df) #"data.frame"



################【7 Intersection IFN-I DEGs list (adjust.p<0.05 & |log2FC|<1)】#################
#----{8 Intersection protein coding-DEGs  list}
my_external_protein.coding.DEGs <- c(
  "ENSG00000187608",
  "ENSG00000126709",
  "ENSG00000134326",
  "ENSG00000131016",
  "ENSG00000160932",
  "ENSG00000135114",
  "ENSG00000159958"
)
library(clusterProfiler)
my_external_protein.coding.DEGs_names<-bitr(my_external_protein.coding.DEGs,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')

#######################【"My data" Pheatmap】#########################
#----{my}
my_protein.codong.tpm<-my_protein.codong.tpm_df[my_protein.codong.tpm_df$gene_id %in% my_external_protein.coding.DEGs,]
my_protein.codong.tpm <- merge(
  my_protein.codong.tpm,
  my_external_protein.coding.DEGs_names[, c("ENSEMBL", "SYMBOL")],
  by.x = "gene_id",  # 左邊用 gene_id
  by.y = "ENSEMBL",  # 右邊用 ENSEMBL
  all.x = TRUE       # 保留所有 my_protein.codong.tpm 的列
)
my_protein.codong.DEG.tpm_matrix<-as.matrix(my_protein.codong.tpm[ , -c(1, ncol(my_protein.codong.tpm))])#去除第一欄為ENSG名 &最後一欄Symbol
rownames(my_protein.codong.DEG.tpm_matrix) <- my_protein.codong.tpm$SYMBOL  #設定列名為gene id
#----{Sample & Group}
sample_df <- data.frame(
  Group = ifelse(
    grepl("^SRR", colnames(my_protein.codong.DEG.tpm_matrix)),
    "Control", "NMOSD"),
  row.names = colnames(my_protein.codong.DEG.tpm_matrix)
)

#----{樣本畫樹狀圖}
#install.packages('dendextend')
library(dendextend)
col<-colorRampPalette(c('blue','white','red'))(100)
euclidean_dist_row<-dist(my_protein.codong.DEG.tpm_matrix,method='euclidean')
# dist(t()) => distances among cols in Euclidean metric
complete_clusters_euclidean_row<-hclust(euclidean_dist_row,method='ward.D')

#----{基因排序}
mat_z <- t(scale(t(my_protein.codong.DEG.tpm_matrix))) #對gene做z-score
nmosd_idx    <- which(sample_df$Group == "NMOSD")
control_idx  <- which(sample_df$Group == "Control")
score <- rowMeans(mat_z[, nmosd_idx,   drop = FALSE]) -rowMeans(mat_z[, control_idx, drop = FALSE]) #Each gene mean(NMOSD)-mean(Control)
#order_genes <- rownames(mat_z)[order(score, decreasing = TRUE)] #Order gene 大->小

order_genes<-c('LY6E','IFI6','ISG15','OASL','CMPK2','TNFRSF13C','AKAP12') #Order gene 大->小
cluster_dend_rows_order <- as.hclust(rotate(as.dendrogram(complete_clusters_euclidean_row), order = order_genes))
#as.hclust 樹狀結構
#as.dendrogram 樹的物件(透過它才可以做排序轉換)

#----{pheatmap}設定畫布1000*1200
pheatmap(
  my_protein.codong.DEG.tpm_matrix,
  scale               = "row",
  color               = col,
  annotation_col      = sample_df,
  annotation_names_col= FALSE,
  annotation_colors   = list(
    Group = c(Control="#67D6F0", NMOSD="#E3A19F")
  ),
  cluster_rows        = cluster_dend_rows_order ,
  cluster_cols        = FALSE ,
  treeheight_row      = 0, #不畫dendrogram
  cellwidth           = 15,
  cellheight          = 35,
  show_rownames       = TRUE,
  show_colnames       = FALSE,
  legend_breaks       = c(-3, 3),
  legend_labels       = c("Low", "High"),
  legend_title        = NULL,
  main = "Hierarchical clustering of DEGs TPM (My dataset)",
  fontsize = 13,
  fontsize_row = 15
)

#scale='row'對於row的參數進行歸一化
#annotation_col註解
#col顏色
#annotation_col=dfh,
#fontsize_col col字體大小
#main 主題名
#treeheight_col =>col的樹狀圖高度(預設50)

#######################【"External data" Pheatmap】#########################
#----{external}
external_protein.codong.tpm<-external_protein.codong.tpm_df[external_protein.codong.tpm_df$gene_id %in% my_external_protein.coding.DEGs,]


# 篩選出欄位名稱含 "NA" 或 "no"
keep_cols <- grep("NA|no", colnames(external_protein.codong.tpm), value = TRUE)
external_protein.codong.tpm<- external_protein.codong.tpm[, c("gene_id", keep_cols)]

external_protein.codong.tpm <- merge(
  external_protein.codong.tpm,
  my_external_protein.coding.DEGs_names[, c("ENSEMBL", "SYMBOL")],
  by.x = "gene_id",  # 左邊用 gene_id
  by.y = "ENSEMBL",  # 右邊用 ENSEMBL
  all.x = TRUE       # 保留所有 my_protein.codong.tpm 的列
)
external_protein.codong.DEG.tpm_matrix<-as.matrix(external_protein.codong.tpm[ , -c(1, ncol(external_protein.codong.tpm))])#去除第一欄為ENSG名 &最後一欄Symbol
rownames(external_protein.codong.DEG.tpm_matrix) <- external_protein.codong.tpm$SYMBOL  #設定列名為gene id


#----{Sample & Group}
sample_df <- data.frame(
  Group = ifelse(
    grepl("^NA", colnames(external_protein.codong.DEG.tpm_matrix)),
    "Control", "NMOSD"),
  row.names = colnames(external_protein.codong.DEG.tpm_matrix)
)

#----{樣本畫樹狀圖}
#install.packages('dendextend')
library(dendextend)
col<-colorRampPalette(c('blue','white','red'))(100)
euclidean_dist_row<-dist(external_protein.codong.DEG.tpm_matrix,method='euclidean')
# dist(t()) => distances among cols in Euclidean metric
complete_clusters_euclidean_row<-hclust(euclidean_dist_row,method='ward.D')

#----{基因排序}
#mat_z <- t(scale(t(external_protein.codong.DEG.tpm_matrix))) #對gene做z-score
#nmosd_idx    <- which(sample_df$Group == "NMOSD")
#control_idx  <- which(sample_df$Group == "Control")
#score <- rowMeans(mat_z[, nmosd_idx,   drop = FALSE]) -
#rowMeans(mat_z[, control_idx, drop = FALSE]) #Each gene mean(NMOSD)-mean(Control)

#order_genes <- rownames(mat_z)[order(score, decreasing = TRUE)] #Order gene 大->小
order_genes<-c('LY6E','IFI6','ISG15','OASL','CMPK2','TNFRSF13C','AKAP12') #順序同my data
cluster_dend_rows_order <- as.hclust(rotate(as.dendrogram(complete_clusters_euclidean_row), order = order_genes))
#as.hclust 樹狀結構
#as.dendrogram 樹的物件(透過它才可以做排序轉換)

#----{pheatmap}設定畫布1000*1200
pheatmap(
  external_protein.codong.DEG.tpm_matrix,
  scale               = "row",
  color               = col,
  annotation_col      = sample_df,
  annotation_names_col= FALSE,
  annotation_colors   = list(
    Group = c(Control="#67D6F0", NMOSD="#E3A19F")
  ),
  cluster_rows        = cluster_dend_rows_order ,
  cluster_cols        = FALSE ,
  treeheight_row      = 0,
  cellwidth           = 15,
  cellheight          = 35,
  show_rownames       = TRUE,
  show_colnames       = FALSE,
  legend_breaks       = c(-3, 3),
  legend_labels       = c("Low", "High"),
  legend_title        = NULL,
  main = "Hierarchical clustering of DEGs TPM (External dataset)",
  fontsize = 13,
  fontsize_row = 15
)

