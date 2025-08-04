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

#################################【讀取log2cpm matrix】#################################
#----{Before}
log2cpm_before_df<-read.csv("./RNA_DATA/log2cpm_matrix_before_BatchEffectCorrection.csv")
class(log2cpm_before_df) #"data.frame"

#----{After}
log2cpm_after_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection(11Control Genes).csv")
class(log2cpm_after_df) #"data.frame"



################【1285 Control genes DEGs list (adjust.p<0.01 & |log2FC|<2)】#################
#----{1285 DEGs list}
degs_1285<-read.csv("./RNA_DATA/NMOSD_RNA_DEGS(11Control Genes Batch Effect Correction).csv") #355個DEGs

#----{Before}
degs_log2cpm_before<-log2cpm_before_df[log2cpm_before_df$X %in% degs_1285$x,]
degs_log2cpm_before_matrix<-as.matrix(degs_log2cpm_before[ , -1]) #去除第一欄為ENSG名

#----{After}
degs_log2cpm_after<-log2cpm_after_df[log2cpm_after_df$X %in% degs_1285$x,]
degs_log2cpm_after_matrix<-as.matrix(degs_log2cpm_after[ , -1]) #去除第一欄為ENSG名


#----{Sample & Group}
sample_df <- data.frame(
  Group = ifelse(
    grepl("^SRR", colnames(degs_log2cpm_before_matrix)),
    "Control", "NMOSD"),
  row.names = colnames(degs_log2cpm_before_matrix)
)

#######################【Pheatmap before correction 】#########################
#----{樣本畫樹狀圖}
#install.packages('dendextend')
library(dendextend)
col<-colorRampPalette(c('blue','white','red'))(100)
euclidean_dist_col<-dist(t(degs_log2cpm_before_matrix),method='euclidean')
# dist(t()) => distances among cols in Euclidean metric
complete_clusters_euclidean_col<-hclust(euclidean_dist_col,method='ward.D')


#----{修改dendrogram高}
hc2 <- complete_clusters_euclidean_col
hc2$height <- log(hc2$height)  

#----{pheatmap}設定畫布1000*1200
pheatmap(
  degs_log2cpm_before_matrix,
  scale               = "row",
  color               = col,
  annotation_col      = sample_df,
  annotation_names_col= FALSE,
  annotation_colors   = list(
    Group = c(Control="#67D6F0", NMOSD="#E3A19F")
  ),
  cluster_rows        = FALSE,
  cluster_cols        = hc2 ,
  treeheight_col      = 50,
  cellwidth           = 20,
  cellheight          = 0.15,
  show_rownames       = FALSE,
  show_colnames       = FALSE,
  legend_breaks       = c(-3, 3),
  legend_labels       = c("Low", "High"),
  legend_title        = NULL,
  main = "Hierarchical clustering of DEGs log2cpm before batch effect correction",
  fontsize = 14
)

#scale='row'對於row的參數進行歸一化
#annotation_col註解
#col顏色
#annotation_col=dfh,
#fontsize_col col字體大小
#main 主題名
#treeheight_col =>col的樹狀圖高度(預設50)

#######################【Pheatmap after correction 】#########################
#----{樣本畫樹狀圖}
#install.packages('dendextend')
library(dendextend)
col<-colorRampPalette(c('blue','white','red'))(100)
euclidean_dist_col<-dist(t(degs_log2cpm_after_matrix),method='euclidean')
# dist(t()) => distances among cols in Euclidean metric
complete_clusters_euclidean_col<-hclust(euclidean_dist_col,method='ward.D')


#----{修改dendrogram高}
hc2 <- complete_clusters_euclidean_col
hc2$height <- log(hc2$height)  

#----{pheatmap}
pheatmap(
  degs_log2cpm_after_matrix,
  scale               = "row",
  color               = col,
  annotation_col      = sample_df,
  annotation_names_col= FALSE,
  annotation_colors   = list(
    Group = c(Control="#67D6F0", NMOSD="#E3A19F")
  ),
  cluster_rows        = FALSE,
  cluster_cols        = hc2 ,
  treeheight_col      = 50,
  cellwidth           = 20,
  cellheight          = 0.15,
  show_rownames       = FALSE,
  show_colnames       = FALSE,
  legend_breaks       = c(-3, 3),
  legend_labels       = c("Low", "High"),
  legend_title        = NULL,
  main = "Hierarchical clustering of DEGs log2cpm after batch effect correction",
  fontsize = 14
)

