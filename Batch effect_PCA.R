setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

#-----------------{PCA:Batch Effect Correction"前"}-----------------
log2cpm_df<-read.csv("./RNA_DATA/log2cpm_matrix_before_BatchEffectCorrection.csv")
#dim(log2cpm_df) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_df[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose,center=TRUE, scale.=FALSE) #只去除平均，不做zscore

print(pca_gene_sample)

pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(mat)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(mat[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(mat)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(mat[x,i]*pc2_weight[x])
  }
  pc2<-append(pc2,sample_pc2)
}

pc2_after<-as.numeric(unlist(pc2))

pca.var<-pca_gene_sample$sdev^2
pca.var.percentage<-round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.percentage,main='Scree plot',xlab='Principle Component',
        ylab='Percent Variation')

pca.data_2<-data.frame(Sample=rownames(gene_sample_transpose),PC1=(pc1_after),PC2=(pc2_after))
pca.data_2

pca.data_2$Group <- ifelse(
  grepl("^SRR", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  


#畫PC1和PC2圖
library(ggplot2)
ggplot(pca.data_2,aes(x=PC1,y=PC2,colour=Group))+
  scale_colour_manual(values=c('cyan3','tomato'))+
  #scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  #coord_cartesian(xlim = c(-0.5, 1.5)) +  # 限制 x 軸範圍到 0-2
  #geom_point(size=4, alpha=0.7) + #設定透明度
  geom_jitter(size = 4, width = 0.2, height = 0.2, alpha = 0.6) + #抖動 ->讓pc1&pc2相同的點可以分開
  xlab(paste('PC1 (',pca.var.percentage[1],'%)',sep=''))+
  ylab(paste('PC2 (',pca.var.percentage[2],'%)',sep=''))+
  theme_bw()+theme(axis.text.x=element_text(size=15),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   legend.title=element_blank(),
                   legend.position = c(0.9,0.9),
                   legend.background=element_rect(
                     colour='gray',fill="white"),
                   legend.key.width = unit(0.8, "cm"),
                   legend.key.height = unit(0.6, "cm"),
                   legend.text = element_text(size = 12)
  )


#-----------------{PCA:Batch Effect Correction"後"}-----------------
#{3371 control genes}
log2cpm_corrected_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection(1285Control Genes).csv")
#dim(log2cpm_df) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_corrected_df[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose,center=TRUE, scale.=FALSE)

print(pca_gene_sample)

pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(mat)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(mat[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(mat)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(mat[x,i]*pc2_weight[x])
  }
  pc2<-append(pc2,sample_pc2)
}

pc2_after<-as.numeric(unlist(pc2))

pca.var<-pca_gene_sample$sdev^2
pca.var.percentage<-round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.percentage,main='Scree plot',xlab='Principle Component',
        ylab='Percent Variation')

pca.data_2<-data.frame(Sample=rownames(gene_sample_transpose),PC1=(pc1_after),PC2=(pc2_after))
pca.data_2

pca.data_2$Group <- ifelse(
  grepl("^SRR", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  

X11()
#畫PC1和PC2圖
library(ggplot2)
ggplot(pca.data_2,aes(x=PC1,y=PC2,colour=Group))+
  scale_colour_manual(values=c('cyan3','tomato'))+
  #scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  #coord_cartesian(xlim = c(-0.5, 1.5)) +  # 限制 x 軸範圍到 0-2
  #geom_point(size=4, alpha=0.7) + #設定透明度
  geom_jitter(size = 4, width = 0.2, height = 0.2, alpha = 0.6) + #抖動 ->讓pc1&pc2相同的點可以分開
  xlab(paste('PC1 (',pca.var.percentage[1],'%)',sep=''))+
  ylab(paste('PC2 (',pca.var.percentage[2],'%)',sep=''))+
  theme_bw()+theme(axis.text.x=element_text(size=15),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   legend.title=element_blank(),
                   legend.position = c(0.9,0.9),
                   legend.background=element_rect(
                     colour='gray',fill="white"),
                   legend.key.width = unit(0.8, "cm"),
                   legend.key.height = unit(0.6, "cm"),
                   legend.text = element_text(size = 12)
  )

#-----------------{PCA:Batch Effect Correction"後"}-----------------
#{11control genes}
log2cpm_corrected_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection(11control_genes).csv")
#dim(log2cpm_df) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_corrected_df[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose,center=TRUE, scale.=FALSE)

print(pca_gene_sample)

pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(mat)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(mat[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(mat)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(mat[x,i]*pc2_weight[x])
  }
  pc2<-append(pc2,sample_pc2)
}

pc2_after<-as.numeric(unlist(pc2))

pca.var<-pca_gene_sample$sdev^2
pca.var.percentage<-round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.percentage,main='Scree plot',xlab='Principle Component',
        ylab='Percent Variation')

pca.data_2<-data.frame(Sample=rownames(gene_sample_transpose),PC1=(pc1_after),PC2=(pc2_after))
pca.data_2

pca.data_2$Group <- ifelse(
  grepl("^SRR", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  


#畫PC1和PC2圖
library(ggplot2)
ggplot(pca.data_2,aes(x=PC1,y=PC2,colour=Group))+
  scale_colour_manual(values=c('cyan3','tomato'))+
  #scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  #coord_cartesian(xlim = c(-0.5, 1.5)) +  # 限制 x 軸範圍到 0-2
  #geom_point(size=4, alpha=0.7) + #設定透明度
  geom_jitter(size = 4, width = 0.2, height = 0.2, alpha = 0.6) + #抖動 ->讓pc1&pc2相同的點可以分開
  xlab(paste('PC1 (',pca.var.percentage[1],'%)',sep=''))+
  ylab(paste('PC2 (',pca.var.percentage[2],'%)',sep=''))+
  theme_bw()+theme(axis.text.x=element_text(size=15),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   legend.title=element_blank(),
                   legend.position = c(0.9,0.9),
                   legend.background=element_rect(
                     colour='gray',fill="white"),
                   legend.key.width = unit(0.8, "cm"),
                   legend.key.height = unit(0.6, "cm"),
                   legend.text = element_text(size = 12)
  )


#-----------------{DEGs:Batch Effect Correction前後gene expression變化}-----------------
####{DEGsENSG00000077549}####
DEG_id<-"ENSG00000077549"
log2cpm_df<-read.csv("./RNA_DATA/log2cpm_matrix_before_BatchEffectCorrection.csv")
DEG_before_info<-log2cpm_df[log2cpm_df$X==DEG_id,]

#----轉成長df
library(tidyr)
library(dplyr) 
DEG_before_long <- DEG_before_info %>%
  pivot_longer(
    cols      = -X,           # 除了 X 之外的所有欄位都 pivot
    names_to  = "Sample",     # 原本的欄位名稱會變成 Sample 欄位
    values_to = "log2CPM"     # 原本的值會集中到 log2CPM 欄位
  )

DEG_before_long$Group <- ifelse(
  grepl("^SRR", DEG_before_long$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
DEG_before_df<-as.data.frame(DEG_before_long)

#----batch effect correction info
log2cpm_corrected_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection(1285Control Genes).csv")
DEG_after_info<-log2cpm_corrected_df[log2cpm_corrected_df$X==DEG_id,]

#----轉成長df
library(tidyr)
library(dplyr) 
DEG_after_long <- DEG_after_info %>%
  pivot_longer(
    cols      = -X,           # 除了 X 之外的所有欄位都 pivot
    names_to  = "Sample",     # 原本的欄位名稱會變成 Sample 欄位
    values_to = "log2CPM_corrected"     # 原本的值會集中到 log2CPM 欄位
  )

DEG_after_long$Group <- ifelse(
  grepl("^SRR", DEG_after_long$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
DEG_after_df<-as.data.frame(DEG_after_long)

#----merge 2個df
DEG_merged_df <- merge(
  DEG_before_df,
  DEG_after_df,
  by     = c("X", "Sample", "Group"),
)
DEG_merged_df$Group <- factor(DEG_merged_df$Group , levels = c("NMOSD", "Control")) #加入此>>>設定NMOSD在左邊
#校正前
library(ggplot2)

boxplot_log2cpm_before <- ggplot(DEG_merged_df, aes(x = Group, y =log2CPM, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =paste0("Expression of ", DEG_id, " before correction") ,
       x=NULL,y = "log2 CPM",fill="Group") +
  theme(plot.title = element_text(size = 16, hjust = 0.5,face='bold'),
        axis.text.x = element_text(size=13, angle=0,hjust = 0.5,vjust =1), #hjust=0.5置中對齊
        axis.title.y = element_text(size = 13),
        legend.background=element_rect(
          colour='white',fill="white"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))+
  scale_fill_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))

print(boxplot_log2cpm_before)
# 校正後

boxplot_log2cpm_after<- ggplot(DEG_merged_df, aes(x = Group, y = log2CPM_corrected, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =paste0("Expression of ", DEG_id, " after correction") ,
       x=NULL,y = "log2 CPM",fill="Group") +
  theme(plot.title = element_text(size = 16, hjust = 0.5,face='bold'),
        axis.text.x = element_text(size=13, angle=0,hjust = 0.5,vjust =1), #hjust=0.5置中對齊
        axis.title.y = element_text(size = 13),
        legend.background=element_rect(
          colour='white',fill="white"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))+
  scale_fill_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))

print(boxplot_log2cpm_after)


#-----------------------------------------------------------------------------------------
####{DEGs:ENSG00000180644}####
DEG_id<-"ENSG00000180644"
log2cpm_df<-read.csv("./RNA_DATA/log2cpm_matrix_before_BatchEffectCorrection.csv")
DEG_before_info<-log2cpm_df[log2cpm_df$X==DEG_id,]

#----轉成長df
library(tidyr)
library(dplyr) 
DEG_before_long <- DEG_before_info %>%
  pivot_longer(
    cols      = -X,           # 除了 X 之外的所有欄位都 pivot
    names_to  = "Sample",     # 原本的欄位名稱會變成 Sample 欄位
    values_to = "log2CPM"     # 原本的值會集中到 log2CPM 欄位
  )

DEG_before_long$Group <- ifelse(
  grepl("^SRR", DEG_before_long$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
DEG_before_df<-as.data.frame(DEG_before_long)

#----batch effect correction info
log2cpm_corrected_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection(1285Control Genes).csv")
DEG_after_info<-log2cpm_corrected_df[log2cpm_corrected_df$X==DEG_id,]

#----轉成長df
library(tidyr)
library(dplyr) 
DEG_after_long <- DEG_after_info %>%
  pivot_longer(
    cols      = -X,           # 除了 X 之外的所有欄位都 pivot
    names_to  = "Sample",     # 原本的欄位名稱會變成 Sample 欄位
    values_to = "log2CPM_corrected"     # 原本的值會集中到 log2CPM 欄位
  )

DEG_after_long$Group <- ifelse(
  grepl("^SRR", DEG_after_long$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
DEG_after_df<-as.data.frame(DEG_after_long)

#----merge 2個df
DEG_merged_df <- merge(
  DEG_before_df,
  DEG_after_df,
  by     = c("X", "Sample", "Group"),
)
DEG_merged_df$Group <- factor(DEG_merged_df$Group , levels = c("NMOSD", "Control")) #加入此>>>設定NMOSD在左邊
#校正前
library(ggplot2)

boxplot_log2cpm_before <- ggplot(DEG_merged_df, aes(x = Group, y =log2CPM, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =paste0("Expression of ", DEG_id, " before correction") ,
       x=NULL,y = "log2 CPM",fill="Group") +
  theme(plot.title = element_text(size = 16, hjust = 0.5,face='bold'),
        axis.text.x = element_text(size=13, angle=0,hjust = 0.5,vjust =1), #hjust=0.5置中對齊
        axis.title.y = element_text(size = 13),
        legend.background=element_rect(
          colour='white',fill="white"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))+
  scale_fill_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))

print(boxplot_log2cpm_before)
# 校正後

boxplot_log2cpm_after<- ggplot(DEG_merged_df, aes(x = Group, y = log2CPM_corrected, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =paste0("Expression of ", DEG_id, " after correction") ,
       x=NULL,y = "log2 CPM",fill="Group") +
  theme(plot.title = element_text(size = 16, hjust = 0.5,face='bold'),
        axis.text.x = element_text(size=13, angle=0,hjust = 0.5,vjust =1), #hjust=0.5置中對齊
        axis.title.y = element_text(size = 13),
        legend.background=element_rect(
          colour='white',fill="white"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))+
  scale_fill_manual(values = c("NMOSD" = "#E3A19F", "Control" = "#67D6F0"))

print(boxplot_log2cpm_after)

#-------------------------------------------------------------------------------------------


