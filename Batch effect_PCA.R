setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

#-----------------{Batch Effect Correction"前"}-----------------
log2cpm_df<-read.csv("./RNA_DATA/log2cpm_matrix_before_BatchEffectCorrection.csv")
#dim(log2cpm_df) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_df[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose)

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


#-----------------{Batch Effect Correction"後"}-----------------
log2cpm_corrected_df<-read.csv("./RNA_DATA/log2cpm_matrix_after_BatchEffectCorrection.csv")
#dim(log2cpm_df) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_corrected_df[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose)

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
