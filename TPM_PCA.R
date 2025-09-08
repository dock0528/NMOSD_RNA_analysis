setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

TPM_df<-read.csv("./RNA_DATA/TPM_gene_samples_matrix.csv")
#dim(TPM_df) gene x sample: 78724 x 22('gene_id+samples)

#-----------【Z-score】-------------
#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(TPM_df[ , -1]) #-1 >>> 把第 1 欄拿掉

#scale()預設是對 column 做標準化
z_mat <- t( scale(t(mat), center = TRUE, scale = TRUE) )
z_mat[ is.na(z_mat) ] <- 0 #把gene在所有samples TPM都為0的設NaN

TPM_gene_normalization<- data.frame(
  gene_id = TPM_df$gene_id,
  z_mat,
  check.names = FALSE #不改變欄位名稱
) #dim(TPM_gene_normalization) 78724 x 22 >>> gene x sample

head(TPM_gene_normalization)

all_gene_sample <- as.matrix(TPM_gene_normalization[,-1])# gene x sample
print(all_gene_sample) #列為gene 欄為sample

gene_sample_transpose<-t(all_gene_sample) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose)

print(pca_gene_sample)

pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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


###############################【篩"My data"的 56657 genes PCA】###############################
my_tpm_df<-read.csv("./RNA_DATA/My_TPM_gene_samples_matrix(intersection gene).csv")
#dim(my_tpm_df) 56657*22

all_gene_sample <- as.matrix(my_tpm_df[, -1])  # gene x sample #-1 >>> 把第 1 欄拿掉



rownames(all_gene_sample) <- my_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重



#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  ggtitle("My dataset PCA ")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
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

###############【篩"External data"的 56657 genes PCA】###############
external_tpm_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene).csv")
#dim(external_tpm_df) 56657*45

all_gene_sample <- as.matrix(external_tpm_df[, -1]) # gene x sample


rownames(all_gene_sample) <- external_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重
#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "NA" 開頭
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
  ggtitle("External dataset PCA ")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
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

#----{Protein coding genes list}
inter_protein_coding_genes<-read.csv("./RNA_DATA/External_with_My_protein_coding_genes.csv") #19722 intersect protein coding genes

###############【篩"My data"的 19722 Protein coding genes PCA】###############
my_tpm_df<-read.csv("./RNA_DATA/My_TPM_gene_samples_matrix(intersection gene).csv")
my_prtein_coding_tpm_df<-my_tpm_df[my_tpm_df$gene_id %in% inter_protein_coding_genes$x,]
#dim(my_prtein_coding_tpm_df) 19722 x 22

all_gene_sample <- as.matrix(my_prtein_coding_tpm_df[, -1]) # gene x sample


rownames(all_gene_sample) <- my_prtein_coding_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#######################【pc1 weighting 找elbow】#######################
#---{pc1 weighting排序}
pc1_sorted <- sort(abs(pc1_weight), decreasing = TRUE)

# ----{use 19722 genes find elbow}
topN <- length(pc1_sorted)
y    <- pc1_sorted[1:topN] #abs(pc1_weighting)
x    <- 1:topN

# ----{每個點正規化到 [0,1]}
x_norm <- (x - min(x)) / (max(x) - min(x))
y_norm <- (y - min(y)) / (max(y) - min(y))

# ----{計算點到斜直線距離}
x1 <- 0; y1 <- 1 #左上點
x2 <- 1; y2 <- 0 #右下點
num   <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + (x2 * y1 - y2 * x1))  #分母
denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2) #分子
dist_to_line <- num / denom

# ----{elbow information}
elbow_rank  <- which.max(dist_to_line)
elbow_abs.pc1_value <- y[elbow_rank]
elbow_genes <- names(pc1_sorted)[1:elbow_rank]

cat("Elbow rank:", elbow_rank, "\n")
cat("Elbow threshold (|PC1 weighting|):", elbow_abs.pc1_value, "\n")
cat("Top genes criteria:", length(elbow_genes), "\n")

#----{Visualization}
plot.point.num <- max(1, floor(topN / 2000))  # 最多畫 ~2000 點
idx_plot   <- seq(1, topN, by = plot.point.num )

plot(idx_plot, y[idx_plot], type = "b", pch = 19, cex = 0.6,, col = "#446270",
     xlab = "Ranked genes",
     ylab = "| PC1 weighting |",
     main = sprintf("Elbow plot PC1 weighting (%d genes)", topN),
     cex.main = 1.3) #cex.main標題大小

# 標註 elbow 位置
abline(v = elbow_rank, lty = 3)
abline(h = elbow_abs.pc1_value, lty = 3)
points(elbow_rank, elbow_abs.pc1_value, pch = 21, bg = "#EF777F", col = "#446270", cex = 2.1)
text(elbow_rank, elbow_abs.pc1_value*3,
     labels = paste0("Elbow rank:", elbow_rank),
     pos = 4, cex = 1, offset = 0.5) #pos=4文字放右邊

###############【看pc1 前100個權重較高的genes富集到的pathway】###############
#取前100個PC1 weighting較高的genes
top300_genes <- names(sort(abs(pc1_weight), decreasing = TRUE))[1:238]

#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  ggtitle("PCA filter protein coding genes (My dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
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

###############【篩"External data"的 19722 Protein coding genes PCA】###############
external_tpm_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene).csv")
external_prtein_coding_tpm_df<-external_tpm_df[external_tpm_df$gene_id %in% inter_protein_coding_genes$x,]
#dim(external_prtein_coding_tpm_df) 19722 x 45


all_gene_sample <- as.matrix(external_prtein_coding_tpm_df[, -1]) # gene x sample



rownames(all_gene_sample) <-external_prtein_coding_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重


#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  

#依PC1做排序
pca.data_2[order(pca.data_2$PC1, decreasing = TRUE), ]

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
  ggtitle("PCA filter protein coding genes (External dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
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



#-----------------------【"My data"的 19722 Protein coding genes PCA】-------------------------------
my_prtein_coding_tpm_df<-read.csv("./RNA_DATA/My_TPM_gene_samples_matrix(intersection gene)_v2.csv")
#dim(my_prtein_coding_tpm_df) 19722 x 22

all_gene_sample <- as.matrix(my_prtein_coding_tpm_df[, -1]) # gene x sample

rownames(all_gene_sample) <- my_prtein_coding_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重

#######################【pc1 weighting 找elbow】#######################
#---{pc1 weighting排序}
pc1_sorted <- sort(abs(pc1_weight), decreasing = TRUE)

# ----{use 19722 genes find elbow}
topN <- length(pc1_sorted)
y    <- pc1_sorted[1:topN] #abs(pc1_weighting)
x    <- 1:topN

# ----{每個點正規化到 [0,1]}
x_norm <- (x - min(x)) / (max(x) - min(x))
y_norm <- (y - min(y)) / (max(y) - min(y))

# ----{計算點到斜直線距離}
x1 <- 0; y1 <- 1 #左上點
x2 <- 1; y2 <- 0 #右下點
num   <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + (x2 * y1 - y2 * x1))  #分母
denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2) #分子
dist_to_line <- num / denom

# ----{elbow information}
elbow_rank  <- which.max(dist_to_line)
elbow_abs.pc1_value <- y[elbow_rank]
elbow_genes <- names(pc1_sorted)[1:elbow_rank]

cat("Elbow rank:", elbow_rank, "\n")
cat("Elbow threshold (|PC1 weighting|):", elbow_abs.pc1_value, "\n")
cat("Top genes criteria:", length(elbow_genes), "\n")

#----{Visualization}
plot.point.num <- max(1, floor(topN / 2000))  # 最多畫 ~2000 點
idx_plot   <- seq(1, topN, by = plot.point.num )

plot(idx_plot, y[idx_plot], type = "b", pch = 19, cex = 0.6,, col = "#446270",
     xlab = "Ranked genes",
     ylab = "| PC1 weighting |",
     main = sprintf("Elbow plot PC1 weighting (%d genes)", topN),
     cex.main = 1.3) #cex.main標題大小

# 標註 elbow 位置
abline(v = elbow_rank, lty = 3)
abline(h = elbow_abs.pc1_value, lty = 3)
points(elbow_rank, elbow_abs.pc1_value, pch = 21, bg = "#EF777F", col = "#446270", cex = 2.1)
text(elbow_rank, elbow_abs.pc1_value*3,
     labels = paste0("Elbow rank:", elbow_rank),
     pos = 4, cex = 1, offset = 0.5) #pos=4文字放右邊

###############【看pc1 前100個權重較高的genes富集到的pathway】###############
#取前100個PC1 weighting較高的genes
top238_genes <- names(sort(abs(pc1_weight), decreasing = TRUE))[1:238]

#----{top100 genes pathway}
library(clusterProfiler)
library(org.Hs.eg.db)

top238_genes_df<-bitr(top238_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#write.csv(top238_genes_df,file='./RNA_DATA/pc1_top238genes.csv',row.names = F)

#--------------【GO(BP)】------------------
#----{DEGs}
#進行 enrichGO 富集分析
ego_BP<-enrichGO(gene=top238_genes ,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.05, #adjusted pvalue cutoff #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)

ego_BP_df<-ego_BP@result
significant_BP_result<-ego_BP@result[ego_BP@result$p.adjust<0.05,]
nrow(significant_BP_result)#228
View(significant_BP_result)
#dotplot
dotplot(ego_BP,title='PC1 top238 genes enrichment GO(BP)',label_format = 100,showCategory=25,color='p.adjust') #color='pvalue'根據pvalue的顏色分
#label_format讓文字不會重疊

#---歸類pathway
library(clusterProfiler)
ego_BP_simple <- simplify(ego_BP, cutoff = 0.6, by = "p.adjust", select_fun = min)

# 查看簡化後的結果
View(ego_BP_simple@result)
#dotplot
dotplot(ego_BP_simple,title='PC1 top238 genes enrichment GO(BP simplify)',label_format = 100,showCategory=25,color='p.adjust')

#--------------【KEGG pathway富集分析】------------------
kk<-enrichKEGG(gene=top238_genes_df$ENTREZID,#你的基因列表
               organism='hsa', #指定物種 #hsa人類
               pvalueCutoff = 0.05 #1不進行過濾
)

View(kk@result)
significant_kk_result<-kk@result[kk@result$p.adjust<0.05,]
nrow(significant_kk_result)#16
View(significant_kk_result)

#KEGG dotplot
dotplot(kk,title='PC1 top238 genes enrichment KEGG',font.size = 12,label_format = 100,showCategory = 16
        ,color='p.adjust')


#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  ggtitle("PCA filter protein coding genes (My dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   #legend.title=element_blank(),
                   #legend.position = c(0.9,0.15),
                   legend.background=element_rect(
                    colour='gray',fill="white")
                   #legend.key.width = unit(0.6, "cm"),
                   #legend.key.height = unit(0.6, "cm"),
                   #legend.text = element_text(size = 12)
  )

###############【"External data"的 19722 Protein coding genes PCA】###############
external_prtein_coding_tpm_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene)_v2.csv")
#dim(external_prtein_coding_tpm_df) 19722 x 45

# 篩選出欄位名稱含 "NA" 或 "no"
#keep_cols <- grep("NA|no", colnames(external_prtein_coding_tpm_df), value = TRUE)
#external_prtein_coding_tpm_df<- external_prtein_coding_tpm_df[, c("gene_id", keep_cols)]
#dim(external_prtein_coding_tpm_df) 19722 x 24

#排除2個OUTLIER(no.1 & NA.15)
#external_prtein_coding_tpm_df <- external_prtein_coding_tpm_df[ ,setdiff(names(external_prtein_coding_tpm_df), c("no.1", "NA.15"))]
#dim(external_prtein_coding_tpm_df) 19722 x 43

all_gene_sample <- as.matrix(external_prtein_coding_tpm_df[, -1]) # gene x sample

rownames(all_gene_sample) <-external_prtein_coding_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重



#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  

pca.data_2[order(pca.data_2$PC1, decreasing = TRUE), ]

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
  ggtitle("PCA filter protein coding genes (External dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   #legend.title=element_blank(),
                   #legend.position = c(0.9,0.15),
                   legend.background=element_rect(
                     colour='gray',fill="white")
                   #legend.key.width = unit(0.6, "cm"),
                   #legend.key.height = unit(0.6, "cm"),
                   #legend.text = element_text(size = 12)
  )

###############【"External data"的 19720 Protein coding genes PCA (No two genes)】###############
external_prtein_coding_tpm_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene)(No two genes).csv")
#dim(external_prtein_coding_tpm_df) 19720 x 45

# 篩選出欄位名稱含 "NA" 或 "no"
#keep_cols <- grep("NA|no", colnames(external_prtein_coding_tpm_df), value = TRUE)
#external_prtein_coding_tpm_df<- external_prtein_coding_tpm_df[, c("gene_id", keep_cols)]
#dim(external_prtein_coding_tpm_df) 19722 x 24

#排除2個OUTLIER(no.1 & NA.15)
#external_prtein_coding_tpm_df <- external_prtein_coding_tpm_df[ ,setdiff(names(external_prtein_coding_tpm_df), c("no.1", "NA.15"))]
#dim(external_prtein_coding_tpm_df) 19722 x 43

all_gene_sample <- as.matrix(external_prtein_coding_tpm_df[, -1]) # gene x sample

rownames(all_gene_sample) <-external_prtein_coding_tpm_df$gene_id  # 設定欄名為基因 ID

gene_sample_transpose <- t(all_gene_sample)  # 列為 sample，欄為 gene
pca_gene_sample <- prcomp(gene_sample_transpose)
print(pca_gene_sample)


pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重



#原先數值x權重
#找pc1新的值
pc1<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc1<-0
  for(x in 1:length(pc1_weight)){
    sample_pc1<-sample_pc1+(all_gene_sample[x,i]*pc1_weight[x])
  }
  pc1<-append(pc1,sample_pc1)
}

pc1_after<-as.numeric(unlist(pc1))#as.numeric()轉成數值向量類型

#找pc2新的值
pc2<-list()
for(i in 1:ncol(all_gene_sample)){
  sample_pc2<-0
  for(x in 1:length(pc2_weight)){
    sample_pc2<-sample_pc2+(all_gene_sample[x,i]*pc2_weight[x])
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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  

pca.data_2[order(pca.data_2$PC1, decreasing = TRUE), ]

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
  ggtitle("PCA filter protein coding genes (External dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.y=element_text(size=15),
                   axis.title.x=element_text(size=15),
                   axis.title.y=element_text(size=15),
                   panel.grid.major=element_line(colour='gray'),
                   panel.grid.minor=element_line(colour='gray'),
                   #legend.title=element_blank(),
                   #legend.position = c(0.9,0.15),
                   legend.background=element_rect(
                     colour='gray',fill="white")
                   #legend.key.width = unit(0.6, "cm"),
                   #legend.key.height = unit(0.6, "cm"),
                   #legend.text = element_text(size = 12)
  )

#-----------------{DEGs:Gene expression change between NMOSD vs Control}-----------------
####{DEGs:ENSG00000198712}####
DEG_id<-"ENSG00000198712"
DEG_id<-"ENSG00000198712"
my_tpm_df<-read.csv("./RNA_DATA/My_TPM_gene_samples_matrix(intersection gene)_v2.csv")
DEG_info<-my_tpm_df[my_tpm_df$gene_id==DEG_id,]

#----轉成長df
library(tidyr)
library(dplyr) 
DEG_long <- DEG_info %>%
  pivot_longer(
    cols      = -gene_id,           # 除了 X 之外的所有欄位都 pivot
    names_to  = "Sample",     # 原本的欄位名稱會變成 Sample 欄位
    values_to = "TPM"     # 原本的值會集中到 log2CPM 欄位
  )

DEG_long$Group <- ifelse(
  grepl("^SRR", DEG_long$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
DEG_df<-as.data.frame(DEG_long)


DEG_df$Group <- factor(DEG_df$Group , levels = c("NMOSD", "Control")) #加入此>>>設定NMOSD在左邊

#畫boxplot
library(ggplot2)

My_boxplot_TPM <- ggplot(DEG_df, aes(x = Group, y =TPM, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =paste0("Expression of ", DEG_id, " (My data)") ,
       x=NULL,y = "TPM",fill="Group") +
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

print(My_boxplot_TPM)
