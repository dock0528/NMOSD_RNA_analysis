setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

###############################【DEGs:Limma+Vomm】###############################
#載入套件
#BiocManager::install(c("edgeR","limma"))
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=count_df, group=sample_df$Group) #counts:row->genes,col->samples


#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)

# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

#---【取log2cpm Matrix】
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix 
#write.csv(log2cpm_matrix,file='RNA_DATA/log2cpm_matrix_reference(no treat vs healthy).csv')

#-----------------{PCA:Batch Effect Correction"前"}-----------------

#dim(log2cpm_matrix) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_matrix[ , -1]) #-1 >>> 把第 1 欄拿掉

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
  ggtitle("PCA before batch effect correction")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.x=element_text(size=15),
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

#----------------------------------------------------------------------------------------------------
#############################【PCA "filter protein coding genes"】#############################
setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

###############################【DEGs:Limma+Vomm】###############################
#載入套件
#BiocManager::install(c("edgeR","limma"))
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

#----連到Ensembl
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----My gene matrix
All_genes<-rownames(count_df)


#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = All_genes,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #20091個
my_protein_coding_genes<-protein_coding_genes$ensembl_gene_id
#write.csv(my_protein_coding_genes,file='./RNA_DATA/My_protein_coding_genes.csv',row.names = FALSE)

#----proteing coding df
protein_coding_count_df<-count_df[rownames(count_df) %in% protein_coding_genes$ensembl_gene_id,]


########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_coding_count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=protein_coding_count_df, group=sample_df$Group) #counts:row->genes,col->samples

#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

#---【取log2cpm Matrix】
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix 
#write.csv(log2cpm_matrix,file='RNA_DATA/log2cpm_matrix(protein coding genes)).csv')

#-----------------{PCA:Batch Effect Correction"前"}-----------------

#dim(log2cpm_matrix) gene x sample: 9198 x 22('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_matrix[ , -1]) #-1 >>> 把第 1 欄拿掉

gene_sample_transpose<-t(mat) #列為sample 欄為gene
pca_gene_sample<-prcomp(gene_sample_transpose,center=TRUE, scale.=FALSE) #只去除平均，不做zscore

print(pca_gene_sample)

pca_gene_sample$rotation #權重
pc1_weight<-pca_gene_sample$rotation[,1]#pc1的權重
pc2_weight<-pca_gene_sample$rotation[,2]#pc2的權重



###############【看pc1 前100個權重較高的genes富集到的pathway】###############
#取前100個PC1 weighting較高的genes
top100_genes <- names(sort(abs(pc1_weight), decreasing = TRUE))[1:100]

#----{top100 genes pathway}
library(clusterProfiler)
library(org.Hs.eg.db)

top100_genes_df<-bitr(top100_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')

#--------------【GO(BP)】------------------
#----{DEGs}
#進行 enrichGO 富集分析
ego_BP<-enrichGO(gene=top100_genes ,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 1, #adjusted pvalue cutoff #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)

ego_BP_df<-ego_BP@result
#dotplot
dotplot(ego_BP,title='PC1 top100 genes enrichment GO(BP)',label_format = 100,showCategory=20,color='pvalue') #color='pvalue'根據pvalue的顏色分
#label_format讓文字不會重疊

#--------------【KEGG pathway富集分析】------------------
kk<-enrichKEGG(gene=top100_genes_df$ENTREZID,#你的基因列表
               organism='hsa', #指定物種 #hsa人類
               pvalueCutoff = 1 #1不進行過濾
)

View(kk@result)

#KEGG dotplot
dotplot(kk,title='PC1 top100 genes enrichment KEGG',font.size = 12,label_format = 100,showCategory = 20
        ,color='pvalue')

##################################################################################################

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
  ggtitle("PCA filter protein coding genes before batch effect correction")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.x=element_text(size=15),
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


#############################【External Dataset - PCA】#############################

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/External_RNA_raw_count_matrix(nmosd vs healthy).csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged_filter<-Raw_count_merged[,!Raw_count_merged['Antibody status',]=='MOG']
Raw_count_merged_filter <-Raw_count_merged_filter[!rownames(Raw_count_merged_filter) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged_filter) 58302 x 44

#Raw count 
count_df<-as.matrix(Raw_count_merged_filter)
storage.mode(count_df) <- "numeric"

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("NA", sample_df$Sample),  # 如果 Sample 以 "NA" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=count_df, group=sample_df$Group) #counts:row->genes,col->samples


#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降


#---【取log2cpm Matrix】
log2cpm_matrix_reference  <- v$E                           # voom 轉出的 log2-CPM matrix 


#dim(log2cpm_matrix_reference) gene x sample: 58302 x 44('gene_id+samples)

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(log2cpm_matrix_reference[ , -1]) #-1 >>> 把第 1 欄拿掉

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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
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
  ggtitle("PCA (External dataset)")+
  theme_bw()+theme(plot.title= element_text(face = "bold", size = 18, hjust = 0.5),
                   axis.text.x=element_text(size=15),
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

############################【篩"Protein Coding Genes"畫PCA】############################

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/External_RNA_raw_count_matrix(nmosd vs healthy).csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged_filter<-Raw_count_merged[,!Raw_count_merged['Antibody status',]=='MOG']
Raw_count_merged_filter <-Raw_count_merged_filter[!rownames(Raw_count_merged_filter) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged_filter) 58302 x 44

#Raw count 
count_df<-as.matrix(Raw_count_merged_filter)
storage.mode(count_df) <- "numeric"


#----連到Ensembl
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----My gene matrix
All_genes<-rownames(count_df)


#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = All_genes,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #20091個
external_protein_coding_genes<-protein_coding_genes$ensembl_gene_id
#write.csv(external_protein_coding_genes,file='./RNA_DATA/External_protein_coding_genes.csv',row.names = FALSE)

#----proteing coding df
protein_coding_count_df<-count_df[rownames(count_df) %in% protein_coding_genes$ensembl_gene_id,]

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_coding_count_df))
sample_df$Group <- ifelse(
  grepl("NA", sample_df$Sample),  # 如果 Sample 以 "NA" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=protein_coding_count_df, group=sample_df$Group) #counts:row->genes,col->samples


#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降


#---【取log2cpm Matrix】
protein_coding_log2cpm_df <- v$E                           # voom 轉出的 log2-CPM matrix 

#dim(protein_coding_log2cpm_df) 14828 x 24

#把原矩陣提出並轉成數值矩陣
mat <- as.matrix(protein_coding_log2cpm_df[ , -1]) #-1 >>> 把第 1 欄拿掉

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
  grepl("^NA", pca.data_2$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)

#顏色排序
pca.data_2$Group <- factor(pca.data_2$Group)
levels(pca.data_2$Group) #"Control" "NMOSD"  


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

############################【My & External Protein coding genes venn diagram】###################
library(VennDiagram)
my_protein_coding_genes<-read.csv('RNA_DATA/My_protein_coding_genes.csv',header=T)$x
external_protein_coding_genes<-read.csv('RNA_DATA/External_protein_coding_genes.csv',header=T)$x 
external_with_ori_protein_coding_genes<-intersect(my_protein_coding_genes,external_protein_coding_genes)#19722
#write.csv(external_with_ori_protein_coding_genes,file='./RNA_DATA/External_with_My_protein_coding_genes.csv',row.names = FALSE)
venn.list2<-list(`My dataset`=my_protein_coding_genes,`External dataset`=external_protein_coding_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='My data vs External data protein coding genes',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.5,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
